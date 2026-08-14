#!/usr/bin/env perl
# eval_combos.pl — run lasm over a grid of graph-cleaning / join feature
# combinations and score each with QUAST.
#
# Grid axes (all current experimental knobs, all default-off in lasm):
#   LASM_ADJ_JOIN          crumb adjacency joins (off|on)
#   --open-simple-loops    loop unrolling w/ multiplicity oracle (off|on)
#   LASM_UNROLL_TOL        unroll integer-fit tolerance (0.25|0.5; only
#                          meaningful with loops on)
#   --max-bubble-arm-length 0   disables bubble popping (bub on|off)
#
# Usage:
#   ./eval_combos.pl <input1 [input2 ...]> <reference.fa> <outdir> [-- extra lasm args...]
#
# Examples:
#   ./eval_combos.pl ~/Workspace/lasm-data/unmap18.bam \
#       ~/Workspace/lasm-data/chr18_chm13.fa /tmp/eval18 \
#       -- -l5 --dbg-lcmer-l 6 --olc-lcmer-l 6 --syncmer-k 11 \
#          --submer-s 3 --submer-t 7
#
#   ./eval_combos.pl ~/Workspace/lasm-data/chr22.tksm.hifi.fa \
#       ~/Workspace/lasm-data/chr22_hg002pat.fa /tmp/evalfa \
#       -- -l4 --syncmer-k 11 --submer-s 3 --submer-t 3
#
# Env knobs:
#   COMBO_ONLY=<name>   run a single combo (see names in summary)
#   JOBS=<n>            concurrent lasm runs (default 2; raise on big-RAM
#                       machines — whole-genome runs are memory-heavy)
#   LASM=<path>         lasm binary (default ./build/lasm)
#   THREADS="<flags>"   lasm thread flags (default "-t10 -@8")
#
# The summary includes lasm_s and quast_s wall times per combo, plus a
# host/input/date header line. Runs are resume-safe: re-running skips combos
# with an existing .metrics file, so an interrupted whole-genome grid
# continues where it stopped.
# Final ranked table: <outdir>/summary.txt (sorted by misassemblies asc,
# then NGA50 desc).
use strict;
use warnings;
use POSIX qw(strftime);

my $usage = "usage: eval_combos.pl '<input1 [input2 ...]>' <ref> <outdir> [-- extra lasm args...]";
sub usage_die {
    print STDERR "$usage\n";
    exit 2;
}
@ARGV >= 3 or usage_die();
my ($INPUTS, $REF, $OUT) = splice @ARGV, 0, 3;
my @extra;
if (@ARGV) {
    usage_die() unless $ARGV[0] eq '--';
    shift @ARGV;
    @extra = @ARGV;
}
my @inputs = map { s{^~/}{$ENV{HOME}/}; $_ } split ' ', $INPUTS;

my $LASM    = $ENV{LASM}    || './build/lasm';
my $JOBS    = $ENV{JOBS}    || 2;
my $THREADS = $ENV{THREADS} || '-t10 -@8';
my @threads = split ' ', $THREADS;

mkdir $OUT unless -d $OUT;

# Build the combo grid: { name, args => [], env => {} }.
my @combos;
for my $adj ('off', 'on') {
    for my $loops ('off', 'on') {
        for my $bub ('on', 'off') {
            my @tols = $loops eq 'on' ? ('025', '05') : ('025');
            for my $tol (@tols) {
                my %c = (
                    name => "adj_${adj}_loops_${loops}_bub_${bub}_tol${tol}",
                    args => [],
                    env  => {},
                );
                push @{$c{args}}, '--open-simple-loops' if $loops eq 'on';
                push @{$c{args}}, '--max-bubble-arm-length', '0' if $bub eq 'off';
                $c{env}{LASM_ADJ_JOIN} = 1 if $adj eq 'on';
                $c{env}{LASM_UNROLL_TOL} = '0.5' if $loops eq 'on' && $tol eq '05';
                push @combos, \%c;
            }
        }
    }
}

if (my $only = $ENV{COMBO_ONLY}) {
    my @keep = grep { $_->{name} eq $only } @combos;
    if (!@keep) {
        my @valid = map { $_->{name} } @combos;
        print STDERR "unknown combo: $only (valid: @valid)\n";
        exit 2;
    }
    @combos = @keep;
}

sub redirect_stdout_to {
    my ($fh) = @_;
    open my $old, '>&', \*STDOUT or die "dup stdout: $!";
    open STDOUT, '>&', $fh or die "redirect stdout: $!";
    return $old;
}

sub redirect_stderr_to {
    my ($fh) = @_;
    open my $old, '>&', \*STDERR or die "dup stderr: $!";
    open STDERR, '>&', $fh or die "redirect stderr: $!";
    return $old;
}

sub restore_fd {
    my ($old, $slot) = @_;
    open $slot, '>&', $old or die "restore fd: $!";
}

sub run_quiet {
    # system(@cmd) with stdout/stderr routed to the given files.
    # When $outf == $errf a single handle serves both (two independent
    # opens would truncate at offset 0 and interleave).
    # Returns the wait status, like system().
    my ($outf, $errf, @cmd) = @_;
    my ($out, $err);
    if ($outf eq $errf) {
        open $out, '>', $outf or die "open $outf: $!";
        $err = $out;
    } else {
        open $out, '>', $outf or die "open $outf: $!";
        open $err, '>', $errf or die "open $errf: $!";
    }
    my $oldout = redirect_stdout_to($out);
    my $olderr = redirect_stderr_to($err);
    my $rc = system(@cmd);
    restore_fd($oldout, \*STDOUT);
    restore_fd($olderr, \*STDERR);
    close $out;
    close $err if $outf ne $errf;
    return $rc;
}

# Time/RSS detection (once, in the parent), output silenced.
my $time_mode = 'none';
if (run_quiet('/dev/null', '/dev/null', '/usr/bin/time', '-v', 'true') == 0) {
    $time_mode = 'gnu';
} elsif (run_quiet('/dev/null', '/dev/null', '/usr/bin/time', '-l', 'true') == 0) {
    $time_mode = 'bsd';
}

sub rss_kb {
    my ($name) = @_;
    my $timef = "$OUT/$name.time";
    open my $fh, '<', $timef or return 0;
    local $/;
    my $txt = <$fh>;
    close $fh;
    if ($time_mode eq 'gnu' && $txt =~ /Maximum resident set size \(kbytes\): (\d+)/) {
        return $1;
    }
    if ($time_mode eq 'bsd' && $txt =~ /(\d+)\s+maximum resident set size/) {
        return int($1 / 1024);
    }
    return 0;
}

sub write_metrics {
    my ($name, @fields) = @_;
    my $tmp = "$OUT/$name.metrics.tmp$$";
    open my $fh, '>', $tmp or die "open $tmp: $!";
    print $fh join("\t", @fields), "\n";
    close $fh or die "write $tmp: $!";
    rename $tmp, "$OUT/$name.metrics" or die "rename $tmp: $!";
}

sub fail_row {
    my ($name, $lasm_s, $quast_s, $mem_gb) = @_;
    write_metrics($name, $name, 'FAIL', (0) x 8, $lasm_s, $quast_s, $mem_gb, 0, 0, 0, 0);
}

sub run_combo {
    my ($combo) = @_;
    my $name = $combo->{name};
    $| = 1;

    if (-e "$OUT/$name.metrics") {
        print "skip $name (done)\n";
        return;
    }
    print "run  $name\n";

    my $t0 = time;
    my @cmd = ($LASM, 'hifi', @inputs, '-o', "$OUT/$name.gfa", @extra, @{$combo->{args}}, @threads);
    my $rc;
    {
        # Scoped to this combo's lasm call only: a worker's later combos
        # start from the child's clean baseline %ENV again.
        local %ENV = (%ENV, %{$combo->{env}});
        if ($time_mode eq 'gnu') {
            $rc = run_quiet("$OUT/$name.log", "$OUT/$name.log", '/usr/bin/time', '-v', '-o', "$OUT/$name.time", @cmd);
        } elsif ($time_mode eq 'bsd') {
            $rc = run_quiet("$OUT/$name.log", "$OUT/$name.time", '/usr/bin/time', '-l', @cmd);
        } else {
            $rc = run_quiet("$OUT/$name.log", "$OUT/$name.log", @cmd);
        }
    }
    my $lasm_s = time - $t0;
    my $mem_gb = sprintf '%.1f', rss_kb($name) / 1048576;

    if ($rc != 0) {
        fail_row($name, $lasm_s, 0, $mem_gb);
        return;
    }
    if (!-s "$OUT/$name.gfa") {
        write_metrics($name, $name, 'EMPTY', (0) x 8, $lasm_s, 0, $mem_gb, 0, 0, 0, 0);
        return;
    }

    # GFA structure: M/C segment counts and their base totals.
    my ($m_nodes, $m_bases, $c_nodes, $c_bases) = (0, 0, 0, 0);
    open my $gfa, '<', "$OUT/$name.gfa" or die "open $OUT/$name.gfa: $!";
    while (my $line = <$gfa>) {
        chomp $line;
        next unless $line =~ /^S\t/;
        my @f = split /\t/, $line;
        my $seg = $f[1];
        my $seq = $f[2];
        if ($seg =~ /^M/) { $m_nodes++; $m_bases += length($seq); }
        elsif ($seg =~ /^C/) { $c_nodes++; $c_bases += length($seq); }
    }
    close $gfa;

    run_quiet("$OUT/$name.fa", "$OUT/$name.quast.log", 'gfatools', 'gfa2fa', "$OUT/$name.gfa");
    my $t2 = time;
    my $qrc = run_quiet("$OUT/$name.quast.log", "$OUT/$name.quast.log", 'quast.py', '-r', $REF,
        "$OUT/$name.fa", '-o', "$OUT/quast_$name", '-t', '4',
        '--min-contig', '500', '--no-plots', '--no-html', '--no-icarus');
    my $quast_s = time - $t2;

    my $tsv = "$OUT/quast_$name/report.tsv";
    if ($qrc != 0 || !-f $tsv) {
        fail_row($name, $lasm_s, $quast_s, $mem_gb);
        return;
    }
    my %v;
    open my $tsvfh, '<', $tsv or die "open $tsv: $!";
    while (my $line = <$tsvfh>) {
        chomp $line;
        my @f = split /\t/, $line;
        $v{$f[0]} = $f[1] if @f >= 2;
    }
    close $tsvfh;

    write_metrics($name, $name,
        $v{'# contigs'} // 0, $v{N50} // 0, $v{NGA50} // 0,
        $v{'# misassemblies'} // 0, $v{'# local misassemblies'} // 0,
        $v{'Genome fraction (%)'} // 0, $v{'Duplication ratio'} // 0,
        $v{'# mismatches per 100 kbp'} // 0, $v{'Largest contig'} // 0,
        $lasm_s, $quast_s, $mem_gb, $m_nodes, $m_bases, $c_nodes, $c_bases);
}

# Fork-stride parallelism: child i processes combos i, i+JOBS, ...
my $nchildren = $JOBS < @combos ? $JOBS : @combos;
for my $i (0 .. $nchildren - 1) {
    my $pid = fork;
    die "fork: $!" unless defined $pid;
    if ($pid == 0) {
        for (my $k = $i; $k < @combos; $k += $nchildren) {
            eval { run_combo($combos[$k]) };
            if ($@) { warn "child $i combo $combos[$k]{name}: $@"; exit 1; }
        }
        exit 0;
    }
}
my $bad = 0;
for my $i (1 .. $nchildren) {
    my $pid = wait;
    $bad++ if $? != 0;
}
die "some combos failed; check run logs\n" if $bad;

# Summary.
my $host = `hostname 2>/dev/null`;
chomp $host;
$host ||= 'unknown';
my $in_label;
my @in_toks = split ' ', $INPUTS;
if (@in_toks == 1) {
    ($in_label = $in_toks[0]) =~ s{.*/}{};
} else {
    my $first = $in_toks[0];
    $first =~ s{.*/}{};
    $in_label = scalar(@in_toks) . " inputs: $first...";
}
my $hdr = sprintf "# host=%s  input=%s  date=%s\n", $host, $in_label,
    strftime('%Y-%m-%dT%H:%MZ', gmtime);
my @cols = qw(name contigs N50 NGA50 misasm locals genfrac dup mism_100kbp
    largest lasm_s quast_s mem_gb M_nodes M_bases C_nodes C_bases);

my @rows;
for my $mf (glob "$OUT/*.metrics") {
    open my $fh, '<', $mf or die "open $mf: $!";
    while (my $line = <$fh>) {
        chomp $line;
        next unless $line =~ /\t/;
        push @rows, [split /\t/, $line, -1];
    }
    close $fh;
}
@rows = sort { $a->[4] <=> $b->[4] || $b->[3] <=> $a->[3] } @rows;

my @widths = map { length($_) } @cols;
for my $row (@rows) {
    for my $j (0 .. $#cols) {
        my $w = length($row->[$j] // '');
        $widths[$j] = $w if $w > $widths[$j];
    }
}
open my $sum, '>', "$OUT/summary.txt" or die "open $OUT/summary.txt: $!";
print $sum $hdr;
print $sum join('', map { sprintf '%-*s', $widths[$_] + 2, $cols[$_] } 0 .. $#cols), "\n";
for my $row (@rows) {
    print $sum join('', map { sprintf '%-*s', $widths[$_] + 2, $row->[$_] // '' } 0 .. $#cols), "\n";
}
close $sum;
open $sum, '<', "$OUT/summary.txt" or die "reopen $OUT/summary.txt: $!";
print while <$sum>;
close $sum;
