#ifndef _LAYOUT_GRAPH_H
#define _LAYOUT_GRAPH_H

#include <gtl/vector.hpp>
#include <gtl/phmap_fwd_decl.hpp>
#include "graph.h"

using layout_small_graph = gtl::flat_hash_map<uint64_t, dbg_index::bidbg_val>;

struct lsg_cleaning{
    void cut_tips(layout_small_graph &lsg, gtl::vector<u64> &to_fix){

    }
    void cut_bubbles(layout_small_graph &lsg){

    }
    void simplify_simple_paths(layout_small_graph &lsg){

    }
    void fix_edges(layout_small_graph &lsg, gtl::vector<u64> &to_fix){

    }
    void run(layout_small_graph &lsg){

    }
};

#endif
