# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(lasm)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Use the mold linker when available (much faster links). mold only links ELF,
# so verify the compiler actually accepts it before selecting it -- on macOS
# (Mach-O) the flag is present but unusable.
find_program(MOLD_LINKER mold)
if(MOLD_LINKER)
    include(CheckLinkerFlag)
    check_linker_flag(CXX "-fuse-ld=mold" LASM_HAVE_MOLD)
    if(LASM_HAVE_MOLD)
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
            set(CMAKE_LINKER_TYPE MOLD)
        else()
            add_link_options(-fuse-ld=mold)
        endif()
        message(STATUS "Using mold linker: ${MOLD_LINKER}")
    endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
enable_testing()
# Build options
option(LASM_NATIVE "Optimize for the build host CPU (-march=native); produces non-portable binaries" OFF)

# Compiler flags
# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wall -Wextra -Wpedantic -fsanitize=undefined,null,bounds -fsanitize-trap=all")
# Release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -ftree-vectorize -DNDEBUG -Wall -Wextra -Wpedantic")
# RelWithDebInfo flags — includes gprof instrumentation (-pg, -no-pie)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -pg -DNDEBUG -Wall -Wextra -ftree-vectorize")
set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "-pg -no-pie")
# MinSizeRel flags
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")

#We are not providing a public ABI, no need to fix this
add_compile_options($<$<CXX_COMPILER_ID:GNU>:-Wno-interference-size>)

# Optimize for the host CPU only when explicitly requested (non-portable;
# default off so distributed binaries run on other machines).
if(LASM_NATIVE)
    add_compile_options($<$<CONFIG:Release,RelWithDebInfo>:-march=native>)
    add_compile_options($<$<CONFIG:RelWithDebInfo>:-mtune=native>)
endif()

if(PGO_GENERATE)
    message(STATUS "Building with PGO instrumentation")
    add_compile_options(-fprofile-generate)
    add_link_options(-fprofile-generate)
endif()
if(PGO_USE)
    add_compile_options(-fprofile-use=${PROFILE_DATA_FILE} -fprofile-correction)
    add_link_options(-fprofile-use=${PROFILE_DATA_FILE})
endif()
# # Remove problematic flags in debug mode
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")

#     # Remove optimization that can interfere with debugging
#     string(REPLACE "-Wno-interference-size" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTS REQUIRED IMPORTED_TARGET htslib)

set(SOURCE_FILES lasm.cpp asm2.cpp lcmer.cpp dump_cores.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ./ext)

# Compile-time OLC end-skip tuning (asm2.cpp OLC_SKIP_LAST macro; default
# l/2). Plain-value macro because CMake cannot pass function-like defines on
# the compiler command line (and CMAKE_CXX_FLAGS leaks -D into the ld64 link
# step). asm2.cpp bridges LASM_OLC_SKIP_LAST_VALUE -> OLC_SKIP_LAST(l).
set(LASM_OLC_SKIP_LAST "" CACHE STRING "Override OLC_SKIP_LAST(l) (default l/2)")
if(NOT LASM_OLC_SKIP_LAST STREQUAL "")
    target_compile_definitions(${PROJECT_NAME}
                               PRIVATE LASM_OLC_SKIP_LAST_VALUE=${LASM_OLC_SKIP_LAST})
endif()

# External dependencies
include(FetchContent)


FetchContent_Declare(
    cli11_proj
    QUIET
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(cli11_proj)

FetchContent_Declare(
    gtl
    GIT_REPOSITORY https://github.com/greg7mdp/gtl.git
    GIT_TAG v1.2.0
)
FetchContent_MakeAvailable(gtl)

FetchContent_Declare(
    concurrentqueue
    GIT_REPOSITORY https://github.com/cameron314/concurrentqueue.git
    GIT_TAG v1.0.4
)
FetchContent_MakeAvailable(concurrentqueue)

FetchContent_Declare(
    libbloom
    GIT_REPOSITORY https://github.com/jvirkki/libbloom.git
    GIT_TAG v2.0
)
FetchContent_MakeAvailable(libbloom)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
                      PkgConfig::HTS gtl CLI11::CLI11
                      concurrentqueue pthread deflate z
                      $<$<PLATFORM_ID:Linux>:rt>)

# ── Unit tests (doctest) ────────────────────────────────────────────
# The TEST_CASE bodies live inline in asm2.cpp / lcmer.cpp because they
# exercise TU-internal functions. The production `lasm` binary strips them via
# DOCTEST_CONFIG_DISABLE (see compile definitions below); this target recompiles
# the same TUs together with doctest's main (tests/test_main.cpp) so the tests
# run for real, without shipping any test code in the production binary.
add_executable(lasm_tests
    tests/test_main.cpp
    asm2.cpp
    lcmer.cpp
    dump_cores.cpp
)
target_include_directories(lasm_tests PRIVATE . ./ext)
target_link_libraries(lasm_tests PRIVATE
                      PkgConfig::HTS gtl CLI11::CLI11
                      concurrentqueue pthread deflate z
                      $<$<PLATFORM_ID:Linux>:rt>)
target_compile_definitions(lasm_tests PRIVATE
    __STDC_FORMAT_MACROS=1
)
add_test(NAME unit COMMAND lasm_tests)

# Convenience target: run the assembler on a local sample (data not in repo).
add_custom_target(rtest COMMAND lasm hifi test/human2.mapped.sorted_small_set.bam test.gfa -l3 -t8)
add_dependencies(rtest ${PROJECT_NAME})

# Debug targets
add_custom_target(debug
    COMMENT "Building in debug mode"
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${PROJECT_NAME}
)


# Compile definitions
# DOCTEST_CONFIG_DISABLE strips the inline TEST_CASE bodies (in asm2.cpp /
# lcmer.cpp) from the production binary; the tests are built into lasm_tests.
target_compile_definitions(${PROJECT_NAME} PRIVATE
    __STDC_FORMAT_MACROS=1
    DOCTEST_CONFIG_DISABLE
)



# Debug-specific definitions
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG_MODE=1)

    # Enable additional debug features
    target_compile_definitions(${PROJECT_NAME} PRIVATE ENABLE_DEBUG_LOGGING=1)

endif()
