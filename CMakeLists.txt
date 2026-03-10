# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(lasm)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_USING_LINKER_mold)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
enable_testing()
# Compiler flags
# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -O0 -DDEBUG -Wall -Wextra -Wpedantic -fsanitize=undefined,null,bounds -fsanitize-trap=all")
# Release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -pg -march=native -ftree-vectorize -DNDEBUG -Wall -Wextra -Wpedantic")
# RelWithDebInfo flags
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG -Wall -Wextra -mtune=native -ftree-vectorize -march=native")
# MinSizeRel flags
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")

if(PGO_GENERATE)
    message(STATUS "Building with PGO instrumentation")
    add_compile_options(-fprofile-generate)
    add_link_options(-fprofile-generate)
endif()
if(PGO_USE)
    add_compile_options(-fprofile-use=${PROFILE_DATA_FILE} -fprofile-correction)
    add_link_options(-fprofile-use=${PROFILE_DATA_FILE})
endif()
# Set linker flags for release builds
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie")
endif()

# # Remove problematic flags in debug mode
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     # Remove optimization that can interfere with debugging
#     string(REPLACE "-Wno-interference-size" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTS REQUIRED htslib)
include_directories(${HTS_INCLUDE_DIRS})
link_directories(${HTS_LIBRARY_DIRS})
add_definitions(${HTS_CFLAGS_OTHER})

set(SOURCE_FILES lasm.cpp asm2.cpp lcmer.cpp map.cpp)
include_directories(./ext)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

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

FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(eigen)


# Link libraries
target_link_libraries(${PROJECT_NAME} ${HTS_LIBRARIES} gtl CLI11::CLI11 concurrentqueue pthread deflate z Eigen3::Eigen)

# Unit tests (doctest)
# Lightweight DBG builder test tool (single-thread, no bloom/locks).
#
# IMPORTANT:
# - Rely on stb-style header-only lcptools implementation (included by the tool),
#   so do not link ext/lcptools/*.c here.
# - Do NOT link lcmer.cpp here, because it includes doctest-based tests and
#   will pull doctest runtime symbols at link time.
# Temporarily disabled - source file missing
# add_executable(light_dbg_build
#     tools/light_dbg_build.cpp
# )
# target_include_directories(light_dbg_build PRIVATE . ./ext)
# target_link_libraries(light_dbg_build PRIVATE gtl CLI11::CLI11)
# target_compile_definitions(light_dbg_build PRIVATE
#     SYNCMER_RULE=DensitySyncmerRule
# )

# Core dump tool
add_executable(dump_cores dump_cores.cpp lcmer.cpp)
target_include_directories(dump_cores PRIVATE . ./ext)
target_link_libraries(dump_cores PRIVATE ${HTS_LIBRARIES} gtl CLI11::CLI11 pthread deflate z)
target_compile_definitions(dump_cores PRIVATE SYNCMER_RULE=DensitySyncmerRule DOCTEST_CONFIG_DISABLE)

# Tests and custom targets
add_test(NAME test COMMAND lasm hifi ../tests/human2.mapped.sorted_small_set.bam test.gfa -l3 -t8 --density 0.05)
add_custom_target(rt COMMAND lasm hifi test/human2.mapped.sorted_small_set.bam test.gfa -l3 -t8)
add_dependencies(rt ${PROJECT_NAME})

# Debug targets
add_custom_target(debug
    COMMENT "Building in debug mode"
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${PROJECT_NAME}
)


# Compile definitions
target_compile_definitions(${PROJECT_NAME} PUBLIC
    SYNCMER_RULE=DensitySyncmerRule
)


# Debug-specific definitions
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PUBLIC DEBUG_MODE=1)

    # Enable additional debug features
    target_compile_definitions(${PROJECT_NAME} PUBLIC ENABLE_DEBUG_LOGGING=1)

endif()
