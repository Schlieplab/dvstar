cmake_minimum_required(VERSION 3.4...3.12)

project(kmc_vlmc_tests CXX)

include(CTest)
find_package(Threads REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BENCHMARK_ENABLE_GTEST_TESTS OFF)

macro(dependency name)
    # Download and unpack googletest at configure time
    configure_file(CMakeLists.txt.${name} ${name}-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name}-download)
    if (result)
        message(FATAL_ERROR "CMake step for ${name} failed: ${result}")
    endif ()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${name}-download)
    if (result)
        message(FATAL_ERROR "Build step for ${name} failed: ${result}")
    endif ()

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/${name}-src
            ${CMAKE_CURRENT_BINARY_DIR}/${name}-build
            EXCLUDE_FROM_ALL)
endmacro(dependency)

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

dependency(googletest)
dependency(benchmark)

if (NOT KMC_Libraries)
    add_library(kmc_file SHARED ../KMC/kmc_api/kmc_file.cpp)
    add_library(kmer_api SHARED ../KMC/kmc_api/kmer_api.cpp)
    add_library(mmer SHARED ../KMC/kmc_api/mmer.cpp)
endif ()

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif ()

macro(test_case unit_test)
    get_filename_component(target "${unit_test}" NAME_WE)

    add_executable(${target} ${unit_test})
    target_link_libraries(${target} gtest_main)
    target_link_libraries(${target} Threads::Threads)
    target_link_libraries(${target} ${CountVLMC_LIBRARIES})

    add_test(NAME ${target} COMMAND ${target})

    unset(unit_test)
    unset(target)
endmacro(test_case)

macro(benchmark_case case)
    get_filename_component(target "${case}" NAME_WE)

    add_executable(${target} ${case})
    target_link_libraries(${target} benchmark)
    target_link_libraries(${target} Threads::Threads)
    target_link_libraries(${target} ${CountVLMC_LIBRARIES})

    unset(case)
    unset(target)
endmacro(benchmark_case)

test_case(kmer_tests.cpp)
test_case(support_pruning_tests.cpp)
test_case(similarity_pruning_tests.cpp)
test_case(count_correctness_tests.cpp)
test_case(regression_tests.cpp)
test_case(negative_log_likelihood_tests.cpp)
test_case(method_regression.cpp)

benchmark_case(benchmarks.cpp)
benchmark_case(method_benchmarks.cpp)

add_subdirectory(distances_tests)