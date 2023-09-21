cmake_minimum_required(VERSION 3.4...3.12)

project(CountVLMC CXX)

option(CountVLMC_COMPILE_TESTS "To compile tests, specify -DCountVLMC_COMPILE_TESTS=ON")

set(MAIN_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(MAIN_PROJECT ON)
endif ()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "/usr/share/cmake/Modules/")

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED CXX)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams)

add_library(kmc_file STATIC submodules/KMC/kmc_api/kmc_file.cpp)
add_library(kmer_api STATIC submodules/KMC/kmc_api/kmer_api.cpp)
add_library(mmer STATIC submodules/KMC/kmc_api/mmer.cpp)
set(KMC_Libraries TRUE)

# STXXL
set(USE_BOOST OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/stxxl)
# apply STXXL CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# Parallelization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ltbb -fopenmp -fpermissive")

# CLI11 for command line interface
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/CLI11)

set(SKIP_PERFORMANCE_COMPARISON ON)

# Add include paths
set(INCLUDE_DIRS ${STXXL_INCLUDE_DIRS} "${PROJECT_SOURCE_DIR}/submodules/stxxl/include" "${PROJECT_SOURCE_DIR}/submodules/robin-hood-hashing/src/include" "${PROJECT_SOURCE_DIR}/include" "${CMAKE_CURRENT_SOURCE_DIR}/submodules/KMC/kmc_api" "${PROJECT_SOURCE_DIR}/submodules/cereal/include" "${PROJECT_SOURCE_DIR}/submodules/eigen" "${PROJECT_SOURCE_DIR}/submodules/unordered_dense/include")

include_directories("${INCLUDE_DIRS}")

set(USE_BOOST OFF)
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_EIGEN ON)
set(HIGHFIVE_EXAMPLES FALSE)

find_package(HDF5 REQUIRED COMPONENTS CXX)
add_subdirectory("${CMAKE_SOURCE_DIR}/submodules/HighFive")

if (NOT MAIN_PROJECT)
    set(CountVLMC_INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)
endif ()


set(CountVLMC_LIBRARIES ${STXXL_LIBRARIES} tbb kmc_file kmer_api mmer CLI11::CLI11 HighFive Boost::iostreams)

add_library(CountVLMC INTERFACE)
add_library(CountVLMC::CountVLMC ALIAS CountVLMC)

target_link_libraries(CountVLMC INTERFACE ${CountVLMC_LIBRARIES})
target_include_directories(CountVLMC INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/KMC/kmc_api/>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/>
        $<INSTALL_INTERFACE:include/>)

if (MAIN_PROJECT)
    add_executable(dvstar src/build_vlmc.cpp)
    target_link_libraries(dvstar Threads::Threads)
    target_link_libraries(dvstar ${CountVLMC_LIBRARIES})

    if (CountVLMC_COMPILE_TESTS)
        add_subdirectory(tests)
    endif ()
endif ()
