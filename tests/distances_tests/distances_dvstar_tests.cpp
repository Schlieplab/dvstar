#pragma once 

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

//Our files 
#include "vlmc_container.hpp"
#include "distances/dvstar.hpp"
#include "vlmc_from_kmers/dvstar.hpp"
#include "global_aliases.hpp"

//Original implementation files
#include <kmc_file.h>
#include "vlmc_from_kmers/dvstar.hpp"
#include "vlmc_from_kmers/build_vlmc.hpp"
#include "vlmc_from_kmers/kmer.hpp"

using VLMC_vector = container::VLMC_vector;

extern const out_t error_tolerance;

class DvstarTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_fasta{"NC_028367.1.fa"};
  std::filesystem::path second_fasta{"NC_045512.2.fa"};
  std::filesystem::path third_fasta{"NC_001497.2.fa"};

  std::filesystem::path first_bintree{"../data/test_VLMCs/sequences_1.bintree"};
  std::filesystem::path second_bintree{"../data/test_VLMCs/sequences_2.bintree"};
  std::filesystem::path third_bintree{"../data/test_VLMCs/sequences_3.bintree"};

  int background_order = 0;

  std::function<out_t(VLMC_vector &, VLMC_vector &)> dist_func = [&](auto &left, auto &right) {
      return distance::dvstar<VLMC_vector>(left, right);
  };
};

// Vector Tests
TEST_F(DvstarTests, Identity_vector) {
  VLMC_vector first_vlmc{first_bintree};

  out_t dist_vector = dist_func(first_vlmc, first_vlmc);
  EXPECT_NEAR(0.0, dist_vector, error_tolerance);
}

TEST_F(DvstarTests, EqualDistance_vector) {
  VLMC_vector first_vlmc{first_bintree};
  VLMC_vector second_vlmc{second_bintree};

  out_t dist_vector = dist_func(first_vlmc, second_vlmc);
  out_t old_dvstar_implementation = vlmc::dvstar(first_bintree, second_bintree, background_order);

  EXPECT_NEAR(old_dvstar_implementation, dist_vector, error_tolerance);
}

TEST_F(DvstarTests, Symmetry_vector) {
  VLMC_vector first_vlmc{first_bintree};
  VLMC_vector second_vlmc{second_bintree};

  out_t dist_one_vector = dist_func(first_vlmc, second_vlmc);
  out_t dist_two_vector = dist_func(second_vlmc, first_vlmc);

  EXPECT_NEAR(dist_one_vector, dist_two_vector, error_tolerance);
}

TEST_F(DvstarTests, multiple_runs_vector) {
  VLMC_vector first_vlmc{first_bintree};
  VLMC_vector second_vlmc{second_bintree};
  size_t runs = 10; 
  out_t prev = dist_func(first_vlmc, second_vlmc);
  for (int i = 0; i < runs; i++){
    out_t current = dist_func(first_vlmc, second_vlmc);
    EXPECT_NEAR(prev, current, error_tolerance); 
    prev = current; 
  }
}

TEST_F(DvstarTests, TestBackgroundOrder) {
  for (int order = 0; order < 5; order++){
    container::VLMC_vector first_vlmc_vector{first_bintree, order};
    container::VLMC_vector second_vlmc_vector{second_bintree, order};

    auto dist_vector = distance::dvstar<VLMC_vector>(first_vlmc_vector, second_vlmc_vector);
    auto old_dvstar_implementation = vlmc::dvstar(first_bintree, second_bintree, order);

    EXPECT_NEAR(old_dvstar_implementation, dist_vector, error_tolerance);
  }
}
