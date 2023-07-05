#pragma once
#include <gtest/gtest.h>

#include "vlmc_from_kmers/distances/containers/b_tree_array.hpp"

class BTreeTest : public ::testing::Test {
protected:
  void SetUp() override {}
};
TEST_F(BTreeTest, BTreeSearchKmer) {
  std::vector<vlmc::ReadInKmer> tmp{0, 1, 2,  3,  4,  5,  6,  7,
                                    8, 9, 10, 11, 12, 13, 14, 15};
  vlmc::array::B_Tree arr{tmp};

  std::vector expected_order{13, 2, 5, 11, 12, 14, 15, 0, 1, 3, 4, 6, 7, 9, 10};

  for (int i = 1; i < 16; i++) {
    EXPECT_EQ(arr.a[i].integer_rep, expected_order[i - 1]);
  }
}
TEST_F(BTreeTest, SearchKmer) {
  std::vector<vlmc::ReadInKmer> tmp{1, 2, 3};
  vlmc::array::B_Tree arr{tmp};
  auto kmer1 = vlmc::ReadInKmer(1);
  auto kmer2 = vlmc::ReadInKmer(2);
  auto kmer3 = vlmc::ReadInKmer(3);
  EXPECT_EQ(kmer1, arr.get_from_array(1));
  EXPECT_EQ(kmer2, arr.get_from_array(2));
  EXPECT_EQ(kmer3, arr.get_from_array(3));
}