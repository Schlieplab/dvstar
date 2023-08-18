#include <gtest/gtest.h>

#include "vlmc_from_kmers/distances/containers/veb_array.hpp"
#include "vlmc_from_kmers/kmer.hpp"
#include <cstdlib>
#include <filesystem>
#include <functional>

class VebArrayTest : public ::testing::Test {
protected:
  void SetUp() override {}
};
TEST_F(VebArrayTest, VebSearchKmer) {
  std::vector<vlmc::ReadInKmer> tmp{0, 1, 2,  3,  4,  5,  6,  7,
                                    8, 9, 10, 11, 12, 13, 14, 15};
  vlmc::array::VanEmdeBoasArray arr{tmp};

  std::vector<int> expected_order{7, 14, 3, 11, 13, 15, 1, 0,
                                  2, 5,  4, 6,  9,  8,  10};

  for (int i = 1; i < 16; i++) {
    EXPECT_EQ(arr.a[i].integer_rep, expected_order[i - 1]);
  }
}
TEST_F(VebArrayTest, SearchKmer) {
  std::vector<vlmc::ReadInKmer> tmp{1, 2, 3};
  vlmc::array::VanEmdeBoasArray arr{tmp};
  auto kmer1 = vlmc::ReadInKmer(1);
  auto kmer2 = vlmc::ReadInKmer(2);
  auto kmer3 = vlmc::ReadInKmer(3);
  EXPECT_EQ(kmer1, arr.get_from_array(1));
  EXPECT_EQ(kmer2, arr.get_from_array(2));
  EXPECT_EQ(kmer3, arr.get_from_array(3));
}