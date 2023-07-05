#include <gtest/gtest.h>

#include "vlmc_from_kmers/distances/containers/eytzinger_array.hpp"
#include "vlmc_from_kmers/kmer.hpp"
#include <cstdlib>
#include <filesystem>
#include <functional>

class EytzArrayTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TEST_F(EytzArrayTest, SearchKmer) {
  std::vector<vlmc::ReadInKmer> tmp{1, 2, 3};
  vlmc::array::EytzingerArray arr{tmp};
  auto kmer1 = vlmc::ReadInKmer(1);
  auto kmer2 = vlmc::ReadInKmer(2);
  auto kmer3 = vlmc::ReadInKmer(3);
  EXPECT_EQ(kmer1, arr.get_from_array(1));
  EXPECT_EQ(kmer2, arr.get_from_array(2));
  EXPECT_EQ(kmer3, arr.get_from_array(3));
}
//TEST_F(EytzArrayTest, Construct) {
//  std::vector<vlmc::ReadInKmer> tmp{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//  vlmc::array::EytzingerArray arr(tmp);
//  for (int i = 1; i < 16; i++) {
//    std::cout << arr.b[i].integer_rep << "\n";
//  }
//}