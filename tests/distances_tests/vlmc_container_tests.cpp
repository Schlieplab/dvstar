#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <functional>

#include "../read_helper.hpp"
#include "vlmc_from_kmers/distances/read_in_kmer.hpp"
#include "vlmc_from_kmers/distances/vlmc_container.hpp"

class VlmcContainerTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TEST_F(VlmcContainerTest, IndexRep1) {
  vlmc::container::SortedVector container{};
  std::string kmer_string{"A"};
  auto created = create_kmer(kmer_string);
  vlmc::ReadInKmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 1);
}

TEST_F(VlmcContainerTest, IndexRep2) {
  vlmc::container::SortedVector container{};
  std::string kmer_string{"T"};
  auto created = create_kmer(kmer_string);
  vlmc::ReadInKmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 4);
}
TEST_F(VlmcContainerTest, IndexRep3) {
  vlmc::container::SortedVector container{};
  std::string kmer_string{"AA"};
  auto created = create_kmer(kmer_string);
  vlmc::ReadInKmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 5);
}
TEST_F(VlmcContainerTest, IndexRep4) {
  vlmc::container::SortedVector container{};
  std::string kmer_string{"AT"};
  auto created = create_kmer(kmer_string);
  vlmc::ReadInKmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 8);
}
TEST_F(VlmcContainerTest, IndexRep5) {
  vlmc::container::SortedVector container{};
  std::string kmer_string{"CC"};
  auto created = create_kmer(kmer_string);
  vlmc::ReadInKmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 10);
}
TEST_F(VlmcContainerTest, IndexRep6) {
  vlmc::container::SortedVector container{};
  std::string kmer_string{"TG"};
  auto created = create_kmer(kmer_string);
  vlmc::ReadInKmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 19);
}
TEST_F(VlmcContainerTest, IndexRep7) {
  vlmc::container::SortedVector container{};
  std::string kmer_string{"AAA"};
  auto created = create_kmer(kmer_string);
  vlmc::ReadInKmer new_kmer{created};
  EXPECT_EQ(new_kmer.integer_rep, 21);
}