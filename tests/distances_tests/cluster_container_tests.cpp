#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#include "vlmc_from_kmers/distances/vlmc_container.hpp"
#include "vlmc_from_kmers/distances/cluster_container.hpp"

using vlmc_c = vlmc::container::SortedVector;
using cluster_c = vlmc::container::ClusterContainer<vlmc_c>;

class ClusterContainerTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_bintree{"NC_001497.bintree"};
  std::filesystem::path second_bintree{"NC_045512.2.bintree"};
  std::filesystem::path third_bintree{"NC_001497.2.bintree"};

  vlmc_c first_vlmc{first_bintree};
  vlmc_c second_vlmc{second_bintree};
  vlmc_c third_vlmc{third_bintree};
};


TEST_F(ClusterContainerTest, AddToContainer) {
  cluster_c container {};
  container.push(first_vlmc);
  container.push(second_vlmc);
  container.push(third_vlmc);
  EXPECT_EQ(container.size(), 3);
}

TEST_F(ClusterContainerTest, VlmcSizeNonZeroAfterAddToContainer) {
  cluster_c container {};
  container.push(first_vlmc);
  EXPECT_GT(container.get(0).size(), 0);
}

// TEST_F(ClusterContainerTest, KmerContainerGet){
//   container::KmerCluster container{};
//   container::KmerPair kmer0 = container::KmerPair(container::ReadInKmer(0), 0);
//   container.push(kmer0);
//   std::vector<container::KmerPair> kmer_bucket = container.get(container.get_bucket(kmer0));
//   EXPECT_EQ(kmer0.id, kmer_bucket[0].id);
// }