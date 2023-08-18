#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

// Our files
#include "vlmc_from_kmers/distances/get_cluster.hpp"
#include "vlmc_from_kmers/distances/vlmc_container.hpp"

class GetClusterTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path path_to_bintrees{"."};
  std::filesystem::path path_to_vlmc{"NC_001497.bintree"};
};

TEST_F(GetClusterTest, ClusterGetWithVlmcVector) {
  auto container =
      vlmc::get_cluster<vlmc::container::SortedVector>(path_to_bintrees, 1, 0);

  // std::cout << container.get(0).get(0).to_string() << std::endl;
  EXPECT_GT(container.size(), 0);
  EXPECT_GT(container.get(0).size(), 0);
}
