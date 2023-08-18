#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "vlmc_from_kmers/distances/calc_dists.hpp"
#include "vlmc_from_kmers/distances/get_cluster.hpp"
#include "vlmc_from_kmers/distances/global_aliases.hpp"
#include "vlmc_from_kmers/distances/vlmc_container.hpp"

extern const out_t error_tolerance;

class ParallelTest : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_directory{"."};
};

TEST_F(ParallelTest, SequentialEqParallel) {
  auto cluster =
      vlmc::get_cluster<vlmc::container::SortedVector>(first_directory, 1, 0);

  matrix_t distances_parallel{cluster.size(), cluster.size()};
  matrix_t distances_sequential{cluster.size(), cluster.size()};

  auto fun_parallel = [&](size_t start_index, size_t stop_index) {
    vlmc::calc_dist::calculate_full_slice<vlmc::container::SortedVector>(
        start_index, stop_index, 0, cluster.size(), distances_parallel, cluster,
        cluster);
  };

  auto fun_sequential = [&](size_t start_index, size_t stop_index) {
    vlmc::calc_dist::calculate_full_slice<vlmc::container::SortedVector>(
        start_index, stop_index, 0, cluster.size(), distances_sequential,
        cluster, cluster);
  };

  vlmc::parallel::parallelize(cluster.size(), fun_parallel, 2);

  fun_sequential(0, cluster.size());

  EXPECT_TRUE(distances_parallel.isApprox(distances_sequential));
}
