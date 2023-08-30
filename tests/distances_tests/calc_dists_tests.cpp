#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <string>

#include "vlmc_from_kmers/distances/calc_dists.hpp"
#include "vlmc_from_kmers/distances/cluster_container.hpp"
#include "vlmc_from_kmers/distances/distances/dvstar.hpp"
#include "vlmc_from_kmers/distances/get_cluster.hpp"
#include "vlmc_from_kmers/distances/global_aliases.hpp"
#include "vlmc_from_kmers/distances/vlmc_container.hpp"

// Original implementation files
#include "vlmc_from_kmers/dvstar.hpp"

using cluster_c =
    vlmc::container::ClusterContainer<vlmc::container::SortedVector>;
const out_t error_tolerance = 0.00000000000001;

class CalcDistsTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_bintree{"NC_028367.bintree"};
  std::filesystem::path second_bintree{"NC_045512.bintree"};
  std::filesystem::path third_bintree{"NC_001497.bintree"};

  std::filesystem::path path_to_bintrees{"."};

  vlmc::container::SortedVector first_vlmc{first_bintree};
  vlmc::container::SortedVector second_vlmc{second_bintree};
  vlmc::container::SortedVector third_vlmc{third_bintree};

  size_t background_order = 0;
};

cluster_c create_cluster(vlmc::container::SortedVector &vlmc, size_t size) {
  cluster_c cluster{};
  for (size_t i = 0; i < size; i++) {
    cluster.push(vlmc);
  }
  return cluster;
}

// Tests when comparing two directories
TEST_F(CalcDistsTests, SizeTwoDir) {
  for (size_t x = 1; x < 5; x++) {
    for (size_t y = 1; y < 5; y++) {
      cluster_c cluster_left = create_cluster(first_vlmc, x);
      cluster_c cluster_right = create_cluster(second_vlmc, y);

      matrix_t distances =
          vlmc::calc_dist::calculate_distances<vlmc::container::SortedVector>(
              cluster_left, cluster_right, 1);
      EXPECT_EQ(distances.size(), x * y);
      EXPECT_EQ(distances.rows(), x);
      EXPECT_EQ(distances.cols(), y);
    }
  }
}

TEST_F(CalcDistsTests, AllValsTwoDir) {
  for (size_t x = 1; x < 2; x++) {
    cluster_c cluster_left = create_cluster(first_vlmc, x);
    cluster_c cluster_right = create_cluster(second_vlmc, x);

    matrix_t distances =
        vlmc::calc_dist::calculate_distances<vlmc::container::SortedVector>(
            cluster_left, cluster_right, 1);

    for (size_t row = 0; row < distances.rows(); row++) {
      for (size_t col = 0; col < distances.cols(); col++) {
        EXPECT_GT(distances(row, col), 0);
      }
    }
  }
}

TEST_F(CalcDistsTests, ValueCheckTwoDir) {
  // Vector Implementation
  auto [left_cluster_v, _lids] =
      vlmc::get_cluster<vlmc::container::SortedVector>(path_to_bintrees, 1,
                                                       background_order);
  auto [right_cluster_v, _rids] =
      vlmc::get_cluster<vlmc::container::SortedVector>(path_to_bintrees, 1,
                                                       background_order);
  matrix_t distances_vector =
      vlmc::calc_dist::calculate_distances<vlmc::container::SortedVector>(
          left_cluster_v, right_cluster_v, 1);

  // Sorted Vector Implementation
  auto [left_cluster_s, _lsids] =
      vlmc::get_cluster<vlmc::container::SortedVector>(path_to_bintrees, 1,
                                                       background_order);
  auto [right_cluster_s, _rsids] =
      vlmc::get_cluster<vlmc::container::SortedVector>(path_to_bintrees, 1,
                                                       background_order);
  matrix_t distances_sorted_vector =
      vlmc::calc_dist::calculate_distances<vlmc::container::SortedVector>(
          left_cluster_s, right_cluster_s, 1);

  // B-tree Implementation
  auto [left_cluster_b, _lcids] = vlmc::get_cluster<vlmc::container::BTree>(
      path_to_bintrees, 1, background_order);
  auto [right_cluster_b, _] = vlmc::get_cluster<vlmc::container::BTree>(
      path_to_bintrees, 1, background_order);
  matrix_t distances_b_tree =
      vlmc::calc_dist::calculate_distances<vlmc::container::BTree>(
          left_cluster_b, right_cluster_b, 1);

  // HashMap Implementation
  auto [left_cluster_h, _lhids] = vlmc::get_cluster<vlmc::container::HashMap>(
      path_to_bintrees, 1, background_order);
  auto [right_cluster_h, rhids] = vlmc::get_cluster<vlmc::container::HashMap>(
      path_to_bintrees, 1, background_order);
  matrix_t distances_hashmap =
      vlmc::calc_dist::calculate_distances<vlmc::container::HashMap>(
          left_cluster_h, right_cluster_h, 1);

  // Kmer major implementation
  auto left_cluster_k =
      vlmc::get_kmer_cluster(path_to_bintrees, 1, background_order);
  auto right_cluster_k =
      vlmc::get_kmer_cluster(path_to_bintrees, 1, background_order);
  matrix_t distances_k_major = vlmc::calc_dist::calculate_distance_major(
      left_cluster_k, right_cluster_k, 1);

  // Veb Implementation
  auto [left_cluster_veb, _lvebids] =
      vlmc::get_cluster<vlmc::container::VanEmdeBoasTree>(path_to_bintrees, 1,
                                                          background_order);
  auto [right_cluster_veb, _rvebids] =
      vlmc::get_cluster<vlmc::container::VanEmdeBoasTree>(path_to_bintrees, 1,
                                                          background_order);
  matrix_t distances_veb =
      vlmc::calc_dist::calculate_distances<vlmc::container::VanEmdeBoasTree>(
          left_cluster_veb, right_cluster_veb, 1);

  // Eytzinger Array
  auto [left_cluster_ey, _leyids] =
      vlmc::get_cluster<vlmc::container::EytzingerTree>(path_to_bintrees, 1,
                                                        background_order);
  auto [right_cluster_ey, _reyids] =
      vlmc::get_cluster<vlmc::container::EytzingerTree>(path_to_bintrees, 1,
                                                        background_order);
  matrix_t distances_ey =
      vlmc::calc_dist::calculate_distances<vlmc::container::EytzingerTree>(
          left_cluster_ey, right_cluster_ey, 1);

  // Sorted Search
  auto [left_cluster_sortsearch, _lssids] =
      vlmc::get_cluster<vlmc::container::SortedSearch>(path_to_bintrees, 1,
                                                       background_order);
  auto [right_cluster_sortsearch, _rssids] =
      vlmc::get_cluster<vlmc::container::SortedSearch>(path_to_bintrees, 1,
                                                       background_order);
  matrix_t distances_sortsearch =
      vlmc::calc_dist::calculate_distances<vlmc::container::SortedSearch>(
          left_cluster_sortsearch, right_cluster_sortsearch, 1);

  // Dvstar Original implementation
  matrix_t distances_org_dvstar{distances_vector.cols(),
                                distances_vector.rows()};
  int x = 0;
  for (const auto &dir_entry_x :
       recursive_directory_iterator(path_to_bintrees)) {
    if (dir_entry_x.path().extension() != ".bintree") {
      continue;
    }
    int y = 0;
    for (const auto &dir_entry_y :
         recursive_directory_iterator(path_to_bintrees)) {
      if (dir_entry_y.path().extension() != ".bintree") {
        continue;
      }
      distances_org_dvstar(x, y) =
          vlmc::dvstar(dir_entry_x, dir_entry_y, background_order);
      y++;
    }
    x++;
  }
  for (int x = 0; x < distances_vector.cols(); x++) {
    for (int y = 0; y < distances_vector.rows(); y++) {
      if (x == y) {
        EXPECT_NEAR(0.0, distances_vector(x, y), error_tolerance);
//        EXPECT_NEAR(0.0, distances_org_dvstar(x, y), error_tolerance);
      } else {
        EXPECT_NEAR(distances_org_dvstar(x, y), distances_vector(x, y),
                    error_tolerance) << x << " " << y;
        EXPECT_NEAR(distances_vector(x, y), distances_sorted_vector(x, y),
                    error_tolerance) << x << " " << y;
        EXPECT_NEAR(distances_vector(x, y), distances_b_tree(x, y),
                    error_tolerance) << x << " " << y;
        EXPECT_NEAR(distances_vector(x, y), distances_hashmap(x, y),
                    error_tolerance) << x << " " << y;
        EXPECT_NEAR(distances_vector(x, y), distances_veb(x, y),
                    error_tolerance) << x << " " << y;
        //        EXPECT_NEAR(distances_vector(x, y), distances_k_major(x, y),
        //                    error_tolerance) << x << " " << y;
        EXPECT_NEAR(distances_vector(x, y), distances_ey(x, y),
                    error_tolerance) << x << " " << y;
        EXPECT_NEAR(distances_vector(x, y), distances_sortsearch(x, y),
                    error_tolerance) << x << " " << y;
      }
    }
  }
}

TEST_F(CalcDistsTests, ValueCheckOneDir) {
  // Vector Implementation
  auto [left_cluster_v, _lvids] =
      vlmc::get_cluster<vlmc::container::SortedVector>(path_to_bintrees, 1,
                                                       background_order);
  auto [right_cluster_v, _rvids] =
      vlmc::get_cluster<vlmc::container::SortedVector>(path_to_bintrees, 1,
                                                       background_order);
  matrix_t distances_vector_two_dirs =
      vlmc::calc_dist::calculate_distances<vlmc::container::SortedVector>(
          left_cluster_v, right_cluster_v, 1);

  auto [cluster_v, _] = vlmc::get_cluster<vlmc::container::SortedVector>(
      path_to_bintrees, 1, background_order);

  for (int nr_cores = 1; nr_cores < 9; nr_cores++) {
    matrix_t distances_vector =
        vlmc::calc_dist::calculate_distances<vlmc::container::SortedVector>(
            cluster_v, nr_cores);
    for (int x = 0; x < distances_vector.cols(); x++) {
      for (int y = 0; y < distances_vector.rows(); y++) {
        if (x <= y) {
          EXPECT_NEAR(distances_vector_two_dirs(x, y), distances_vector(x, y),
                      error_tolerance);
        } else {
          EXPECT_NEAR(0.0, distances_vector(x, y), error_tolerance);
        }
      }
    }
  }
}

// Tests for inter-directory comparisons
TEST_F(CalcDistsTests, SizeOneDir) {
  for (size_t x = 1; x < 5; x++) {
    cluster_c cluster_left = create_cluster(first_vlmc, x);
    matrix_t distances = vlmc::calc_dist::calculate_distances(cluster_left, 1);
    EXPECT_EQ(distances.size(), x * x);
    EXPECT_EQ(distances.rows(), x);
    EXPECT_EQ(distances.cols(), x);
  }
}

TEST_F(CalcDistsTests, AllValsOneDir) {
  for (size_t x = 1; x < 2; x++) {
    cluster_c cluster = create_cluster(first_vlmc, x);

    matrix_t distances = vlmc::calc_dist::calculate_distances(cluster, 1);

    for (size_t row = 0; row < distances.rows(); row++) {
      for (size_t col = 1 + row; col < distances.cols(); col++) {
        EXPECT_GT(distances(row, col), 0);
      }
    }
  }
}

TEST_F(CalcDistsTests, RegressionTest) {
  std::filesystem::path first_bintree_path{"NC_028367.bintree"};
  std::filesystem::path second_bintree_path{"NC_045512.bintree"};

  auto left = vlmc::container::SortedSearch{first_bintree_path, 0};
  auto right = vlmc::container::SortedSearch{second_bintree_path, 0};
  auto sorted_vector_result =
      vlmc::distance::dvstar<vlmc::container::SortedSearch>(left, right);

  auto left_kmers = vlmc::get_sorted_kmers(first_bintree_path);
  auto right_kmers = vlmc::get_sorted_kmers(second_bintree_path);

  auto dvstar_result = dvstar(left_kmers, right_kmers, 0);

  EXPECT_NEAR(sorted_vector_result, dvstar_result, error_tolerance);
}

TEST_F(CalcDistsTests, IdentityShouldBeZeroTest) {
  std::filesystem::path first_bintree_path{"EU579861.1.fa.bintree"};
  auto left = vlmc::container::SortedSearch{first_bintree_path, 0};

  auto sorted_vector_result =
      vlmc::distance::dvstar<vlmc::container::SortedSearch>(left, left);

  std::cout << sorted_vector_result;

  EXPECT_NEAR(sorted_vector_result, 0, error_tolerance);
}