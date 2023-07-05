#pragma once

#include <Eigen/Core>
#include <mutex>

#include "cluster_container.hpp"
#include "vlmc_container.hpp"
#include "parallel.hpp"
#include "distances/dvstar.hpp"
#include "global_aliases.hpp"
#include "utils.hpp"

namespace calc_dist {
  using kmer_pair = cluster_container::Kmer_Pair;

  template <typename VC>
  void calculate_triangle_slice(
    int x1, int y1, int x2, int y2, int x3, int y3, matrix_t& distances,
    cluster_container::Cluster_Container<VC>& cluster_left, cluster_container::Cluster_Container<VC>& cluster_right) {

    auto rec_fun = [&](int left, int right) {
      if (distances(left, right) == 0) {
        distances(left, right) = distance::dvstar<VC>(cluster_left.get(left), cluster_right.get(right));
      }
    };

    if (x1 < x2) {
      if (y2 > y3) {
        utils::triangle_recursion(x1, x3, y1, y2, x1, y1, x2, y2, x3, y3, rec_fun);
      }
      else {
        utils::triangle_recursion(x1, x3, y1, y3, x1, y1, x2, y2, x3, y3, rec_fun);
      }
    }
    else {
      if (y2 > y3) {
        utils::triangle_recursion(x2, x3, y1, y2, x1, y1, x2, y2, x3, y3, rec_fun);
      }
      else {
        utils::triangle_recursion(x2, x3, y1, y3, x1, y1, x2, y2, x3, y3, rec_fun);
      }
    }
  }

  template <typename VC>
  void calculate_full_slice(size_t start_index_left, size_t stop_index_left, size_t start_index_right, size_t stop_index_right,
    matrix_t& distances, cluster_container::Cluster_Container<VC>& cluster_left, cluster_container::Cluster_Container<VC>& cluster_right) {

    auto rec_fun = [&](size_t left, size_t right) {
      distances(left, right) = distance::dvstar<VC>(cluster_left.get(left), cluster_right.get(right));
    };

    utils::matrix_recursion(start_index_left, stop_index_left, start_index_right, stop_index_right, rec_fun);
  }

  //--------------------------------//
  // For inter-directory comparison //
  //--------------------------------//
  template <typename VC>
  matrix_t calculate_distances(cluster_container::Cluster_Container<VC>& cluster, size_t requested_cores) {

    matrix_t distances = matrix_t::Constant(cluster.size(), cluster.size(), 0);

    auto fun = [&](int x1, int y1, int x2, int y2, int x3, int y3) {
      calculate_triangle_slice<VC>(x1, y1, x2, y2, x3, y3, distances,
        cluster, cluster);
    };

    parallel::parallelize_triangle(cluster.size(), fun, requested_cores);
    return distances;
  }

  //-------------------------------//
  // For comparing two directories //
  //-------------------------------//
  template <typename VC>
  matrix_t calculate_distances(
    cluster_container::Cluster_Container<VC>& cluster_left,
    cluster_container::Cluster_Container<VC>& cluster_right, size_t requested_cores) {

    matrix_t distances{ cluster_left.size(), cluster_right.size() };

    auto fun = [&](size_t start_index_left, size_t stop_index_left, size_t start_index_right, size_t stop_index_right) {
      calculate_full_slice<VC>(start_index_left, stop_index_left, start_index_right, stop_index_right, std::ref(distances),
        std::ref(cluster_left), std::ref(cluster_right));
    };

    parallel::parallelize(cluster_left.size(), cluster_right.size(), fun, requested_cores);
    return distances;
  }

  void calculate_kmer_buckets(
    cluster_container::Kmer_Cluster& cluster_left, cluster_container::Kmer_Cluster& cluster_right,
    int left_offset, int right_offset, matrix_t& distances) {
    matrix_t dot_prod = matrix_t::Zero(cluster_left.size(), cluster_right.size());
    matrix_t left_norm = matrix_t::Zero(cluster_left.size(), cluster_right.size());
    matrix_t right_norm = matrix_t::Zero(cluster_left.size(), cluster_right.size());

    auto left_it = cluster_left.get_begin();
    auto left_end = cluster_left.get_end();
    auto left_size = cluster_left.experimental_bucket_count();
    for (int i = 0; i < left_size; i++) {
      auto idx = left_it->first;
      auto right_it = cluster_right.find(idx);
      if (right_it != cluster_right.get_end()) {
        distance::dvstar_kmer_major(left_it->second, right_it->second, dot_prod, left_norm, right_norm);
      }
      left_it++;
    }

    for (int x = 0; x < dot_prod.rows(); x++) {
      for (int y = 0; y < dot_prod.cols(); y++) {
        distances(x + left_offset, y + right_offset) = distance::normalise_dvstar(dot_prod(x, y), left_norm(x, y), right_norm(x, y));
      }
    }
  }

  //---------------------------//
  // Kmer-major implementation //
  //---------------------------//
  matrix_t calculate_distance_major(
    std::vector<cluster_container::Kmer_Cluster>& cluster_left,
    std::vector<cluster_container::Kmer_Cluster>& cluster_right, size_t nr_cores_to_use) {

    auto cluster_left_size = 0;
    std::vector<int> cluster_left_offsets{};
    for (int i = 0; i < cluster_left.size(); i++) {
      cluster_left_offsets.push_back(cluster_left_size);
      cluster_left_size += cluster_left[i].size();
    }

    auto cluster_right_size = 0;
    std::vector<int> cluster_right_offsets{};
    for (int i = 0; i < cluster_right.size(); i++) {
      cluster_right_offsets.push_back(cluster_right_size);
      cluster_right_size += cluster_right[i].size();
    }

    matrix_t distances = matrix_t::Zero(cluster_left_size, cluster_right_size);

    auto fun = [&](size_t start_index_left, size_t stop_index_left, size_t start_index_right, size_t stop_index_right) {
      for (auto left_i = start_index_left; left_i < stop_index_left; left_i++) {
        for (auto right_i = start_index_right; right_i < stop_index_right; right_i++) {
          calculate_kmer_buckets(cluster_left[left_i], cluster_right[right_i], cluster_left_offsets[left_i], cluster_right_offsets[right_i], distances);
        }
      }
    };

    parallel::parallelize(cluster_left.size(), cluster_right.size(), fun, nr_cores_to_use);

    return distances;
  }
}