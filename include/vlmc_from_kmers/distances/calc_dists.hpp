#pragma once

#include <Eigen/Core>
#include <mutex>

#include "cluster_container.hpp"
#include "distances/dvstar.hpp"
#include "global_aliases.hpp"
#include "parallel.hpp"
#include "utils.hpp"
#include "vlmc_container.hpp"

namespace vlmc::calc_dist {

template <typename VC>
void calculate_triangle_slice(
    int x1, int y1, int x2, int y2, int x3, int y3, matrix_t &distances,
    vlmc::container::ClusterContainer<VC> &cluster_left,
    vlmc::container::ClusterContainer<VC> &cluster_right) {

  auto rec_fun = [&](int left, int right) {
    if (distances(left, right) == 0) {
      distances(left, right) = distance::dvstar<VC>(cluster_left.get(left),
                                                    cluster_right.get(right));
    }
  };

  if (x1 < x2) {
    if (y2 > y3) {
      utils::triangle_recursion(x1, x3, y1, y2, x1, y1, x2, y2, x3, y3,
                                rec_fun);
    } else {
      utils::triangle_recursion(x1, x3, y1, y3, x1, y1, x2, y2, x3, y3,
                                rec_fun);
    }
  } else {
    if (y2 > y3) {
      utils::triangle_recursion(x2, x3, y1, y2, x1, y1, x2, y2, x3, y3,
                                rec_fun);
    } else {
      utils::triangle_recursion(x2, x3, y1, y3, x1, y1, x2, y2, x3, y3,
                                rec_fun);
    }
  }
}

template <typename VC>
void calculate_isocles_triangle_slice(
    size_t start_height, size_t end_height, size_t size, matrix_t &distances,
    vlmc::container::ClusterContainer<VC> &cluster) {

  for (size_t i = start_height; i < end_height; i++) {
    for (size_t j = 0; j < i; j++) {
      distances(i, j) = distance::dvstar<VC>(cluster.get(i), cluster.get(j));
    }
  }
}

template <typename VC>
void calculate_full_slice(size_t start_index_left, size_t stop_index_left,
                          size_t start_index_right, size_t stop_index_right,
                          matrix_t &distances,
                          container::ClusterContainer<VC> &cluster_left,
                          container::ClusterContainer<VC> &cluster_right) {

  auto rec_fun = [&](size_t left, size_t right) {
    distances(left, right) =
        distance::dvstar<VC>(cluster_left.get(left), cluster_right.get(right));
  };

  utils::matrix_recursion(start_index_left, stop_index_left, start_index_right,
                          stop_index_right, rec_fun);
}

//--------------------------------//
// For inter-directory comparison //
//--------------------------------//
template <typename VC>
matrix_t calculate_distances(vlmc::container::ClusterContainer<VC> &cluster,
                             size_t requested_cores) {

  matrix_t distances = matrix_t::Constant(cluster.size(), cluster.size(), 0);

  auto fun = [&](const size_t start_height, const size_t end_height,
                 const size_t size) {
    calculate_isocles_triangle_slice<VC>(start_height, end_height, size,
                                         distances, cluster);
  };

  parallel::parallelize_isocles_triangle(cluster.size(), fun, requested_cores);

  // We've only calculated the lower triangle.
  distances.triangularView<Eigen::Upper>() =
      distances.triangularView<Eigen::Lower>().transpose();

  return distances;
}

//-------------------------------//
// For comparing two directories //
//-------------------------------//
template <typename VC>
matrix_t
calculate_distances(vlmc::container::ClusterContainer<VC> &cluster_left,
                    vlmc::container::ClusterContainer<VC> &cluster_right,
                    size_t requested_cores) {

  matrix_t distances{cluster_left.size(), cluster_right.size()};

  auto fun = [&](size_t start_index_left, size_t stop_index_left,
                 size_t start_index_right, size_t stop_index_right) {
    calculate_full_slice<VC>(
        start_index_left, stop_index_left, start_index_right, stop_index_right,
        std::ref(distances), std::ref(cluster_left), std::ref(cluster_right));
  };

  parallel::parallelize(cluster_left.size(), cluster_right.size(), fun,
                        requested_cores);
  return distances;
}
} // namespace vlmc::calc_dist