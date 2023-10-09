#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <thread>

#include "cluster_container.hpp"
#include "global_aliases.hpp"
#include "read_in_kmer.hpp"

namespace utils {
void matrix_recursion(
    size_t start_index_left, size_t stop_index_left, size_t start_index_right,
    size_t stop_index_right,
    const std::function<void(size_t &left, size_t &right)> &fun) {

  std::queue<std::tuple<size_t, size_t, size_t, size_t>> queue{};
  queue.push(
      {start_index_left, stop_index_left, start_index_right, stop_index_right});

  while (!queue.empty()) {
    auto [start_left, stop_left, start_right, stop_right] = queue.front();
    auto diff_left = stop_left - start_left;
    auto diff_right = stop_right - start_right;

    if (diff_left == 1 && diff_right == 1) {
      fun(start_left, start_right);

    } else if (diff_right > diff_left) {
      auto new_right = (stop_right + start_right) / 2;
      queue.push({start_left, stop_left, start_right, new_right});
      queue.push({start_left, stop_left, new_right, stop_right});
    } else {
      auto new_left = (stop_left + start_left) / 2;
      queue.push({start_left, new_left, start_right, stop_right});
      queue.push({new_left, stop_left, start_right, stop_right});
    }

    queue.pop();
  }
}

size_t get_used_cores(size_t requested_cores, size_t size) {
  const size_t processor_count = std::thread::hardware_concurrency();
  size_t used_cores = 1;
  if (requested_cores > size) {
    used_cores = size;
  } else if (requested_cores <= processor_count) {
    used_cores = requested_cores;
  } else {
    used_cores = processor_count;
  }
  return used_cores;
}

void print_matrix(const matrix_t &distance_matrix,
                  const std::vector<std::string> &ids_from,
                  const std::vector<std::string> &ids_to) {
  for (size_t j = 0; j < distance_matrix.cols(); j++) {
    std::cout << ids_to[j] << "\t";
  }
  std::cout << std::endl;

  for (size_t i = 0; i < distance_matrix.rows(); i++) {
    std::cout << ids_from[i] << "\t";
    for (size_t j = 0; j < distance_matrix.cols(); j++) {
      std::cout << distance_matrix(i, j) << "\t";
    }
    std::cout << std::endl;
  }
}

int sign(int p1x, int p1y, int p2x, int p2y, int p3x, int p3y) {
  return (p1x - p3x) * (p2y - p3y) - (p2x - p3x) * (p1y - p3y);
}

bool is_point_in_triangle(int ptx, int pty, int v1x, int v1y, int v2x, int v2y,
                          int v3x, int v3y) {

  int d1 = sign(ptx, pty, v1x, v1y, v2x, v2y);
  int d2 = sign(ptx, pty, v2x, v2y, v3x, v3y);
  int d3 = sign(ptx, pty, v3x, v3y, v1x, v1y);

  bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
  bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

  return !(has_neg && has_pos);
}

void triangle_recursion(int start_index_left, int stop_index_left,
                        int start_index_right, int stop_index_right, int x1,
                        int y1, int x2, int y2, int x3, int y3,
                        const std::function<void(int &left, int &right)> &fun) {
  auto diff_left = stop_index_left - start_index_left;
  auto diff_right = stop_index_right - start_index_right;
  if (diff_left == 1 && diff_right == 1) {
    if (is_point_in_triangle(start_index_left, start_index_right, x1, y1, x2,
                             y2, x3, y3)) {
      fun(start_index_left, start_index_right);
    }
  } else if (diff_right > diff_left) {
    auto new_right_index = (stop_index_right + start_index_right) / 2;
    triangle_recursion(start_index_left, stop_index_left, start_index_right,
                       new_right_index, x1, y1, x2, y2, x3, y3, fun);
    triangle_recursion(start_index_left, stop_index_left, new_right_index,
                       stop_index_right, x1, y1, x2, y2, x3, y3, fun);
  } else {
    auto new_left_index = (stop_index_left + start_index_left) / 2;
    triangle_recursion(start_index_left, new_left_index, start_index_right,
                       stop_index_right, x1, y1, x2, y2, x3, y3, fun);
    triangle_recursion(new_left_index, stop_index_left, start_index_right,
                       stop_index_right, x1, y1, x2, y2, x3, y3, fun);
  }
}
} // namespace utils