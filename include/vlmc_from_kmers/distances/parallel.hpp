#pragma once

#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

#include "utils.hpp"

namespace vlmc::parallel {
std::vector<std::tuple<size_t, size_t>>
get_x_bounds(size_t size, const size_t requested_cores) {
  size_t used_cores = utils::get_used_cores(requested_cores, size);
  std::vector<std::tuple<size_t, size_t>> bounds_per_thread{};
  float values_per_thread = float(size) / float(used_cores);

  auto limit = std::min(used_cores, size);
  for (size_t i = 0; i < limit; i++) {
    size_t start_index = std::floor(values_per_thread * i);
    size_t stop_index = std::floor(values_per_thread * (i + 1.0));
    if (i == (limit - 1)) {
      stop_index = size;
    }
    bounds_per_thread.emplace_back(start_index, stop_index);
  }

  return bounds_per_thread;
}

void parallelize_isocles_triangle(
    size_t size, const std::function<void(size_t, size_t, size_t)> &fun,
    const size_t requested_cores) {
  std::vector<std::thread> threads{};

  int used_cores =
      static_cast<int>(utils::get_used_cores(requested_cores, size));

  std::vector<size_t> triangle_coords(used_cores);

  double triangle_area = double(size) * double(size) / 2.0;
  double area_per_thread = triangle_area / double(used_cores);

  // The computation of the triangle is divided into #threads amount of ~equal
  // area (rounded to whole integers). The idea is that the triangle is divided
  // into a top triangle and then stacked trapezoids.

  triangle_coords[0] = std::round(std::sqrt(2 * area_per_thread));

  for (int i = 1; i < used_cores - 1; i++) {
    triangle_coords[i] = std::round(std::sqrt(
        2 * (area_per_thread + std::pow(triangle_coords[i - 1], 2.0) / 2)));
  }
  triangle_coords[triangle_coords.size() - 1] = size;

  for (int i = 0; i < triangle_coords.size(); i++) {
    if (i == 0) {
      threads.emplace_back(fun, 0, triangle_coords[i], size);
    } else {
      threads.emplace_back(fun, triangle_coords[i - 1], triangle_coords[i],
                           size);
    }
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void parallelize(size_t size, const std::function<void(size_t, size_t)> &fun,
                 const size_t requested_cores) {
  std::vector<std::thread> threads{};

  auto bounds = get_x_bounds(size, requested_cores);
  for (auto &[start_index, stop_index] : bounds) {
    threads.emplace_back(fun, start_index, stop_index);
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void parallelize(size_t size_left, size_t size_right,
                 const std::function<void(size_t, size_t, size_t, size_t)> &fun,
                 const size_t requested_cores) {
  std::vector<std::thread> threads{};
  if (size_left > size_right) {
    auto bounds = get_x_bounds(size_left, requested_cores);
    for (auto &[start_index, stop_index] : bounds) {
      threads.emplace_back(fun, start_index, stop_index, 0, size_right);
    }
  } else {
    auto bounds = get_x_bounds(size_right, requested_cores);
    for (auto &[start_index, stop_index] : bounds) {
      threads.emplace_back(fun, 0, size_left, start_index, stop_index);
    }
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void parallelize_kmer_major(
    size_t size, const std::function<void(size_t, size_t, size_t)> &fun,
    size_t nr_cores_to_use) {
  std::vector<std::thread> threads{};
  std::vector<std::tuple<size_t, size_t>> bounds{};
  float values_per_thread = (size + nr_cores_to_use - 1) / nr_cores_to_use;

  auto start_index = 0;
  while (start_index < size) {
    if (start_index + values_per_thread > size) {
      bounds.emplace_back(start_index, size);
      break;
    }
    bounds.emplace_back(start_index, start_index + values_per_thread);
    start_index += values_per_thread;
  }

  int idx = 0;
  for (auto &[start_index, stop_index] : bounds) {
    threads.emplace_back(fun, start_index, stop_index, idx);
    idx++;
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}
} // namespace vlmc::parallel