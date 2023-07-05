#pragma once

#include <stdlib.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

#include "utils.hpp"

namespace parallel {
  std::vector<std::tuple<size_t, size_t>> get_x_bounds(size_t size, const size_t requested_cores) {
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

  void recursive_get_triangle_coords(std::vector<std::array<int, 6>>& triangle_coords,
    int x1, int y1, int x2, int y2, int x3, int y3, const int remaining_cores) {
    if (remaining_cores == 1) {
      triangle_coords.push_back({ x1, y1, x2, y2, x3, y3 });
      return;
    }
    auto one_to_three = x1 == x3 || y1 == y3;
    auto two_to_three = x2 == x3 || y2 == y3;
    auto one_to_two = x1 == x2 || y1 == y2;

    if (two_to_three && one_to_two) {
      auto new_x = (x3 + x1 + 1) / 2;
      auto new_y = (y3 + y1 + 1) / 2;
      if (remaining_cores <= 2) {
        triangle_coords.push_back({ x1, y1, x2, y2, new_x, new_y });
        triangle_coords.push_back({ new_x, new_y, x2, y2, x3, y3 });
      }
      else {
        auto cores = remaining_cores - (remaining_cores / 2);
        recursive_get_triangle_coords(triangle_coords, x1, y1, x2, y2, new_x, new_y, remaining_cores / 2);
        recursive_get_triangle_coords(triangle_coords, new_x, new_y, x2, y2, x3, y3, cores);
      }
    }
    else if ((two_to_three && not one_to_three) || (one_to_three && one_to_two)) {
      auto new_x = (x3 + x2 + 1) / 2;
      auto new_y = (y3 + y2 + 1) / 2;
      if (remaining_cores <= 2) {
        triangle_coords.push_back({ x1, y1, x2, y2, new_x, new_y });
        triangle_coords.push_back({ x1, y1, new_x, new_y, x3, y3 });
      }
      else {
        auto cores = remaining_cores - (remaining_cores / 2);
        recursive_get_triangle_coords(triangle_coords, x1, y1, x2, y2, new_x, new_y, remaining_cores / 2);
        recursive_get_triangle_coords(triangle_coords, x1, y1, new_x, new_y, x3, y3, cores);
      }
    }
    else {
      auto new_x = (x2 + x1 + 1) / 2;
      auto new_y = (y2 + y1 + 1) / 2;
      if (remaining_cores <= 2) {
        triangle_coords.push_back({ x1, y1, new_x, new_y, x3, y3 });
        triangle_coords.push_back({ new_x, new_y, x2, y2, x3, y3 });
      }
      else {
        auto cores = remaining_cores - (remaining_cores / 2);
        recursive_get_triangle_coords(triangle_coords, x1, y1, new_x, new_y, x3, y3, remaining_cores / 2);
        recursive_get_triangle_coords(triangle_coords, new_x, new_y, x2, y2, x3, y3, cores);
      }
    }
  }

  void parallelize_triangle(size_t size, const std::function<void(int, int, int, int, int, int)>& fun, const size_t requested_cores) {
    std::vector<std::thread> threads{};
    std::vector<std::array<int, 6>> triangle_coords;
    size_t used_cores = utils::get_used_cores(requested_cores, size);
    recursive_get_triangle_coords(triangle_coords, 0, 0, 0, size, size, size, used_cores);

    for (auto& coords : triangle_coords) {
      threads.emplace_back(fun, coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]);
    }

    for (auto& thread : threads) {
      if (thread.joinable()) {
        thread.join();
      }
    }
  }

  void parallelize(size_t size, const std::function<void(size_t, size_t)>& fun, const size_t requested_cores) {
    std::vector<std::thread> threads{};

    auto bounds = get_x_bounds(size, requested_cores);
    for (auto& [start_index, stop_index] : bounds) {
      threads.emplace_back(fun, start_index, stop_index);
    }

    for (auto& thread : threads) {
      if (thread.joinable()) {
        thread.join();
      }
    }
  }

  void parallelize(size_t size_left, size_t size_right, const std::function<void(size_t, size_t, size_t, size_t)>& fun, const size_t requested_cores) {
    std::vector<std::thread> threads{};
    if (size_left > size_right) {
      auto bounds = get_x_bounds(size_left, requested_cores);
      for (auto& [start_index, stop_index] : bounds) {
        threads.emplace_back(fun, start_index, stop_index, 0, size_right);
      }
    }
    else {
      auto bounds = get_x_bounds(size_right, requested_cores);
      for (auto& [start_index, stop_index] : bounds) {
        threads.emplace_back(fun, 0, size_left, start_index, stop_index);
      }
    }

    for (auto& thread : threads) {
      if (thread.joinable()) {
        thread.join();
      }
    }
  }

  void parallelize_kmer_major(size_t size, const std::function<void(size_t, size_t, size_t)>& fun, size_t nr_cores_to_use) {
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
    for (auto& [start_index, stop_index] : bounds) {
      threads.emplace_back(fun, start_index, stop_index, idx);
      idx++;
    }

    for (auto& thread : threads) {
      if (thread.joinable()) {
        thread.join();
      }
    }
  }
}