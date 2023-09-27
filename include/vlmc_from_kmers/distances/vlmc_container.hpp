#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <climits>
#include <cmath>
#include <exception>
#include <execution>
#include <filesystem>
#include <functional>

#include <ankerl/unordered_dense.h>

#include "../kmer.hpp"
#include "global_aliases.hpp"
#include "read_in_kmer.hpp"

#include "containers/b_tree_array.hpp"
#include "containers/eytzinger_array.hpp"
#include "containers/veb_array.hpp"

namespace vlmc::container {

int load_VLMCs_from_file(const std::filesystem::path &path_to_bintree,
                         eigenx_t &cached_context,
                         const std::function<void(const ReadInKmer &kmer)> f,
                         const size_t background_order = 0,
                         const double pseudo_count_amount = 1.0) {

  int offset_to_remove = 0;
  for (int i = 0; i < background_order; i++) {
    offset_to_remove += std::pow(4, i);
  }

  iterate_archive(path_to_bintree, [&](const VLMCKmer &input_kmer) {
    ReadInKmer ri_kmer{input_kmer, pseudo_count_amount};
    if (input_kmer.length <= background_order) {
      if (input_kmer.length + 1 > background_order) {
        int offset = ri_kmer.integer_rep - offset_to_remove;
        for (int x = 0; x < 4; x++) {
          cached_context(offset, x) = ri_kmer.next_char_prob[x];
        }
        // cached_context.row(offset) = ri_kmer.next_char_prob;
      }
    } else {
      f(ri_kmer);
    }
  });

  return offset_to_remove;
}

/*
  Storing Kmers in a sorted vector.
*/
class SortedVector {

public:
  std::vector<ReadInKmer> vector{};
  SortedVector() = default;
  ~SortedVector() = default;

  explicit SortedVector(const std::filesystem::path &path_to_bintree,
                        const size_t background_order = 0,
                        const double pseudo_count_amount = 1.0,
                        bool use_new = false) {
    eigenx_t cached_context((int)std::pow(4, background_order), 4);

    auto fun = [&](const ReadInKmer &kmer) { push(kmer); };

    int offset_to_remove =
        load_VLMCs_from_file(path_to_bintree, cached_context, fun,
                             background_order, pseudo_count_amount);

    std::sort(std::execution::seq, vector.begin(), vector.end());
    for (size_t i = 0; i < size(); i++) {
      ReadInKmer kmer = get(i);
      int background_idx = ReadInKmer::background_order_index(kmer.integer_rep,
                                                              background_order);
      int offset = background_idx - offset_to_remove;
      for (int x = 0; x < 4; x++) {
        get(i).next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
      }
    }
  }

  size_t size() const { return vector.size(); }

  void push(const ReadInKmer &kmer) { vector.push_back(kmer); }

  std::vector<ReadInKmer>::iterator begin() { return vector.begin(); };
  std::vector<ReadInKmer>::iterator end() { return vector.end(); };

  ReadInKmer &get(const int i) { return vector[i]; }
};

void iterate_kmers(SortedVector &left_kmers, SortedVector &right_kmers,
                   const std::function<void(const ReadInKmer &left_kmer,
                                            const ReadInKmer &right_kmer)> &f) {
  auto right_it = right_kmers.begin();
  auto right_end = right_kmers.end();
  auto left_it = left_kmers.begin();
  auto left_end = left_kmers.end();

  while (left_it != left_end && right_it != right_end) {
    auto left_kmer = *left_it;
    auto right_kmer = *right_it;
    if (left_kmer == right_kmer) {
      f(left_kmer, right_kmer);
      ++left_it;
      ++right_it;
    } else if (left_kmer < right_kmer) {
      ++left_it;
    } else
      ++right_it;
  }
}

/*
  Storing Kmers in an unordered map (HashMap).
*/
class HashMap {

public:
  ankerl::unordered_dense::map<int, ReadInKmer> map{};
  HashMap() = default;
  ~HashMap() = default;

  explicit HashMap(const std::filesystem::path &path_to_bintree,
                   const size_t background_order = 0,
                   const double pseudo_count_amount = 1.0) {
    // cached_context : pointer to array which for each A, C, T, G has the next
    // char probs
    eigenx_t cached_context((int)std::pow(4, background_order), 4);

    auto fun = [&](const ReadInKmer &kmer) { push(kmer); };

    int offset_to_remove =
        load_VLMCs_from_file(path_to_bintree, cached_context, fun,
                             background_order, pseudo_count_amount);

    for (auto &[i_rep, kmer] : map) {
      int background_idx =
          kmer.background_order_index(kmer.integer_rep, background_order);
      int offset = background_idx - offset_to_remove;
      for (int x = 0; x < 4; x++) {
        kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
      }
    }
  }

  size_t size() const { return map.size(); }

  void push(const ReadInKmer &kmer) { map[kmer.integer_rep] = kmer; }

  ReadInKmer &get(const int i) { return map[i]; }
};

void iterate_kmers(HashMap &left_kmers, HashMap &right_kmers,
                   const std::function<void(const ReadInKmer &left_kmer,
                                            const ReadInKmer &right_kmer)> &f) {
  for (auto &[i_rep, left_kmer] : left_kmers.map) {
    auto res = right_kmers.map.find(i_rep);
    if (res != right_kmers.map.end()) {
      auto right_kmer = res->second;
      f(left_kmer, right_kmer);
    }
  }
}

class VanEmdeBoasTree {

public:
  array::VanEmdeBoasArray *veb;
  VanEmdeBoasTree() = default;
  ~VanEmdeBoasTree() = default;

  explicit VanEmdeBoasTree(const std::filesystem::path &path_to_bintree,
                           const size_t background_order = 0,
                           const double pseudo_count_amount = 1.0) {
    // cached_context : pointer to array which for each A, C, T, G has the next
    // char probs
    eigenx_t cached_context((int)std::pow(4, background_order), 4);

    auto tmp_container = std::vector<ReadInKmer>{};
    auto fun = [&](const ReadInKmer &kmer) { tmp_container.push_back(kmer); };

    int offset_to_remove =
        load_VLMCs_from_file(path_to_bintree, cached_context, fun,
                             background_order, pseudo_count_amount);

    std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
    for (auto &kmer : tmp_container) {
      int background_idx =
          kmer.background_order_index(kmer.integer_rep, background_order);
      int offset = background_idx - offset_to_remove;
      for (int x = 0; x < 4; x++) {
        kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
      }
    }
    veb = new array::VanEmdeBoasArray(tmp_container);
  }

  size_t size() const { return veb->n + 1; }

  ReadInKmer &get(const int i) const { return veb->get_from_array(i); }
};

void iterate_kmers(VanEmdeBoasTree &left_kmers, VanEmdeBoasTree &right_kmers,
                   const std::function<void(const ReadInKmer &left_kmer,
                                            const ReadInKmer &right_kmer)> &f) {
  int i = 0;
  while (i < left_kmers.veb->n) {
    ReadInKmer &left_kmer = left_kmers.veb->a[i];
    ReadInKmer &right_kmer = right_kmers.get(left_kmer.integer_rep);
    if (left_kmer == right_kmer) {
      f(left_kmer, right_kmer);
    }
    i++;
  }
}

class EytzingerTree {

public:
  array::EytzingerArray *arr;
  EytzingerTree() = default;
  ~EytzingerTree() = default;

  explicit EytzingerTree(const std::filesystem::path &path_to_bintree,
                         const size_t background_order = 0,
                         const double pseudo_count_amount = 1.0) {
    // cached_context : pointer to array which for each A, C, T, G has the next
    // char probs
    eigenx_t cached_context((int)std::pow(4, background_order), 4);

    auto tmp_container = std::vector<ReadInKmer>{};
    auto fun = [&](const ReadInKmer &kmer) { tmp_container.push_back(kmer); };

    int offset_to_remove =
        load_VLMCs_from_file(path_to_bintree, cached_context, fun,
                             background_order, pseudo_count_amount);

    std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
    for (auto &kmer : tmp_container) {
      int background_idx =
          kmer.background_order_index(kmer.integer_rep, background_order);
      int offset = background_idx - offset_to_remove;
      for (int x = 0; x < 4; x++) {
        kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
      }
    }
    arr = new array::EytzingerArray(tmp_container);
  }

  size_t size() const { return arr->size + 1; }

  ReadInKmer &get(const int i) const { return arr->get_from_array(i); }
};

void iterate_kmers(EytzingerTree &left_kmers, EytzingerTree &right_kmers,
                   const std::function<void(const ReadInKmer &left_kmer,
                                            const ReadInKmer &right_kmer)> &f) {
  int i = 0;
  while (i <= left_kmers.arr->size) {
    ReadInKmer &left_kmer = left_kmers.arr->ey_sorted_kmers[i];
    ReadInKmer &right_kmer = right_kmers.get(left_kmer.integer_rep);
    if (left_kmer == right_kmer) {
      f(left_kmer, right_kmer);
    }
    i++;
  }
}

class BTree {
public:
  array::B_Tree *arr;
  BTree() = default;
  ~BTree() = default;

  explicit BTree(const std::filesystem::path &path_to_bintree,
                 const size_t background_order = 0,
                 const double pseudo_count_amount = 1.0) {
    // cached_context : pointer to array which for each A, C, T, G has the next
    // char probs
    eigenx_t cached_context((int)std::pow(4, background_order), 4);

    auto tmp_container = std::vector<ReadInKmer>{};
    auto fun = [&](const ReadInKmer &kmer) { tmp_container.push_back(kmer); };

    int offset_to_remove =
        load_VLMCs_from_file(path_to_bintree, cached_context, fun,
                             background_order, pseudo_count_amount);

    std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
    for (auto &kmer : tmp_container) {
      int background_idx =
          kmer.background_order_index(kmer.integer_rep, background_order);
      int offset = background_idx - offset_to_remove;
      for (int x = 0; x < 4; x++) {
        kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
      }
    }
    arr = new array::B_Tree(tmp_container);
  }

  size_t size() const { return arr->size + 1; }

  ReadInKmer &get(const int i) const { return arr->get_from_array(i); }
};

void iterate_kmers(BTree &left_kmers, BTree &right_kmers,
                   const std::function<void(const ReadInKmer &left_kmer,
                                            const ReadInKmer &right_kmer)> &f) {
  int i = 0;
  while (i < left_kmers.arr->size) {
    ReadInKmer &left_kmer = left_kmers.arr->a[i];
    ReadInKmer &right_kmer = right_kmers.get(left_kmer.integer_rep);
    if (left_kmer == right_kmer) {
      f(left_kmer, right_kmer);
    }
    i++;
  }
}

/*
  Storing Kmers in a sorted vector with a summary structure to skip past misses.
*/

struct Min_max_node {
  int block_start;
  int max;

  Min_max_node(int idx, int max) {
    this->block_start = idx;
    this->max = max;
  }

  Min_max_node() = default;
  ~Min_max_node() = default;
};

class SortedSearch {

private:
  ReadInKmer null_kmer{};
  int skip_size;

public:
  std::vector<ReadInKmer> container_{};
  std::vector<Min_max_node> summary{};
  int place_in_summary = 0;
  SortedSearch() = default;
  ~SortedSearch() = default;

  explicit SortedSearch(const std::filesystem::path &path_to_bintree,
                        const size_t background_order = 0,
                        const double pseudo_count_amount = 1.0,
                        bool use_new = false) {
    eigenx_t cached_context((int)std::pow(4, background_order), 4);

    auto fun = [&](const ReadInKmer &kmer) { push(kmer); };

    int offset_to_remove =
        load_VLMCs_from_file(path_to_bintree, cached_context, fun,
                             background_order, pseudo_count_amount);

    std::sort(std::execution::seq, container_.begin(), container_.end());
    for (size_t i = 0; i < size(); i++) {
      ReadInKmer kmer = get(i);
      int background_idx =
          kmer.background_order_index(kmer.integer_rep, background_order);
      int offset = background_idx - offset_to_remove;
      for (int x = 0; x < 4; x++) {
        get(i).next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
      }
    }
    // Build summary
    if (!container_.empty()) {
      skip_size = std::ceil(std::log2(container_.size()));
      if (skip_size == 0) {
        skip_size = 1;
      }
      summary.reserve(container_.size() / skip_size);
      int i = 0;
      for (; i < container_.size() - skip_size; i += skip_size) {
        summary.emplace_back(i, container_[i + skip_size - 1].integer_rep);
      }
      summary.emplace_back(i, container_[size() - 1].integer_rep);
    }
  }

  size_t size() const { return container_.size(); }

  void push(const ReadInKmer &kmer) { container_.push_back(kmer); }

  std::vector<ReadInKmer>::iterator begin() { return container_.begin(); };
  std::vector<ReadInKmer>::iterator end() { return container_.end(); };

  ReadInKmer &get(const int i) { return container_[i]; }

  int find_block_start(int i_rep) {
    for (int i = place_in_summary; i < summary.size(); i++) {
      if (i_rep <= summary[i].max) {
        place_in_summary = i;
        return summary[i].block_start;
      }
    }
    return container_.size();
  }
};

void iterate_kmers(SortedSearch &left_kmers, SortedSearch &right_kmers,
                   const std::function<void(const ReadInKmer &left_kmer,
                                            const ReadInKmer &right_kmer)> &f) {
  left_kmers.place_in_summary = 0;
  right_kmers.place_in_summary = 0;

  auto left_i = 0;
  auto right_i = 0;
  auto left_size = left_kmers.size();
  auto right_size = right_kmers.size();

  while (left_i < left_size && right_i < right_size) {
    ReadInKmer &left_kmer = left_kmers.get(left_i);
    ReadInKmer &right_kmer = right_kmers.get(right_i);
    if (left_kmer == right_kmer) {
      f(left_kmer, right_kmer);
      ++left_i;
      ++right_i;
    } else if (left_kmer < right_kmer) {
      if (left_kmers.summary[left_kmers.place_in_summary].max <
          right_kmer.integer_rep) {
        left_i = left_kmers.find_block_start(right_kmer.integer_rep);
      } else {
        while (true) {
          ++left_i;
          ReadInKmer &left_kmer = left_kmers.get(left_i);
          if (left_kmer >= right_kmer || left_i >= left_size) {
            break;
          }
        }
      }
    } else {
      if (right_kmers.summary[right_kmers.place_in_summary].max <
          left_kmer.integer_rep) {
        right_i = right_kmers.find_block_start(left_kmer.integer_rep);
      } else {
        while (true) {
          ++right_i;
          ReadInKmer &right_kmer = right_kmers.get(right_i);
          if (right_kmer >= left_kmer || right_i >= right_size) {
            break;
          }
        }
      }
    }
  }
}
} // namespace vlmc::container