#pragma once

#include <functional>
#include <filesystem>
#include <limits.h>
#include <exception>
#include <algorithm>
#include <execution>
#include <math.h>
#include <Eigen/Core>

#include "kmer.hpp"
#include "read_in_kmer.hpp"
#include "global_aliases.hpp"
#include "unordered_dense.h"

#include "vlmc_containers/veb_array.hpp"
#include "vlmc_containers/eytzinger_array.hpp"
#include "vlmc_containers/b_tree_array.hpp"

namespace vlmc_container {
  using RI_Kmer = kmers::RI_Kmer;

  int load_VLMCs_from_file(const std::filesystem::path& path_to_bintree, eigenx_t& cached_context,
    const std::function<void(const RI_Kmer& kmer)> f, const size_t background_order = 0) {
    std::ifstream ifs(path_to_bintree, std::ios::binary);
    cereal::BinaryInputArchive archive(ifs);
    kmers::VLMCKmer input_kmer{};

    int offset_to_remove = 0;
    for (int i = 0; i < background_order; i++) {
      offset_to_remove += std::pow(4, i);
    }

    while (ifs.peek() != EOF) {
      archive(input_kmer);
      RI_Kmer ri_kmer{ input_kmer };
      if (input_kmer.length <= background_order) {
        if (input_kmer.length + 1 > background_order) {
          int offset = ri_kmer.integer_rep - offset_to_remove;
          for (int x = 0; x < 4; x++) {
            cached_context(offset, x) = ri_kmer.next_char_prob[x];
          }
          // cached_context.row(offset) = ri_kmer.next_char_prob;
        }
      }
      else {
        f(ri_kmer);
      }
    }
    ifs.close();

    return offset_to_remove;
  }

  /*
    Storing Kmers in a sorted vector.
  */
  class VLMC_sorted_vector {

  public:
    std::vector<RI_Kmer> container{};
    VLMC_sorted_vector() = default;
    ~VLMC_sorted_vector() = default;

    VLMC_sorted_vector(const std::filesystem::path& path_to_bintree, const size_t background_order = 0, bool use_new = false) {
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer& kmer) { push(kmer); };

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      std::sort(std::execution::seq, container.begin(), container.end());
      for (size_t i = 0; i < size(); i++) {
        RI_Kmer kmer = get(i);
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        for (int x = 0; x < 4; x++) {
          get(i).next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
        }
      }
    }

    size_t size() const { return container.size(); }

    void push(const RI_Kmer& kmer) { container.push_back(kmer); }

    std::vector<RI_Kmer>::iterator begin() { return container.begin(); };
    std::vector<RI_Kmer>::iterator end() { return container.end(); };

    RI_Kmer& get(const int i) { return container[i]; }
  };

  void iterate_kmers(VLMC_sorted_vector& left_kmers, VLMC_sorted_vector& right_kmers, const std::function<void(const RI_Kmer& left_kmer, const RI_Kmer& right_kmer)>& f) {
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
      }
      else if (left_kmer < right_kmer) {
        ++left_it;
      }
      else
        ++right_it;
    }
  }

  /*
    Storing Kmers in a unordered map (HashMap).
  */
  class VLMC_hashmap {

  public:
    ankerl::unordered_dense::map<int, RI_Kmer> container{};
    VLMC_hashmap() = default;
    ~VLMC_hashmap() = default;

    VLMC_hashmap(const std::filesystem::path& path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer& kmer) { push(kmer); };

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      for (auto& [i_rep, kmer] : container) {
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        for (int x = 0; x < 4; x++) {
          kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
        }
      }
    }

    size_t size() const { return container.size(); }

    void push(const RI_Kmer& kmer) { container[kmer.integer_rep] = kmer; }

    RI_Kmer& get(const int i) { return container[i]; }
  };

  void iterate_kmers(VLMC_hashmap& left_kmers, VLMC_hashmap& right_kmers, const std::function<void(const RI_Kmer& left_kmer, const RI_Kmer& right_kmer)>& f) {
    for (auto& [i_rep, left_kmer] : left_kmers.container) {
      auto res = right_kmers.container.find(i_rep);
      if (res != right_kmers.container.end()) {
        auto right_kmer = res->second;
        f(left_kmer, right_kmer);
      }
    }
  }

  class VLMC_Veb {

  public:
    array::Veb_array* veb;
    VLMC_Veb() = default;
    ~VLMC_Veb() = default;

    VLMC_Veb(const std::filesystem::path& path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto tmp_container = std::vector<RI_Kmer>{};
      auto fun = [&](const RI_Kmer& kmer) { tmp_container.push_back(kmer); };

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
      for (auto& kmer : tmp_container) {
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        for (int x = 0; x < 4; x++) {
          kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
        }
      }
      veb = new array::Veb_array(tmp_container);
    }

    size_t size() const { return veb->n + 1; }

    RI_Kmer& get(const int i) {
      return veb->get_from_array(i);
    }
  };

  void iterate_kmers(VLMC_Veb& left_kmers, VLMC_Veb& right_kmers, const std::function<void(const RI_Kmer& left_kmer, const RI_Kmer& right_kmer)>& f) {
    int i = 0;
    while (i < left_kmers.veb->n) {
      RI_Kmer& left_kmer = left_kmers.veb->a[i];
      RI_Kmer& right_kmer = right_kmers.get(left_kmer.integer_rep);
      if (left_kmer == right_kmer) {
        f(left_kmer, right_kmer);
      }
      i++;
    }
  }

  class VLMC_Eytzinger {

  public:
    array::Ey_array* arr;
    VLMC_Eytzinger() = default;
    ~VLMC_Eytzinger() = default;

    VLMC_Eytzinger(const std::filesystem::path& path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto tmp_container = std::vector<RI_Kmer>{};
      auto fun = [&](const RI_Kmer& kmer) { tmp_container.push_back(kmer); };

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
      for (auto& kmer : tmp_container) {
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        for (int x = 0; x < 4; x++) {
          kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
        }
      }
      arr = new array::Ey_array(tmp_container);
    }

    size_t size() const { return arr->size + 1; }

    RI_Kmer& get(const int i) {
      ;
      return arr->get_from_array(i);
    }
  };

  void iterate_kmers(VLMC_Eytzinger& left_kmers, VLMC_Eytzinger& right_kmers, const std::function<void(const RI_Kmer& left_kmer, const RI_Kmer& right_kmer)>& f) {
    int i = 0;
    while (i <= left_kmers.arr->size) {
      RI_Kmer& left_kmer = left_kmers.arr->ey_sorted_kmers[i];
      RI_Kmer& right_kmer = right_kmers.get(left_kmer.integer_rep);
      if (left_kmer == right_kmer) {
        f(left_kmer, right_kmer);
      }
      i++;
    }
  }

  class VLMC_B_tree {
  public:
    array::B_Tree* arr;
    VLMC_B_tree() = default;
    ~VLMC_B_tree() = default;

    VLMC_B_tree(const std::filesystem::path& path_to_bintree, const size_t background_order = 0) {
      // cached_context : pointer to array which for each A, C, T, G has the next char probs
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto tmp_container = std::vector<RI_Kmer>{};
      auto fun = [&](const RI_Kmer& kmer) { tmp_container.push_back(kmer); };

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      std::sort(std::execution::seq, tmp_container.begin(), tmp_container.end());
      for (auto& kmer : tmp_container) {
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        for (int x = 0; x < 4; x++) {
          kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
        }
      }
      arr = new array::B_Tree(tmp_container);
    }

    size_t size() const { return arr->size + 1; }

    RI_Kmer& get(const int i) {
      ;
      return arr->get_from_array(i);
    }
  };

  void iterate_kmers(VLMC_B_tree& left_kmers, VLMC_B_tree& right_kmers, const std::function<void(const RI_Kmer& left_kmer, const RI_Kmer& right_kmer)>& f) {
    int i = 0;
    while (i < left_kmers.arr->size) {
      RI_Kmer& left_kmer = left_kmers.arr->a[i];
      RI_Kmer& right_kmer = right_kmers.get(left_kmer.integer_rep);
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

  class VLMC_sorted_search {

  private:
    RI_Kmer null_kmer{};
    int skip_size;

  public:
    std::vector<RI_Kmer> container{};
    std::vector<Min_max_node> summary{};
    int place_in_summary = 0;
    VLMC_sorted_search() = default;
    ~VLMC_sorted_search() = default;

    VLMC_sorted_search(const std::filesystem::path& path_to_bintree, const size_t background_order = 0, bool use_new = false) {
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const RI_Kmer& kmer) { push(kmer); };

      int offset_to_remove = load_VLMCs_from_file(path_to_bintree, cached_context, fun, background_order);

      std::sort(std::execution::seq, container.begin(), container.end());
      for (size_t i = 0; i < size(); i++) {
        RI_Kmer kmer = get(i);
        int background_idx = kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        for (int x = 0; x < 4; x++) {
          get(i).next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
        }
      }
      // Build summary
      if (container.size() > 0) {
        skip_size = std::ceil(std::log2(container.size()));
        if (skip_size == 0) {
          skip_size = 1;
        }
        summary.reserve(container.size() / skip_size);
        int i = 0;
        for (; i < container.size() - skip_size; i += skip_size) {
          summary.push_back(Min_max_node(i, container[i + skip_size - 1].integer_rep));
        }
        summary.push_back(Min_max_node(i, container[size() - 1].integer_rep));
      }
    }

    size_t size() const { return container.size(); }

    void push(const RI_Kmer& kmer) { container.push_back(kmer); }

    std::vector<RI_Kmer>::iterator begin() { return container.begin(); };
    std::vector<RI_Kmer>::iterator end() { return container.end(); };

    RI_Kmer& get(const int i) { return container[i]; }

    int find_block_start(int i_rep) {
      for (int i = place_in_summary; i < summary.size(); i++) {
        if (i_rep <= summary[i].max) {
          place_in_summary = i;
          return summary[i].block_start;
        }
      }
      return container.size();
    }
  };

  void iterate_kmers(VLMC_sorted_search& left_kmers, VLMC_sorted_search& right_kmers, const std::function<void(const RI_Kmer& left_kmer, const RI_Kmer& right_kmer)>& f) {
    left_kmers.place_in_summary = 0;
    right_kmers.place_in_summary = 0;

    auto left_i = 0;
    auto right_i = 0;
    auto left_size = left_kmers.size();
    auto right_size = right_kmers.size();

    while (left_i < left_size && right_i < right_size) {
      RI_Kmer& left_kmer = left_kmers.get(left_i);
      RI_Kmer& right_kmer = right_kmers.get(right_i);
      if (left_kmer == right_kmer) {
        f(left_kmer, right_kmer);
        ++left_i;
        ++right_i;
      }
      else if (left_kmer < right_kmer) {
        if (left_kmers.summary[left_kmers.place_in_summary].max < right_kmer.integer_rep) {
          left_i = left_kmers.find_block_start(right_kmer.integer_rep);
        }
        else {
          while (true) {
            ++left_i;
            RI_Kmer& left_kmer = left_kmers.get(left_i);
            if (left_kmer >= right_kmer || left_i >= left_size) {
              break;
            }
          }
        }
      }
      else {
        if (right_kmers.summary[right_kmers.place_in_summary].max < left_kmer.integer_rep) {
          right_i = right_kmers.find_block_start(left_kmer.integer_rep);
        }
        else {
          while (true) {
            ++right_i;
            RI_Kmer& right_kmer = right_kmers.get(right_i);
            if (right_kmer >= left_kmer || right_i >= right_size) {
              break;
            }
          }
        }
      }
    }
  }
}