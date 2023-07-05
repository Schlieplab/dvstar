#pragma once

#include <array>
#include <functional>
#include <unordered_map>

#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/array.hpp>

namespace kmers {
  using uint64 = unsigned long;
  using uint32 = unsigned int;

  struct VLMCKmer {
    VLMCKmer() = default;
    VLMCKmer(uint32 length_, uint64 count_,
      std::array<uint64, 4> next_symbol_counts_)
      : length(length_), count(count_), next_symbol_counts(next_symbol_counts_),
      divergence(-1.0), kmer_data(), is_terminal(false), has_children(false),
      to_be_removed(true) {
      this->n_rows = (length_ / 32.0) + 1; // std::ceil, but no float conversion
      if (length_ == 0) {
        this->n_rows = 0;
      }
    }

    ~VLMCKmer() = default;

    std::array<uint64, 4> kmer_data;
    uint64 count;
    std::array<uint64, 4> next_symbol_counts;
    double divergence;
    uint32 length;
    uint32 n_rows;
    bool is_terminal;
    bool has_children;
    bool to_be_removed;

    // This method lets cereal know which data members to serialize
    template <class Archive>
    void serialize(Archive& archive) {
      archive(kmer_data, length, n_rows, count, next_symbol_counts, divergence,
        is_terminal);
    }
  };
}