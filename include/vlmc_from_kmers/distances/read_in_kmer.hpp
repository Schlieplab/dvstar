#pragma once

#include <Eigen/Core>
#include <filesystem>
#include <functional>
#include <memory>
#include <numeric>

#include "../kmer.hpp"
#include "global_aliases.hpp"

namespace vlmc {
struct ReadInKmer {
  int integer_rep{};
  std::array<out_t, 4> next_char_prob{};

  ReadInKmer() = default;
  ~ReadInKmer() = default;

  explicit ReadInKmer(const VLMCKmer &old_kmer,
                      const double pseudo_count_amount = 1.0) {
    out_t child_count =
        old_kmer.next_symbol_counts[0] + old_kmer.next_symbol_counts[1] +
        old_kmer.next_symbol_counts[2] + old_kmer.next_symbol_counts[3] +
        pseudo_count_amount * 4;
    this->next_char_prob = {
        (old_kmer.next_symbol_counts[0] + pseudo_count_amount) / child_count,
        (old_kmer.next_symbol_counts[1] + pseudo_count_amount) / child_count,
        (old_kmer.next_symbol_counts[2] + pseudo_count_amount) / child_count,
        (old_kmer.next_symbol_counts[3] + pseudo_count_amount) / child_count};
    this->integer_rep = get_index_rep(old_kmer);
  }

  ReadInKmer(const int vlmc_rep) { this->integer_rep = vlmc_rep; }

  int get_index_rep(const VLMCKmer &kmer) const {
    int integer_value = 0;
    int offset = 1;
    for (int i = kmer.length - 1; i >= 0; i--) {
      auto kmer_2_bits = extract2bits(kmer, i) + 1;
      integer_value += (kmer_2_bits * offset);
      offset *= 4;
    }
    return integer_value;
  }

  inline char extract2bits(const VLMCKmer &kmer, unsigned int pos) const {
    char row = pos >> 5;
    char pos_in_row = pos & 31;
    char n_shift_pos_to_end = (62 - pos_in_row * 2);
    return (kmer.kmer_data[row] >> n_shift_pos_to_end) & 3;
  }

  static int background_order_index(int integer_rep_, int order) {
    if (integer_rep_ < std::pow(4, order))
      return integer_rep_;
    int back_rep = 0;
    int i = 1;
    for (int o = 0; o < order; o++) {
      int r = integer_rep_ % 4;
      if (r == 0)
        r = 4;
      integer_rep_ = (integer_rep_ - r) / 4;
      back_rep += r * i;
      i *= 4;
    }
    return back_rep;
  }

  inline bool operator<(const ReadInKmer &kmer) const {
    return this->integer_rep < kmer.integer_rep;
  };
  inline bool operator>(const ReadInKmer &kmer) const {
    return this->integer_rep > kmer.integer_rep;
  };
  inline bool operator>=(const ReadInKmer &kmer) const {
    return this->integer_rep >= kmer.integer_rep;
  };
  inline bool operator<=(const ReadInKmer &kmer) const {
    return this->integer_rep <= kmer.integer_rep;
  };
  inline bool operator==(const ReadInKmer &kmer) const {
    return this->integer_rep == kmer.integer_rep;
  };
  inline bool operator!=(const ReadInKmer &kmer) const {
    return this->integer_rep != kmer.integer_rep;
  };
};
} // namespace vlmc
