#pragma once

#include <Eigen/Core>
#include <functional>
#include <memory>
#include <filesystem>
#include <numeric>

#include "kmer.hpp"
#include "global_aliases.hpp"

constexpr out_t pseudo_count_amount = 1.0;

namespace kmers {
  struct RI_Kmer {
    int integer_rep;
    std::array<out_t, 4> next_char_prob;

    RI_Kmer() = default;
    ~RI_Kmer() = default;

    RI_Kmer(const kmers::VLMCKmer& old_kmer) {
      out_t child_count = old_kmer.next_symbol_counts[0] + old_kmer.next_symbol_counts[1] + old_kmer.next_symbol_counts[2] + old_kmer.next_symbol_counts[3] + 4;
      this->next_char_prob = { (old_kmer.next_symbol_counts[0] + pseudo_count_amount) / child_count,
                     (old_kmer.next_symbol_counts[1] + pseudo_count_amount) / child_count,
                     (old_kmer.next_symbol_counts[2] + pseudo_count_amount) / child_count,
                     (old_kmer.next_symbol_counts[3] + pseudo_count_amount) / child_count };
      this->integer_rep = get_index_rep(old_kmer);
    }

    RI_Kmer(const int vlmc_rep) {
      this->integer_rep = vlmc_rep;
    }

    int get_index_rep(const kmers::VLMCKmer& kmer) {
      int integer_value = 0;
      int offset = 1;
      for (int i = kmer.length - 1; i >= 0; i--) {
        auto kmer_2_bits = extract2bits(kmer, i) + 1;
        integer_value += (kmer_2_bits * offset);
        offset *= 4;
      }
      return integer_value;
    }

    inline char extract2bits(const kmers::VLMCKmer& kmer, unsigned int pos) const {
      char row = pos >> 5;
      char pos_in_row = pos & 31;
      char n_shift_pos_to_end = (62 - pos_in_row * 2);
      return (kmer.kmer_data[row] >> n_shift_pos_to_end) & 3;
    }

    int background_order_index(int integer_rep, int order) {
      if (integer_rep < std::pow(4, order))
        return integer_rep;
      int back_rep = 0;
      int i = 1;
      for (int o = 0; o < order; o++) {
        int r = integer_rep % 4;
        if (r == 0)
          r = 4;
        integer_rep = (integer_rep - r) / 4;
        back_rep += r * i;
        i *= 4;
      }
      return back_rep;
    }

    inline bool operator<(const RI_Kmer& kmer) const {
      return this->integer_rep < kmer.integer_rep;
    };
    inline bool operator>(const RI_Kmer& kmer) const {
      return this->integer_rep > kmer.integer_rep;
    };
    inline bool operator>=(const RI_Kmer& kmer) const {
      return this->integer_rep >= kmer.integer_rep;
    };
    inline bool operator<=(const RI_Kmer& kmer) const {
      return this->integer_rep <= kmer.integer_rep;
    };
    inline bool operator==(const RI_Kmer& kmer) const {
      return this->integer_rep == kmer.integer_rep;
    };
    inline bool operator!=(const RI_Kmer& kmer) const {
      return this->integer_rep != kmer.integer_rep;
    };
  };
}
