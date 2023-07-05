#pragma once

#include <math.h>

#include "vlmc_container.hpp"
#include "cluster_container.hpp"
#include "read_in_kmer.hpp"
#include "utils.hpp"
#include "global_aliases.hpp"

namespace distance {

  using bucket_t = std::vector<cluster_container::Kmer_Pair>;
  using RI_Kmer = kmers::RI_Kmer;

  out_t normalise_dvstar(out_t dot_product, out_t left_norm, out_t right_norm) {

    left_norm = std::sqrt(left_norm);
    right_norm = std::sqrt(right_norm);
    if (left_norm == 0 || right_norm == 0) {
      return 1.0;
    }
    else {
      out_t Dvstar = dot_product / (left_norm * right_norm);

      out_t dvstar = 0.5 * (1 - Dvstar);

      out_t angular_distance = 2 * std::acos(Dvstar) / M_PI;
      if (isnan(angular_distance)) {
        return 0.0;
      }
      else {
        return angular_distance;
      }
    }
  }

  template <typename VC>
  out_t dvstar(VC& left, VC& right) {
    out_t dot_product = 0.0;
    out_t left_norm = 0.0;
    out_t right_norm = 0.0;

    auto f = [&](const RI_Kmer& left_kmer, const RI_Kmer& right_kmer) {
      dot_product += ((left_kmer.next_char_prob[0] * right_kmer.next_char_prob[0]) +
        (left_kmer.next_char_prob[1] * right_kmer.next_char_prob[1]) +
        (left_kmer.next_char_prob[2] * right_kmer.next_char_prob[2]) +
        (left_kmer.next_char_prob[3] * right_kmer.next_char_prob[3]));
      left_norm += ((left_kmer.next_char_prob[0] * left_kmer.next_char_prob[0]) +
        (left_kmer.next_char_prob[1] * left_kmer.next_char_prob[1]) +
        (left_kmer.next_char_prob[2] * left_kmer.next_char_prob[2]) +
        (left_kmer.next_char_prob[3] * left_kmer.next_char_prob[3]));
      right_norm += ((right_kmer.next_char_prob[0] * right_kmer.next_char_prob[0]) +
        (right_kmer.next_char_prob[1] * right_kmer.next_char_prob[1]) +
        (right_kmer.next_char_prob[2] * right_kmer.next_char_prob[2]) +
        (right_kmer.next_char_prob[3] * right_kmer.next_char_prob[3]));
    };

    if (left.size() < right.size()) {
      vlmc_container::iterate_kmers(left, right, f);
    }
    else {
      vlmc_container::iterate_kmers(right, left, f);
    }

    return normalise_dvstar(dot_product, left_norm, right_norm);
  }

  void dvstar_kmer_major(bucket_t& left_vector, bucket_t& right_vector,
    matrix_t& dot_prod, matrix_t& left_norm, matrix_t& right_norm) {
    auto rec_fun = [&](size_t& left, size_t& right) {
      auto left_id = left_vector[left].id;
      auto right_id = right_vector[right].id;
      auto left_kmer = left_vector[left].kmer;
      auto right_kmer = right_vector[right].kmer;
      dot_prod(left_id, right_id) += ((left_kmer.next_char_prob[0] * right_kmer.next_char_prob[0]) +
        (left_kmer.next_char_prob[1] * right_kmer.next_char_prob[1]) +
        (left_kmer.next_char_prob[2] * right_kmer.next_char_prob[2]) +
        (left_kmer.next_char_prob[3] * right_kmer.next_char_prob[3]));
      left_norm(left_id, right_id) += ((left_kmer.next_char_prob[0] * left_kmer.next_char_prob[0]) +
        (left_kmer.next_char_prob[1] * left_kmer.next_char_prob[1]) +
        (left_kmer.next_char_prob[2] * left_kmer.next_char_prob[2]) +
        (left_kmer.next_char_prob[3] * left_kmer.next_char_prob[3]));
      right_norm(left_id, right_id) += ((right_kmer.next_char_prob[0] * right_kmer.next_char_prob[0]) +
        (right_kmer.next_char_prob[1] * right_kmer.next_char_prob[1]) +
        (right_kmer.next_char_prob[2] * right_kmer.next_char_prob[2]) +
        (right_kmer.next_char_prob[3] * right_kmer.next_char_prob[3]));
    };

    utils::matrix_recursion(0, left_vector.size(), 0, right_vector.size(), rec_fun);
  }
}