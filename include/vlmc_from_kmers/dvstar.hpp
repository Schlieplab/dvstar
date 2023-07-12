#pragma once

#include <filesystem>

#include "context_archive.hpp"
#include "kmer_container.hpp"
#include "similarity_pruning.hpp"

namespace vlmc::details {

// Have to decide a couple of things:
// 1. Do we allow for a separate VLMC for the  background model or not?
// 2. Do we disregard $k$-mers that are shorter than the background order?

std::array<std::array<double, 4>, 2>
get_components(const VLMCKmer &left, const VLMCKmer &left_background,
               const VLMCKmer &right, const VLMCKmer &right_background,
               const double &pseudo_count_amount) {
  auto left_probs = get_next_symbol_probabilities(left, pseudo_count_amount);
  auto left_probs_background =
      get_next_symbol_probabilities(left_background, pseudo_count_amount);
  auto right_probs = get_next_symbol_probabilities(right, pseudo_count_amount);
  auto right_probs_background =
      get_next_symbol_probabilities(right_background, pseudo_count_amount);

  return {std::array<double, 4>{
              left_probs[0] / std::sqrt(left_probs_background[0]),
              left_probs[1] / std::sqrt(left_probs_background[1]),
              left_probs[2] / std::sqrt(left_probs_background[2]),
              left_probs[3] / std::sqrt(left_probs_background[3])},
          std::array<double, 4>{
              right_probs[0] / std::sqrt(right_probs_background[0]),
              right_probs[1] / std::sqrt(right_probs_background[1]),
              right_probs[2] / std::sqrt(right_probs_background[2]),
              right_probs[3] / std::sqrt(right_probs_background[3])}};
}

VLMCKmer &find_background(std::ifstream &fs,
                          cereal::BinaryInputArchive &archive,
                          const VLMCKmer &kmer, VLMCKmer &background_kmer,
                          int background_order) {
  if (!has_next(fs)) {
    return background_kmer;
  }

  auto diff_pos = VLMCKmer::get_first_differing_position(
      kmer, background_kmer, kmer.length - background_order);

  while (diff_pos != -1 || background_kmer.length != background_order) {
    //    background_kmer.output(std::cout);
    if (has_next(fs)) {
      next(archive, background_kmer);
    } else {
      // TODO So this must be the last background, lets just use that for now...
      return background_kmer;
    }
    diff_pos = VLMCKmer::get_first_differing_position(
        kmer, background_kmer, kmer.length - background_order);
  }

  return background_kmer;
}

VLMCKmer &find_background(std::vector<VLMCKmer> &kmers, const VLMCKmer &kmer,
                          uint &background_i, int background_order) {
  auto diff_pos = VLMCKmer::get_first_differing_position(
      kmer, kmers[background_i], kmer.length - background_order);

  while (diff_pos != -1 || kmers[background_i].length != background_order) {
    if (background_i < kmers.size()) {
      background_i++;
    } else {
      // TODO So this must be the last background, lets just use that for now...
      return kmers[background_i];
    }
    diff_pos = VLMCKmer::get_first_differing_position(
        kmer, kmers[background_i], kmer.length - background_order);
  }

  return kmers[background_i];
}

void advance(std::ifstream &fs, cereal::BinaryInputArchive &archive,
             VLMCKmer &kmer, bool &iterating, uint &prev_length) {
  prev_length = kmer.length;
  if (has_next(fs)) {
    next(archive, kmer);
  } else {
    iterating = false;
  }
}

double transform_to_distance(double dot_product, double left_norm,
                             double right_norm) {
  left_norm = std::sqrt(left_norm);
  right_norm = std::sqrt(right_norm);

  if (left_norm == 0 || right_norm == 0) {
    return 1.0;
  } else {
    double Dvstar = dot_product / (left_norm * right_norm);

    double dvstar = 0.5 * (1 - Dvstar);

    double angular_distance = 2 * std::acos(Dvstar) / M_PI;
    return angular_distance;
  }
}
} // namespace vlmc::details

namespace vlmc {
std::vector<VLMCKmer> get_sorted_kmers(const std::filesystem::path &path) {
  std::vector<VLMCKmer> kmers{};
  iterate_archive(path, [&](const VLMCKmer &kmer) { kmers.push_back(kmer); });

  std::sort(std::execution::par_unseq, kmers.begin(), kmers.end(),
            ReverseKMerComparator<255>());
  return kmers;
}

/**
 * Applies f on the shared k-mers from the two lists of k-mers given.
 * @param left_kmers List of k-mers in a VLMC.
 * @param right_kmers List of k-mers in a VLMC.
 * @param f function applied to the shared k-mers of the two VLMCs.
 */
void iterate_shared_kmers_sorted(
    std::vector<VLMCKmer> &left_kmers, std::vector<VLMCKmer> &right_kmers,
    const std::function<void(const VLMCKmer &left, const VLMCKmer &right)> &f,
    const std::function<void(const VLMCKmer &left, const VLMCKmer &right)>
        &f_not_shared) {
  uint left_i = 0;
  uint right_i = 0;

  while (left_i < left_kmers.size() && right_i < right_kmers.size()) {
    const VLMCKmer &left_kmer = left_kmers[left_i];
    const VLMCKmer &right_kmer = right_kmers[right_i];
    if (right_kmer == left_kmer) {
      f(left_kmer, right_kmer);
      left_i++;
      right_i++;
    } else if (left_kmer.reverse_less_than(right_kmer)) {
      f_not_shared(left_kmer, right_kmer);
      left_i++;
    } else {
      f_not_shared(right_kmer, left_kmer);
      right_i++;
    }
  }
}

/**
 * Applies f on the shared k-mers from the two .bintree files stored at
 * left_path and right_path.
 * @param left_path Path to a .bintree file.
 * @param right_path Path to a .bintree file.
 * @param f function applied to the shared k-mers of the two VLMCs.
 */
void iterate_shared_kmers_sorted(
    const std::filesystem::path &left_path,
    const std::filesystem::path &right_path,
    const std::function<void(const VLMCKmer &left, const VLMCKmer &right)> &f,
    const std::function<void(const VLMCKmer &left, const VLMCKmer &right)>
        &f_not_shared) {

  auto left_kmers = get_sorted_kmers(left_path);
  auto right_kmers = get_sorted_kmers(right_path);
  iterate_shared_kmers_sorted(left_kmers, right_kmers, f, f_not_shared);
}

/**
 * Does not work.  Idea is to iterate the two lists of k-mers stored at
 * left_path and right_path without loading the two lists into memory and
 * re-sorting them. Currently doesn't work because of the funky output order
 * from similarity pruning.
 * @param left_path Path to a .bintree file.
 * @param right_path Path to a .bintree file.
 * @param f function applied to the shared k-mers of the two VLMCs.
 */
void iterate_shared_kmers(
    const std::filesystem::path &left_path,
    const std::filesystem::path &right_path,
    const std::function<void(const VLMCKmer &left, const VLMCKmer &right)> &f) {
  // TODO: Don't use this, it doesn't work.
  // But I do think there should be a way to do this without re-sorting the
  // output.

  std::ifstream left_fs(left_path, std::ios::binary);
  std::ifstream right_fs(right_path, std::ios::binary);

  cereal::BinaryInputArchive left_archive(left_fs);
  cereal::BinaryInputArchive right_archive(right_fs);

  VLMCKmer left_kmer{};
  VLMCKmer right_kmer{};
  VLMCKmer suffix_kmer{};
  bool iterating = true;

  if (!details::has_next(left_fs) || !details::has_next(right_fs)) {
    return;
  }

  details::next(left_archive, left_kmer);
  details::next(right_archive, right_kmer);

  auto last_equal_length = left_kmer.length;
  auto prev_left_length = left_kmer.length;
  auto prev_right_length = right_kmer.length;

  while (iterating) {
    std::cout << left_kmer.to_string() << ", " << right_kmer.to_string();

    // Finding the shared k-mers would be super easy if the lists were
    // sorted.  Then, it would be as simple as checking which of each not equal
    // pair is smallest and popping accordingly.
    // However, since the similarity pruning does not output
    // k-mers in a sorted order it gets a little bit trickier.

    // 1. Check if k-mers are equal, if they are, call f and proceed.
    // 2. If the k-mers are not equal, we could be in a number of situations.
    //    If one is longer than the last equal k-mer that would indicate that
    //    it is part of a branch that the other is missing. This entire
    //    branch needs to be popped before we can proceed.  We detect a popped
    //    branch the length being the same as the other k-mer.
    //    If
    if (right_kmer == left_kmer) {
      f(left_kmer, right_kmer);

      prev_left_length = left_kmer.length;
      prev_right_length = right_kmer.length;

      details::advance(right_fs, right_archive, right_kmer, iterating,
                       prev_right_length);
      details::advance(left_fs, left_archive, left_kmer, iterating,
                       prev_left_length);

      if (left_kmer.length < right_kmer.length) {

        while (left_kmer.length < right_kmer.length) {
          details::advance(right_fs, right_archive, right_kmer, iterating,
                           prev_right_length);
        }
      } else if (right_kmer.length < left_kmer.length) {
        while (right_kmer.length < left_kmer.length) {
          details::advance(left_fs, left_archive, left_kmer, iterating,
                           prev_left_length);
        }
      }
    } else if (left_kmer.reverse_less_than(right_kmer)) {
      details::advance(left_fs, left_archive, left_kmer, iterating,
                       prev_left_length);

    } else {
      details::advance(right_fs, right_archive, right_kmer, iterating,
                       prev_right_length);
    }
  }

  left_fs.close();
  right_fs.close();
}

/**
 * Calculate the dvstar metric between the two VLMCs consisting of the two given
 * lists of k-mers.  Background order determines how long the background is.
 * @param left_kmers List of k-mers in a VLMC.
 * @param right_kmers List of k-mers in a VLMC.
 * @param background_order Background order of the measure, 0 corresponds to
 * GC-adjusted, 1 to dinucleotide-adjusted etc.
 * @return the dvstar dissimilarity between the two VLMCs
 */
double dvstar(std::vector<VLMCKmer> &left_kmers,
              std::vector<VLMCKmer> &right_kmers, int background_order = 0,
              const double &pseudo_count_amount = 1) {
  double dot_product = 0.0;
  double left_norm = 0.0;
  double right_norm = 0.0;

  uint left_background_i = 0;
  uint right_background_i = 0;

  iterate_shared_kmers_sorted(
      left_kmers, right_kmers,
      [&](const VLMCKmer &left_kmer, const VLMCKmer &right_kmer) {
        if (left_kmer.length <= background_order) {
          return;
        }
        auto &left_kmer_background = details::find_background(
            left_kmers, left_kmer, left_background_i, background_order);
        auto &right_kmer_background = details::find_background(
            right_kmers, right_kmer, right_background_i, background_order);

        auto [left_comp, right_comp] =
            details::get_components(left_kmer, left_kmer_background, right_kmer,
                                    right_kmer_background, pseudo_count_amount);

        for (int i = 0; i < 4; i++) {
          dot_product += left_comp[i] * right_comp[i];
          left_norm += std::pow(left_comp[i], 2.0);
          right_norm += std::pow(right_comp[i], 2.0);
        }
      },
      [](auto l, auto r) {});

  return details::transform_to_distance(dot_product, left_norm, right_norm);
}

/**
 * Calculate the dvstar metric between the two VLMCs from the two given paths.
 * Background order determines how long the background is.
 * @param left_path Path to a .bintree file.
 * @param right_path Path to a .bintree file.
 * @param background_order Background order of the measure, 0 corresponds to
 * GC-adjusted, 1 to dinucleotide-adjusted etc.
 * @return the dvstar dissimilarity between the two VLMCs
 */
std::tuple<double, double>
dvstar(const std::filesystem::path &left_path,
       const std::filesystem::path &right_path,
       const std::function<void(const VLMCKmer &left, const VLMCKmer &right)>
           &f_not_shared,
       int background_order = 0, const double &pseudo_count_amount = 1) {
  std::ifstream left_fs_background(left_path, std::ios::binary);
  std::ifstream right_fs_background(right_path, std::ios::binary);

  double dot_product = 0.0;
  double left_norm = 0.0;
  double right_norm = 0.0;
  double n_shared_kmers = 0.0;

  {
    cereal::BinaryInputArchive left_archive_background(left_fs_background);
    cereal::BinaryInputArchive right_archive_background(right_fs_background);

    VLMCKmer left_kmer_background{};
    VLMCKmer right_kmer_background{};
    VLMCKmer suffix_kmer{};

    if (!details::has_next(left_fs_background) ||
        !details::has_next(right_fs_background)) {
      return {1.0, 0.0};
    }

    details::next(left_archive_background, left_kmer_background);
    details::next(right_archive_background, right_kmer_background);

    iterate_shared_kmers_sorted(
        left_path, right_path,
        [&](const VLMCKmer &left_kmer, const VLMCKmer &right_kmer) {
          if (left_kmer.length <= background_order) {
            return;
          }
          n_shared_kmers++;
          details::find_background(left_fs_background, left_archive_background,
                                   left_kmer, left_kmer_background,
                                   background_order);
          details::find_background(right_fs_background,
                                   right_archive_background, right_kmer,
                                   right_kmer_background, background_order);

          auto [left_comp, right_comp] = details::get_components(
              left_kmer, left_kmer_background, right_kmer,
              right_kmer_background, pseudo_count_amount);

          for (int i = 0; i < 4; i++) {
            dot_product += left_comp[i] * right_comp[i];
            left_norm += std::pow(left_comp[i], 2.0);
            right_norm += std::pow(right_comp[i], 2.0);
          }
        },
        f_not_shared);
  }

  left_fs_background.close();
  left_fs_background.close();

  return {details::transform_to_distance(dot_product, left_norm, right_norm),
          n_shared_kmers};
}

/**
 * Calculate the dvstar metric between the two VLMCs from the two given paths.
 * Background order determines how long the background is.
 * @param left_path Path to a .bintree file.
 * @param right_path Path to a .bintree file.
 * @param background_order Background order of the measure, 0 corresponds to
 * GC-adjusted, 1 to dinucleotide-adjusted etc.
 * @return the dvstar dissimilarity between the two VLMCs
 */
double dvstar(const std::filesystem::path &left_path,
              const std::filesystem::path &right_path, int background_order = 0,
              const double &pseudo_count_amount = 1) {
  auto [diss, _n_shared] = dvstar(
      left_path, right_path, [](auto l, auto r) {}, background_order,
      pseudo_count_amount);
  return diss;
}

/**
 * Calculate the penalized dvstar metric between the two VLMCs from the two
 * given paths. THe penalizsation refers to the added dissimilarity based on how
 * many of the contexts are missing from the models. Background order determines
 * how long the background is.
 * @param left_path Path to a .bintree file.
 * @param right_path Path to a .bintree file.
 * @param background_order Background order of the measure, 0 corresponds to
 * GC-adjusted, 1 to dinucleotide-adjusted etc.
 * @return the dvstar dissimilarity between the two VLMCs
 */
double dvstar_missing_penalized(const std::filesystem::path &left_path,
                                const std::filesystem::path &right_path,
                                int background_order = 0,
                                const double &pseudo_count_amount = 1) {
  double missing_contexts = 0;
  double total_contexts = 0;
  auto [dvstar_dist, n_shared] = dvstar(
      left_path, right_path,
      [&](auto l, auto r) {
        missing_contexts++;
        total_contexts++;
      },
      background_order, pseudo_count_amount);

  double missing_frac = missing_contexts / (total_contexts + n_shared);
  // Ensure the scale is from 0-1 with equal weight on missing as dvstar.
  return dvstar_dist * 0.75 + missing_frac * 0.25;
}

} // namespace vlmc
