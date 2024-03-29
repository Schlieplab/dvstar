#pragma once

#include <array>
#include <filesystem>
#include <numeric>
#include <vector>

#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/array.hpp>

#include "estimators.hpp"
#include "kmer.hpp"
#include "kmer_container.hpp"
#include "kmers_per_level.hpp"

namespace vlmc {

// Need to keep track of which kmers have children, as those with children
// can't be removed.
struct PstKmer {
  VLMCKmer kmer;
  bool has_children = false;
  bool real_child = false;
  bool is_terminal = true;
};

std::tuple<bool, bool> process_parent(const VLMCKmer &kmer,
                                      KMersPerLevel<PstKmer> &kmers_per_level,
                                      cereal::BinaryOutputArchive &oarchive,
                                      const estimator_f &remove_node) {
  auto &children = kmers_per_level[kmer.length + 1];

  int removed_children = 0;
  int n_real_children = 0;

  for (auto &[child, has_children, real_child, is_terminal] : children) {
    if (!real_child) {
      continue;
    }
    n_real_children++;

    if (has_children) {
      // Can't remove nodes with children.
      child.divergence = -1.0;
      child.is_terminal = is_terminal;
      oarchive(child);
    } else {
      auto [remove_child, divergence] = remove_node(child, kmer);

      if (remove_child) {
        removed_children++;
      } else {
        child.divergence = divergence;
        child.is_terminal = is_terminal;
        oarchive(child);
      }
    }
  }

  // Reset children
  kmers_per_level.reset_depth(kmer.length + 1);

  return {removed_children != n_real_children,
          n_real_children - removed_children != 4};
}

void similarity_prune(const VLMCKmer &prev_kmer, const VLMCKmer &kmer,
                      KMersPerLevel<PstKmer> &kmers_per_level,
                      cereal::BinaryOutputArchive &oarchive,
                      const estimator_f &remove_node) {
  bool has_children = true;
  bool is_terminal = true;

  if (prev_kmer.length > kmer.length) {
    // Has different length, kmer must be parent of the previous nodes
    auto [has_children_, is_terminal_] =
        process_parent(kmer, kmers_per_level, oarchive, remove_node);
    has_children = has_children_;
    is_terminal = is_terminal_;
  } else {
    // Has same length or shorter, so the kmers have to be children
    // of the same node, or the first kmer in a new set of children.
    has_children = false;
    is_terminal = true;
  }

  if (kmer.length == 0) {
    // Root node, should be the last node, and always outputted.
    oarchive(kmer);
  } else {
    // 3 - char_pos which ensures the ordering of the outputted k-mers remain
    // in lexicographically descending order.
    // Only important for post-processing for e.g. dvstar.
    auto start_char_pos = 3 - kmer.char_pos(0);
    kmers_per_level[kmer.length][start_char_pos] = {kmer, has_children, true,
                                                    is_terminal};
  }
}

template <int kmer_size>
void similarity_pruning(std::shared_ptr<KmerContainer<>> &container,
                        cereal::BinaryOutputArchive &oarchive,
                        const estimator_f &remove_node) {

  KMersPerLevel<PstKmer> kmers_per_level{4, kmer_size + 1};

  VLMCTranslator kmer_api(kmer_size);
  VLMCKmer prev_kmer = kmer_api.construct_vlmc_kmer();

  // Example container output:
  // TTTTTTTT
  // GTTTTTTT
  // CTTTTTTT
  // ATTTTTTT
  // TTTTTTT
  // TGTTTTTT

  container->for_each([&](VLMCKmer &kmer) {
    similarity_prune(prev_kmer, kmer, kmers_per_level, oarchive, remove_node);

    prev_kmer = kmer;
  });
}

template <int kmer_size>
void similarity_pruning_hash(std::shared_ptr<KmerContainer<>> &container,
                             cereal::BinaryOutputArchive &oarchive,
                             const estimator_f &remove_node) {

  KMersPerLevel<PstKmer> kmers_per_level{4, kmer_size + 1};

  VLMCKmer fake_parent_kmer{kmer_size, 0, {}};

  bool found_one = false;

  container->for_each([&](VLMCKmer &kmer) {
    if (kmer.length != 0) {
      VLMCKmer::create_suffix_kmer(kmer, fake_parent_kmer);
      auto &parent_kmer = container->get(fake_parent_kmer);
      auto [remove_child, divergence] = remove_node(kmer, parent_kmer);

      kmer.divergence = divergence;

      if (!kmer.has_children && remove_child) {
        kmer.to_be_removed = true;
      } else {
        parent_kmer.has_children = true;
        parent_kmer.to_be_removed = false;
        kmer.is_terminal = true;
        kmer.to_be_removed = false;

        for (int i = parent_kmer.length - 1; i > 0; i--) {
          VLMCKmer::create_suffix_kmer(fake_parent_kmer, fake_parent_kmer);
          auto &grandparent_kmer = container->get(fake_parent_kmer);

          grandparent_kmer.has_children = true;
          grandparent_kmer.to_be_removed = false;
        }
      }
    }
  });

  container->for_each([&](VLMCKmer &kmer) {
    if (kmer.length == 0) {
      oarchive(kmer);
    } else {
      int substr_pos = 0;
      if (!kmer.to_be_removed || kmer.has_children) {
        oarchive(kmer);
      }
    }
  });
}

} // namespace vlmc
