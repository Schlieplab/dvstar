#pragma once

#include <filesystem>
#include <format>
#include <mutex>
#include <thread>

#include "cluster_container.hpp"
#include "global_aliases.hpp"
#include "parallel.hpp"
#include "vlmc_from_kmers/io_helper.h"

namespace vlmc {


template <typename VC>
std::tuple<container::ClusterContainer<VC>, std::vector<std::string>>
get_cluster(const std::filesystem::path &path, size_t nr_cores_to_use,
            const size_t background_order,
            const double pseudo_count_amount = 1.0, const int set_size = -1) {
  auto paths = get_recursive_paths(path, {".bintree"});

  size_t paths_size = paths.size();
  if ((set_size != -1) && (set_size < paths_size)) {
    paths_size = set_size;
  }

  if (nr_cores_to_use > paths_size) {
    nr_cores_to_use = paths_size;
  }
  if (nr_cores_to_use > 4) {
    nr_cores_to_use = 4;
  }

  container::ClusterContainer<VC> cluster{paths_size};
  std::vector<std::string> ids{paths_size};

  auto fun = [&](size_t start_index, size_t stop_index) {
    for (int index = start_index; index < stop_index; index++) {
      cluster[index] = VC(paths[index], background_order, pseudo_count_amount);
      ids[index] = paths[index].stem();
    }
  };

  parallel::parallelize(paths_size, fun, nr_cores_to_use);

  return {cluster, ids};
}

std::vector<container::KmerCluster>
get_kmer_cluster(const std::filesystem::path &directory, size_t nr_cores_to_use,
                 const size_t background_order = 0, const int set_size = -1,
                 const double pseudo_count_amount = 1.0) {
  std::vector<std::filesystem::path> paths{};

  for (const auto &dir_entry : recursive_directory_iterator(directory)) {
    if (dir_entry.path().extension() == ".bintree") {
      paths.push_back(dir_entry.path());
    }
  }

  size_t paths_size = paths.size();
  if ((set_size != -1) && (set_size < paths_size)) {
    paths_size = set_size;
  }

  if (nr_cores_to_use > paths_size)
    nr_cores_to_use = paths_size;

  std::vector<container::KmerCluster> clusters(nr_cores_to_use);

  auto fun = [&](size_t start_index, size_t stop_index, size_t idx) {
    for (int index = start_index; index < stop_index; index++) {
      std::vector<ReadInKmer> input_vector{};
      eigenx_t cached_context((int)std::pow(4, background_order), 4);

      auto fun = [&](const ReadInKmer &kmer) { input_vector.push_back(kmer); };

      int offset_to_remove = container::load_VLMCs_from_file(
          paths[index], cached_context, fun, background_order,
          pseudo_count_amount);

      for (auto kmer : input_vector) {
        int background_idx =
            kmer.background_order_index(kmer.integer_rep, background_order);
        int offset = background_idx - offset_to_remove;
        for (int x = 0; x < 4; x++) {
          kmer.next_char_prob[x] *= 1.0 / std::sqrt(cached_context(offset, x));
        }
        clusters[idx].push(container::KmerPair{kmer, index - start_index});
      }
    }
    clusters[idx].set_size(stop_index - start_index);
  };

  parallel::parallelize_kmer_major(paths_size, fun, nr_cores_to_use);

  return clusters;
}
} // namespace vlmc
