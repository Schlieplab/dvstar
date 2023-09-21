#pragma once

#include <chrono>
#include <execution>
#include <filesystem>
#include <iostream>
#include <string>

#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>

#include <kmc_file.h>

#include "cli_helper.hpp"
#include "context_archive.hpp"
#include "estimators.hpp"
#include "kmc_runner.hpp"
#include "kmer_container.hpp"
#include "negative_log_likelihood.hpp"
#include "similarity_pruning.hpp"
#include "support_pruning.hpp"

bool LOGGING = true;

namespace vlmc {

void configure_stxxl(const std::filesystem::path &tmp_path) {
  std::string random_name = get_random_name("stxxl_");
  // get uninitialized config singleton
  stxxl::config *cfg = stxxl::config::get_instance();
  std::filesystem::path stxxl_disk_path = tmp_path / random_name;
  // create a disk_config structure.
  // allocate 6GB
  stxxl::disk_config disk1{stxxl_disk_path, uint64(6000) * 1024 * 1024,
                           "linuxaio unlink"};

  cfg->add_disk(disk1);
}

void similarity_pruning_steps(std::shared_ptr<KmerContainer<>> &container,
                              const Core &in_or_out_of_core,
                              const std::filesystem::path &out_path,
                              const estimator_f &remove_node) {
  std::ofstream file_stream(out_path, std::ios::binary);
  auto similarity_pruning_start = std::chrono::steady_clock::now();
  {
    cereal::BinaryOutputArchive oarchive(file_stream);

    if (in_or_out_of_core == Core::hash) {
      similarity_pruning_hash<max_k>(container, oarchive, remove_node);
    } else {
      similarity_pruning<max_k>(container, oarchive, remove_node);
    }
  }
  if (!out_path.empty()) {
    file_stream.close();
  }
}

int build_vlmc_from_kmc_db(
    const std::filesystem::path &kmc_db_path, const int max_depth,
    const int min_count, const double threshold,
    const std::filesystem::path &out_path, const Core &in_or_out_of_core,
    const double pseudo_count_amount = 1.0,
    const Estimator estimator = Estimator::kullback_leibler,
    const SequencingParameters sequencing_parameters = SequencingParameters{
        false}) {
  auto start = std::chrono::steady_clock::now();

  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(kmc_db_path);

  if (!status) {
    std::cerr << "ERROR: Opening kmc db not successful.  Try removing the file "
                 "extension."
              << std::endl;
    return EXIT_FAILURE;
  }

  auto include_node = [&](int length, size_t count) -> bool {
    if (sequencing_parameters.adjust_for_sequencing_errors) {
      const double adjusted_min_count =
          double(min_count) * adjusted_value(sequencing_parameters, length);

      return length <= max_depth && double(count) >= adjusted_min_count;
    } else {
      return length <= max_depth && count >= min_count;
    }
  };

  auto container =
      parse_kmer_container<ReverseKMerComparator<max_k>>(in_or_out_of_core);

  auto support_pruning_start = std::chrono::steady_clock::now();
  auto root = sequential_support_pruning<max_k>(kmer_database, container,
                                                max_depth + 1, include_node);
  auto support_pruning_done = std::chrono::steady_clock::now();

  kmer_database.Close();

  auto sorting_start = std::chrono::steady_clock::now();

  double sequence_length = root.count;

  container->sort();

  auto sorting_done = std::chrono::steady_clock::now();

  auto out_path_directory = out_path.parent_path();
  if (!out_path_directory.empty()) {
    std::filesystem::create_directories(out_path_directory);
  }

  auto similarity_pruning_start = std::chrono::steady_clock::now();

  estimator_f remove_node;
  if (estimator == Estimator::kullback_leibler) {
    remove_node =
        kl_estimator(threshold, pseudo_count_amount, sequencing_parameters);
  } else {
    remove_node = peres_shield_estimator(sequence_length, pseudo_count_amount,
                                         sequencing_parameters);
  }

  similarity_pruning_steps(container, in_or_out_of_core, out_path, remove_node);
  auto similarity_pruning_done = std::chrono::steady_clock::now();

  if (LOGGING) {
    std::chrono::duration<double> support_seconds =
        support_pruning_done - support_pruning_start;
    std::cout << "Support pruning time: " << support_seconds.count() << "s\n";

    std::chrono::duration<double> sort_seconds = sorting_done - sorting_start;
    std::cout << "Sorting time: " << sort_seconds.count() << "s\n";

    std::chrono::duration<double> similarity_seconds =
        similarity_pruning_done - sorting_done;
    std::cout << "Similarity pruning time: " << similarity_seconds.count()
              << "s\n";

    std::chrono::duration<double> total_seconds =
        similarity_pruning_done - start;
    std::cout << "VLMC steps total time: " << total_seconds.count() << "s\n";
  }

  return EXIT_SUCCESS;
}

int build_vlmc(const std::filesystem::path &fasta_path, const int max_depth,
               const int min_count, const double threshold,
               const std::filesystem::path &out_path,
               const std::filesystem::path &tmp_path,
               const Core &in_or_out_of_core,
               const double pseudo_count_amount = 1.0,
               const Estimator estimator = Estimator::kullback_leibler,
               const SequencingParameters sequencing_parameters =
                   SequencingParameters{false}) {
  auto start = std::chrono::steady_clock::now();

  const int kmer_size = max_depth + 1;

  auto kmc_db_path =
      run_kmc(fasta_path, kmer_size, tmp_path, in_or_out_of_core);

  auto kmc_done = std::chrono::steady_clock::now();

  if (LOGGING) {
    std::chrono::duration<double> kmc_seconds = kmc_done - start;
    std::cout << "KMC time: " << kmc_seconds.count() << "s\n";
  }

  auto status = build_vlmc_from_kmc_db(
      kmc_db_path, max_depth, min_count, threshold, out_path, in_or_out_of_core,
      pseudo_count_amount, estimator, sequencing_parameters);
  remove_kmc_files(kmc_db_path);

  return status;
}

int dump_path(const std::filesystem::path &in_path,
              const std::filesystem::path &out_path) {
  std::ifstream file_stream(in_path, std::ios::binary);
  cereal::BinaryInputArchive iarchive(file_stream);

  std::ostream *ofs = &std::cout;
  std::ofstream out_stream(out_path);

  if (out_path.empty()) {
    ofs = &std::cout;
  } else {
    ofs = &out_stream;
  }

  vlmc::iterate_archive(in_path,
                        [&](const VLMCKmer &kmer) { kmer.output(*ofs); });

  out_stream.close();

  return EXIT_SUCCESS;
}

int reprune_vlmc(const std::filesystem::path &in_path,
                 const std::filesystem::path &out_path,
                 const Core &in_or_out_of_core, const double new_threshold,
                 const double pseudo_count_amount = 1.0) {

  auto container =
      parse_kmer_container<ReverseKMerComparator<max_k>>(in_or_out_of_core);

  VLMCKmer root{};
  vlmc::iterate_archive(in_path, [&](const VLMCKmer &kmer) {
    if (kmer.length == 0) {
      root = kmer;
    }
    container->push(kmer);
  });

  container->sort();

  auto kl_remove_node = kl_estimator(new_threshold, pseudo_count_amount);
  auto ps_remove_node = peres_shield_estimator(root.count, pseudo_count_amount);

  similarity_pruning_steps(container, in_or_out_of_core, out_path,
                           kl_remove_node);

  return EXIT_SUCCESS;
}

} // namespace vlmc
