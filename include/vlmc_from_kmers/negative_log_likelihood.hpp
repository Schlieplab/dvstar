#pragma once

#include <filesystem>
#include <memory>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <ankerl/unordered_dense.h>

#include "io_helper.h"
#include "kmc_runner.hpp"
#include "kmer.hpp"
#include "kmer_container.hpp"
#include "kmers_per_level.hpp"
#include "vlmc_from_kmers/distances/parallel.hpp"

namespace vlmc::details {
void load_kmers(const std::filesystem::path &vlmc_path,
                KmerContainer<KMerComparator<max_k>> &container) {

  vlmc::iterate_archive(vlmc_path,
                        [&](const VLMCKmer &kmer) { container.push(kmer); });
}

std::tuple<double, bool, bool> compare_kmers(VLMCKmer &vlmc_kmer,
                                             VLMCKmer &kmc_kmer,
                                             const int prefix_length) {
  // Need to find the vlmc kmer that best matches the given kmer
  auto diff_pos = VLMCKmer::get_first_differing_position(vlmc_kmer, kmc_kmer, 0,
                                                         prefix_length);
  // If the two kmers differ anywhere (that isn't the last position),
  // they don't match.
  // If they don't differ anywhere in the vlmc kmer, we've found a match.
  if (diff_pos == -1) {
    auto char_idx = kmc_kmer.char_pos(kmc_kmer.length - 1);

    // Remember the pseudo-counts.
    double probability = double(vlmc_kmer.next_symbol_counts[char_idx] + 1) /
                         double(vlmc_kmer.count + 4);
    double log_likelihood = std::log(probability) * double(kmc_kmer.count);

    return {log_likelihood, true, false};
  } else {
    auto kmc_char_idx = kmc_kmer.char_pos(diff_pos + prefix_length);
    auto vlmc_char_idx = vlmc_kmer.char_pos(diff_pos);

    if (kmc_char_idx < vlmc_char_idx) {
      // Advance the kmc kmer, but save it as haven't found
      // find its matching vlmc kmer yet.
      return {0.0, true, true};
    } else {
      // Advance to the next vlmc kmer
      return {0.0, false, false};
    }
  }
}

std::tuple<double, stxxl::vector<VLMCKmer>>
score_kmers(KmerContainer<KMerComparator<max_k>> &container,
            const stxxl::vector<VLMCKmer> &kmers, const int kmer_length,
            const int prefix_length) {
  if (kmers.empty()) {
    return {0.0, {}};
  }
  stxxl::vector<VLMCKmer> next_level_kmers{};
  double log_likelihood = 0.0;

  bool next_please = false;
  bool kmc_done = false;
  int kmers_idx = 0;

  VLMCKmer kmc_kmer = kmers[kmers_idx];

  container.for_each([&](VLMCKmer &vlmc_kmer) {
    if (kmer_length != vlmc_kmer.length || kmc_done) {
      // Iterate the vlmc-kmers one depth at a time.
      return;
    }

    do {
      auto [log_likelihood_, next_please_, save_kmer] =
          compare_kmers(vlmc_kmer, kmc_kmer, prefix_length);
      if (save_kmer) {
        next_level_kmers.push_back(kmc_kmer);
      }
      next_please = next_please_;
      log_likelihood += log_likelihood_;

      if (next_please) {
        if (kmers.size() <= kmers_idx + 1) {
          kmc_done = true;
          return;
        } else {
          kmers_idx++;
          kmc_kmer = kmers[kmers_idx];
        }
      }
    } while (next_please);
  });

  return {log_likelihood, next_level_kmers};
}

std::tuple<double, stxxl::vector<VLMCKmer>, size_t>
score_kmers_full_length(KmerContainer<KMerComparator<max_k>> &container,
                        const std::string &kmc_db_name,
                        const int actual_kmer_size) {
  CKMCFile kmer_database;
  auto status = kmer_database.OpenForListing(kmc_db_name);

  if (!status) {
    std::cout << "opening file not successful" << std::endl;
    throw std::invalid_argument("internal: kmc db path does not open.");
  }

  VLMCTranslator kmer_api(actual_kmer_size + 1);
  uint64 counter;
  VLMCKmer kmc_kmer(actual_kmer_size + 1, 0, {});

  status = kmer_database.ReadNextKmer(kmer_api, counter);
  kmc_kmer = kmer_api.construct_vlmc_kmer();
  kmc_kmer.count = counter;

  if (!status) {
    return {0.0, {}, 0};
  }

  VLMCKmer prev_vlmc_kmer{};
  double log_likelihood = 0.0;
  size_t sequence_length = counter;

  stxxl::vector<VLMCKmer> next_level_kmers{};

  bool next_please = false;

  container.for_each([&](VLMCKmer &vlmc_kmer) {
    if (actual_kmer_size != vlmc_kmer.length) {
      // Iterate the vlmc-kmers one depth at a time.
      return;
    }

    do {
      auto [log_likelihood_, next_please_, save_kmer] =
          compare_kmers(vlmc_kmer, kmc_kmer, 0);

      if (save_kmer) {
        next_level_kmers.push_back(kmc_kmer);
      }

      next_please = next_please_;
      log_likelihood += log_likelihood_;

      if (next_please) {
        auto status = kmer_database.ReadNextKmer(kmer_api, counter);
        if (!status) {
          return;
        }
        kmc_kmer = kmer_api.construct_vlmc_kmer();
        kmc_kmer.count = counter;
        sequence_length += counter;
      }
    } while (next_please);
  });

  return {log_likelihood, next_level_kmers, sequence_length};
}

double score(KmerContainer<KMerComparator<max_k>> &container,
             const std::string &kmc_db_name, const int actual_kmer_size) {

  auto [full_log_likelihood, next_level_kmers, sequence_length] =
      score_kmers_full_length(container, kmc_db_name, actual_kmer_size);

  double log_likelihood = full_log_likelihood;

  int kmer_length = actual_kmer_size - 1;
  int prefix_length = 1;

  while (!next_level_kmers.empty()) {
    auto [log_likelihood_, next_level_kmers_] =
        score_kmers(container, next_level_kmers, kmer_length, prefix_length);

    next_level_kmers = next_level_kmers_;

    log_likelihood += log_likelihood_;

    kmer_length--;
    prefix_length++;
  }

  return -log_likelihood / double(sequence_length);
}

double negative_log_likelihood_from_kmc_db(
    const std::filesystem::path &fasta_path,
    const std::filesystem::path &tmp_path,
    const std::filesystem::path &vlmc_path, const Core &in_or_out_of_core,
    const int actual_kmer_size, const std::filesystem::path &kmc_db_path) {
  auto kmer_container =
      parse_kmer_container<KMerComparator<max_k>>(in_or_out_of_core);
  load_kmers(vlmc_path, *kmer_container);

  kmer_container->sort();

  double score_ = score(*kmer_container, kmc_db_path, actual_kmer_size);
  std::cout << "nll: " << score_ << std::endl;
  return score_;
}

std::tuple<ankerl::unordered_dense::map<std::string, std::array<double, 4>>,
           int32>
load_kmers_map(const std::filesystem::path &vlmc_path,
               double pseudo_count_amount) {
  int32 longest_kmer = 0;

  ankerl::unordered_dense::map<std::string, std::array<double, 4>> kmers{};
  vlmc::iterate_archive(vlmc_path, [&](const VLMCKmer &kmer) {
    kmers[kmer.to_string()] =
        get_next_symbol_probabilities(kmer, pseudo_count_amount);
    auto probs = get_next_symbol_probabilities(kmer, pseudo_count_amount);
    longest_kmer = std::max(longest_kmer, int(kmer.length));
  });

  return {kmers, longest_kmer};
}

std::array<double, 4> get_closest_state(
    ankerl::unordered_dense::map<std::string, std::array<double, 4>> &kmers,
    const std::string &subsequence) {

  for (int i = 0; i < subsequence.length(); i++) {
    auto iter = kmers.find(subsequence.substr(i));
    if (iter != kmers.end()) {
      return iter->second;
    }
  }

  // This has to exist, and technically, we should already have found it.
  // If it does not exist, something else is wrong.
  return kmers[""];
}

double score_kmer(
    ankerl::unordered_dense::map<std::string, std::array<double, 4>> &kmers,
    const std::string &subsequence, const char &c) {
  auto kmer_probs = get_closest_state(kmers, subsequence);
  int c_i = code_codes.at(c);
  return std::log(kmer_probs[c_i]);
}

std::tuple<std::vector<std::string>, std::vector<double>> iterate_fasta(
    const std::filesystem::path &fasta_path,
    ankerl::unordered_dense::map<std::string, std::array<double, 4>> &kmers,
    const int max_kmer_length) {
  std::vector<std::string> identifiers{};
  std::vector<double> nlls{};

  std::string subsequence{};
  subsequence.resize(max_kmer_length);

  double log_likelihood = 0.0;
  double length = 0.0;

  auto process_line = [&](std::string &line) -> void {
    if (line.substr(0, 1) == ">") {
      identifiers.emplace_back(line.substr(1));
      if (identifiers.size() > 1) {
        nlls.push_back(-log_likelihood / length);
        log_likelihood = 0.0;
        length = 0.0;
      }
    } else {
      for (char c : line) {
        if (code_codes.contains(c)) {
          log_likelihood += score_kmer(kmers, subsequence, c);
          length += 1;
        }

        // Fill the context from the back. First left shift all positions
        // I'm sure this could be optimized
        // When fully adopting c++20
        // `std::shift_left(begin(subsequence), end(subsequence) - 1, 1);`
        // is faster.
        for (int j = 0; j < max_kmer_length - 1; j++) {
          subsequence[j] = subsequence[j + 1];
        }
        // Then fill the last position.
        if (subsequence.length() > 0) {
          subsequence[subsequence.size() - 1] = c;
        }
      }
    }
  };

  if (fasta_path.extension() == ".gz") {
    std::ifstream file_stream(fasta_path,
                              std::ios_base::in | std::ios_base::binary);
    try {
      boost::iostreams::filtering_istream in;
      in.push(boost::iostreams::gzip_decompressor());
      in.push(file_stream);
      for (std::string line; std::getline(in, line);) {
        process_line(line);
      }
    } catch (const boost::iostreams::gzip_error &e) {
      std::cout << e.what() << '\n';
    }
  } else {
    std::ifstream file_stream(fasta_path);
    if (!file_stream.is_open()) {
      throw std::invalid_argument{"Failed to open file."};
    }
    for (std::string line; std::getline(file_stream, line, '\n');) {
      process_line(line);
    }
  }

  nlls.push_back(-log_likelihood / length);
  return {identifiers, nlls};
}

std::tuple<std::vector<std::string>, std::vector<double>> score_file(
    const std::filesystem::path &fasta_path,
    ankerl::unordered_dense::map<std::string, std::array<double, 4>> &kmers,
    const int max_kmer_length) {
  if (fasta_path.extension() == ".fa" || fasta_path.extension() == ".fasta" ||
      fasta_path.extension() == ".fna" || fasta_path.extension() == ".gz") {
    return iterate_fasta(fasta_path, kmers, max_kmer_length);
  } else {
    throw std::invalid_argument("In fasta path is not .fa, .fasta or .fna.");
  }
}
} // namespace vlmc::details

namespace vlmc {
/**
 * Calculation of the negative log-likelihood of all entries in a fasta path
 * for the given vlmc. Iterates the sequence and scores linearily.
 *
 * The outputted negative log likelihood is normalized by the length of the
 * sequence.
 *
 * @param fasta_path Path to fasta (or multifasta) file.
 * @param vlmc_path Path to .bintree vlmc file.
 * @param pseudo_count_amount Pseudo count amount,
 * @return tuple with identifiers for all entries in the fasta and their
 * corresponding normalized negative log likelihoods.
 */
std::tuple<std::vector<std::string>, std::vector<double>>
negative_log_likelihood(const std::filesystem::path &fasta_path,
                        const std::filesystem::path &vlmc_path,
                        const double pseudo_count_amount) {

  auto [kmers, longest_kmer] =
      details::load_kmers_map(vlmc_path, pseudo_count_amount);

  auto [identifiers, nlls] =
      details::score_file(fasta_path, kmers, longest_kmer);

  //  for (int i = 0; i < identifiers.size(); i++) {
  //    std::cout << identifiers[i] << "\t" << nlls[i] << std::endl;
  //  }

  return {identifiers, nlls};
}

/**
 * Calculation of the negative log-likelihood of all entries of all fasta paths
 * and all vlmcs in path. Iterates the sequence and scores linearily.
 *
 * The outputted negative log likelihood is normalized by the length of the
 * sequence.
 *
 * @param fasta_path Path to fasta (or multifasta) file.
 * @param vlmc_path Path to .bintree vlmc file.
 * @param pseudo_count_amount Pseudo count amount,
 * @return tuple with identifiers for all entries in the fasta and their
 * corresponding normalized negative log likelihoods.
 */
std::tuple<std::vector<std::vector<std::string>>,
           std::vector<std::vector<double>>>
negative_log_likelihood_multiple(const std::filesystem::path &fasta_dir_path,
                                 const std::filesystem::path &vlmc_dir_path,
                                 const double pseudo_count_amount,
                                 const int number_of_cores) {

  auto fasta_paths =
      get_recursive_paths(fasta_dir_path, {".fasta", ".fa", ".fna", ".gz"});
  auto vlmc_paths = get_recursive_paths(vlmc_dir_path, {".bintree"});

  std::vector<std::vector<std::string>> identifiers{};
  std::vector<std::vector<double>> nlls{};
  nlls.resize(vlmc_paths.size());
  identifiers.resize(vlmc_paths.size());

  auto fun = [&](size_t start_index, size_t stop_index) {
    for (size_t index = start_index; index < stop_index; index++) {
      std::vector<std::string> identifiers_{};
      std::vector<double> nlls_{};
      for (auto &fasta_path : fasta_paths) {

        auto [tmp_identifiers, tmp_nlls] = negative_log_likelihood(
            fasta_path, vlmc_paths[index], pseudo_count_amount);
        identifiers_.insert(identifiers_.end(), tmp_identifiers.begin(),
                            tmp_identifiers.end());
        nlls_.insert(nlls_.end(), tmp_nlls.begin(), tmp_nlls.end());
      }
      identifiers[index] = identifiers_;
      nlls[index] = nlls_;
    }
  };

  parallel::parallelize(vlmc_paths.size(), fun, number_of_cores);

  for (const auto &identifier : identifiers[0]) {
    std::cout << identifier << "\t";
  }
  std::cout << std::endl;

  for (int i = 0; i < identifiers.size(); i++) {
    std::cout << vlmc_paths[i].stem() << "\t";
    for (int j = 0; j < identifiers[0].size(); j++) {
      std::cout << nlls[i][j] << "\t";
    }
    std::cout << std::endl;
  }

  return {identifiers, nlls};
}

/**
 * Calculate the negative log-likelihood of the given fasta path using the vlmc.
 * This version first computes all kmers with kmc3 and then matches the
 * corresponding kmers to those in the vlmc. This is likely slower than the
 * `negative_log_likelihood(...)` version.
 *
 * @param fasta_path Path to fasta file.
 * @param tmp_path Path to tmp directory for kmc.
 * @param vlmc_path Path to .bintree vlmc file.
 * @param in_or_out_of_core In RAM or on external memory.
 * @param actual_kmer_size Size of the k-mers in the vlmc.
 * @return The combined negative log-likelihood of the fasta.
 */
double negative_log_likelihood_kmc(const std::filesystem::path &fasta_path,
                                   const std::filesystem::path &tmp_path,
                                   const std::filesystem::path &vlmc_path,
                                   const Core &in_or_out_of_core,
                                   const int actual_kmer_size) {

  auto kmc_db_name =
      run_kmc(fasta_path, actual_kmer_size + 1, tmp_path, in_or_out_of_core);

  return details::negative_log_likelihood_from_kmc_db(
      fasta_path, tmp_path, vlmc_path, in_or_out_of_core, actual_kmer_size,
      kmc_db_name);
}
} // namespace vlmc
