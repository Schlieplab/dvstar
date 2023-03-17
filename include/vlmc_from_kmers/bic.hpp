#pragma once

#include "build_vlmc.hpp"
#include "negative_log_likelihood.hpp"

namespace vlmc {

struct Result {
  int min_count;
  int max_depth;
  double threshold;
  double bic_score;
  double log_likelihood;
  uint64 n_parameters;

  void output_result() const {
    std::cout << "Min count: " << this->min_count
              << ". Max depth: " << this->max_depth
              << ". Threshold: " << this->threshold
              << ". Bic score: " << this->bic_score
              << ". Log likelihood: " << this->log_likelihood
              << ". N parameters: " << this->n_parameters << std::endl;
  }
};

std::tuple<uint64, uint64>
count_terminal_nodes(const std::filesystem::path &tree_path) {
  uint64 n_terminal = 0;
  uint64 sequence_size = 0;
  std::ifstream file_stream(tree_path, std::ios::binary);
  {
    cereal::BinaryInputArchive iarchive(file_stream);

    VLMCKmer kmer{};
    while (file_stream.peek() != EOF) {
      iarchive(kmer);

      if (kmer.is_terminal) {
        n_terminal++;
      }

      if (kmer.length == 0) {
        sequence_size = kmer.count;
      }
    }
  }

  file_stream.close();

  return {n_terminal, sequence_size};
}


std::tuple<uint64, uint64>
terminal_node_sum(const std::filesystem::path &tree_path) {
  uint64 terminal_sum = 0;
  uint64 sequence_size = 0;
  std::ifstream file_stream(tree_path, std::ios::binary);
  {
    cereal::BinaryInputArchive iarchive(file_stream);

    VLMCKmer kmer{};
    while (file_stream.peek() != EOF) {
      iarchive(kmer);

      if (kmer.is_terminal) {
        terminal_sum += kmer.length;
      }

      if (kmer.length == 1) {
        sequence_size += kmer.count;
      }
    }
  }

  file_stream.close();

  return {terminal_sum, sequence_size};
}

std::tuple<int, int, double> find_best_parameters_bic(
    const std::filesystem::path &fasta_path, const int max_max_depth,
    const int min_min_count, const std::filesystem::path &tree_path,
    const std::filesystem::path &tmp_path, const Core &in_or_out_of_core) {

  std::vector<Result> results{};

  Result best_result{0, 0, 0.0, std::numeric_limits<double>::max(), 0.0, 0};

  std::vector<double> thresholds{1.2, 3.9075, 10.0};
  std::vector<int> min_counts{
      1,   2,   3,   4,   5,   6,    7,    8,    9,    10,   20,   30,  40,
      50,  60,  70,  80,  90,  100,  125,  150,  175,  200,  250,  300, 400,
      500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 10000};

  for (int max_depth = 3; max_depth < max_max_depth; max_depth++) {
    const int kmer_size = max_depth + 1;

    auto kmc_db_path =
        run_kmc(fasta_path, kmer_size, tmp_path, in_or_out_of_core);

    for (auto min_count : min_counts) {
      for (auto threshold : thresholds) {
        vlmc::build_vlmc_from_kmc_db(kmc_db_path, max_depth, min_count,
                                     threshold, tree_path,
                                     in_or_out_of_core);

        auto [n_terminal, sequence_size] = count_terminal_nodes(tree_path);

        auto nll_start = std::chrono::steady_clock::now();

        double negative_log_likelihood =
            vlmc::negative_log_likelihood_from_kmc_db(
                fasta_path, tmp_path, tree_path, in_or_out_of_core, max_depth,
                kmc_db_path);
        auto nll_done = std::chrono::steady_clock::now();

        std::chrono::duration<double> kmc_seconds = nll_done - nll_start;
        std::cout << "NLL time: " << kmc_seconds.count() << "s\n";

        double log_likelihood = -negative_log_likelihood * sequence_size;

        uint64 n_parameters = 3 * n_terminal;

        double bic_score =
            n_parameters * std::log(sequence_size) - 2 * log_likelihood;

        if (bic_score < best_result.bic_score) {
          best_result = Result{min_count, max_depth,      threshold,
                               bic_score, log_likelihood, n_parameters};
        }

        results.push_back(Result{min_count, max_depth, threshold, bic_score,
                                 log_likelihood, n_parameters});

        results.back().output_result();
      }
    }
  }

  std::stringstream ss;
  ss << "bic-" << fasta_path.stem() << ".csv";
  std::string s = ss.str();

  std::filesystem::path path{s};

  bool add_header = !std::filesystem::exists(path);

  std::ofstream ofs(path, std::ios::app);

  if (add_header) {
    ofs << "n_params,bic_score,min_count,max_depth,threshold,log_likelihood,"
           "fifth_percentile_frequency"
        << std::endl;
  }

  for (auto [min_count, max_depth, threshold, score, log_likelihood, n_params] :
       results) {
    ofs << n_params << "," << score << "," << min_count << "," << max_depth
        << "," << threshold << "," << log_likelihood << "," << 0 << std::endl;
  }
  ofs.close();

  std::cout << "---------- Best bic ----------" << std::endl;
  best_result.output_result();

  return {best_result.max_depth, best_result.min_count, best_result.threshold};
}

} // namespace vlmc
