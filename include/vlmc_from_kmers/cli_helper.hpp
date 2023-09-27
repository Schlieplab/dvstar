#pragma once

#include <filesystem>
#include <memory>
#include <random>
#include <string>

#include "kmer_container.hpp"
#include "sequencing_adjustment.hpp"

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

namespace vlmc {

enum Mode {
  build,
  score_sequence,
  dump,
  bic,
  build_from_kmc_db,
  dissimilarity,
  reprune,
  size,
  dissimilarity_fasta,
};

enum class DistancesFormat {
  phylip,
  square,
  hdf5,
};

enum VLMCRepresentation {
  vlmc_sorted_search,
  vlmc_sorted_vector,
  vlmc_b_tree,
  vlmc_ey,
  vlmc_hashmap,
  vlmc_kmer_major,
  vlmc_veb
};

enum Dissimilarity {
  dvstar_dissimliarity,
};

enum Estimator { kullback_leibler, peres_shields };

struct cli_arguments {
  Mode mode{Mode::build};
  Dissimilarity dissimilarity{Dissimilarity::dvstar_dissimliarity};
  Estimator estimator{Estimator::kullback_leibler};
  std::filesystem::path fasta_path;
  std::filesystem::path in_path;
  std::filesystem::path to_path;
  std::filesystem::path tmp_path{"./tmp"};
  std::filesystem::path out_path{};
  std::filesystem::path cache_path{"./vlmc"};

  int degree_of_parallelism{
      static_cast<int>(std::thread::hardware_concurrency())};
  VLMCRepresentation vlmc_representation{
      VLMCRepresentation::vlmc_sorted_search};

  int min_count = 2;
  int max_depth = 9;
  double threshold = 3.9075;

  double pseudo_count_amount = 0.001;
  int background_order = 0;

  Core in_or_out_of_core{Core::in};

  SequencingParameters sequencing_parameters{false};

  DistancesFormat distances_output_format{DistancesFormat::phylip};
};

void add_pseudo_count_option(CLI::App &app, cli_arguments &arguments) {
  app.add_option(
      "-a,--pseudo-count-amount", arguments.pseudo_count_amount,
      "Size of pseudo count for probability estimation. See e.g. "
      "https://en.wikipedia.org/wiki/Additive_smoothing . Defaults to 0.001.");
}

void add_memory_model_options(CLI::App &app, cli_arguments &arguments) {
  std::map<std::string, Core> core_map{
      {"internal", Core::in}, {"external", Core::out}, {"hash", Core::hash}};
  app.add_option(
         "-i, --in-or-out-of-core", arguments.in_or_out_of_core,
         "Specify 'internal' for in-core or 'external' for out-of-core memory "
         "model.  Out of core is slower, but is not memory bound. Defaults to "
         "internal.")
      ->transform(CLI::CheckedTransformer(core_map, CLI::ignore_case));
}

void add_tmp_path_option(CLI::App &app, cli_arguments &arguments) {
  app.add_option("-t,--temp-path", arguments.tmp_path,
                 "Path to temporary folder for the external memory algorithms. "
                 " For good performance, this needs to be on a local machine.  "
                 "For sorting, at least 2GB will be allocated to this path.  "
                 "Defaults to ./tmp");
}

void add_shared_build_options(CLI::App &app, cli_arguments &arguments) {
  std::map<std::string, Estimator> estimator_map{
      {"kullback-leibler", Estimator::kullback_leibler},
      {"peres-shields", Estimator::peres_shields},
  };

  app.add_option("--estimator", arguments.estimator,
                 "Estimator for the pruning of the VLMC, either "
                 "'kullback-leibler',  or 'peres-shields'.")
      ->transform(CLI::CheckedTransformer(estimator_map, CLI::ignore_case));

  add_tmp_path_option(app, arguments);

  app.add_option(
      "-c,--min-count", arguments.min_count,
      "Minimum count required for every k-mer in the tree. Defaults to 2.");

  app.add_option("-k,--threshold", arguments.threshold,
                 "Kullback-Leibler threshold. Defaults to 3.9075.");

  app.add_option("-d,--max-depth", arguments.max_depth,
                 "Maximum depth/length for included k-mers. Defaults to 9.");

  add_pseudo_count_option(app, arguments);

  app.add_flag(
      "--adjust-for-sequencing-errors",
      arguments.sequencing_parameters.adjust_for_sequencing_errors,
      "Give this flag to adjust the estimator parameters and min counts for "
      "the sequencing depth and error rates of a sequencing dataset. See "
      "--sequencing-depth and --sequencing-error-rate for parameters.");

  app.add_option("--sequencing-depth", arguments.sequencing_parameters.depth,
                 "If --adjust-for-sequencing-errors is given, this parameter "
                 "is used to alter the estimator parameters to reflect that "
                 "many k-mers will be --sequencing-depth times more frequent.");

  app.add_option("--sequencing-error-rate",
                 arguments.sequencing_parameters.error_rate,
                 "If --adjust-for-sequencing-errors is given, this parameter "
                 "is used to alter to estimate the number of k-mers that will "
                 "be missing due to sequencing errors.");

  add_memory_model_options(app, arguments);
}

void add_build_options(CLI::App &app, cli_arguments &arguments) {
  app.add_option("-p,--fasta-path", arguments.fasta_path, "Path to fasta file.")
      ->required();
  app.add_option("-o,--out-path", arguments.out_path,
                 "Path to output .bintree file. Will add .bintree as extension "
                 "to the file if missing.")
      ->required();

  add_shared_build_options(app, arguments);
}

void add_build_from_kmc_options(CLI::App &app, cli_arguments &arguments) {
  app.add_option("-p,--in-path", arguments.in_path,
                 "Path to kmc db file. The path to the kmc db file needs to be "
                 "supplied without the file extension")
      ->required();
  app.add_option("-o,--out-path", arguments.out_path,
                 "Path to output .bintree file. Will add .bintree as extension "
                 "to the file if missing.")
      ->required();

  add_shared_build_options(app, arguments);
}

void add_dump_options(CLI::App &app, cli_arguments &arguments) {
  app.add_option("-p,--in-path", arguments.in_path,
                 "Path to a .bintree file (the output from `build`).")
      ->required();

  app.add_option(
      "-o,--out-path", arguments.out_path,
      "Path to output text file. Will write to stdout if not provided.");
}

void add_score_options(CLI::App &app, cli_arguments &arguments) {
  app.add_option("--in-path", arguments.in_path,
                 "Path to a .bintree file (the output from `build`) or a "
                 "directory of .bintree files.")
      ->required();

  app.add_option("--fasta-path", arguments.fasta_path,
                 "Path to either a fasta file or a directory of fasta files.")
      ->required();

  add_pseudo_count_option(app, arguments);

  app.add_option("-n,--n-threads", arguments.degree_of_parallelism,
                 "Number of threads for distance computations. Defaults to all "
                 "available cores.");
}

void add_bic_options(CLI::App &app, cli_arguments &arguments) {
  add_memory_model_options(app, arguments);

  app.add_option("-p,--fasta-path", arguments.fasta_path,
                 "Path to fasta file to build VLMC for.")
      ->required();

  app.add_option("-d,--max-max-depth", arguments.max_depth,
                 "The largest max depth to try. Larger max depths drastically "
                 "increases the computation time. Defaults to 9.");

  app.add_option("-t,--min-min-count", arguments.min_count,
                 "The smallest min count to try. Smaller min counts increase "
                 "the computation time. Defaults to 2.");

  app.add_option("-o,--out-path", arguments.out_path,
                 "The .bintree file is written to this location during "
                 "computation of the BIC.")
      ->required();

  add_pseudo_count_option(app, arguments);
}

int parse_degree_of_parallelism(int requested_cores) {
  if (requested_cores < 1) {
    throw std::invalid_argument("Too low degree of parallelism, must be >= 1");
  } else {
    return requested_cores;
  }
}

void add_shared_distance_options(CLI::App &app, cli_arguments &arguments) {
  std::map<std::string, VLMCRepresentation> VLMC_rep_map{
      {"sbs", VLMCRepresentation::vlmc_sorted_search},
      {"sorted-vector", VLMCRepresentation::vlmc_sorted_vector},
      {"b-tree", VLMCRepresentation::vlmc_b_tree},
      {"eytzinger", VLMCRepresentation::vlmc_ey},
      {"HashMap", VLMCRepresentation::vlmc_hashmap},
      {"kmer-major", VLMCRepresentation::vlmc_kmer_major},
      {"veb", VLMCRepresentation::vlmc_veb}};

  app.add_option("-v,--vlmc-rep", arguments.vlmc_representation,
                 "VLMC container representation to use for distances "
                 "computation. Defaults to sorted-search which gives the "
                 "fastest average execution time.")
      ->transform(CLI::CheckedTransformer(VLMC_rep_map, CLI::ignore_case));

  app.add_option("-n,--n-threads", arguments.degree_of_parallelism,
                 "Number of threads for distance computations. Defaults to all "
                 "available cores.");

  app.add_option("-b,--background-order", arguments.background_order,
                 "Background order.");

  std::map<std::string, DistancesFormat> format_map{
      {"phylip", DistancesFormat::phylip},
      {"square", DistancesFormat::square},
      {"hdf5", DistancesFormat::hdf5}};

  app.add_option("--distances-format", arguments.distances_output_format,
                 "Output format of the distances. Either 'phylip', 'square', "
                 "or 'hdf5'.")
      ->transform(CLI::CheckedTransformer(format_map, CLI::ignore_case));

  app.add_option(
      "-o,--out-path", arguments.out_path,
      "Path to output file. If empty, the result will be written to stdout.");

  // No longer have multiple dissimilarities implemented here.
  //  std::map<std::string, Dissimilarity> dissimilarity_map{
  //      {"dvstar", Dissimilarity::dvstar_dissimliarity},
  //  };
  //  app.add_option("--dissimilarity", arguments.dissimilarity,
  //                 "Dissimilarity type, either 'dvstar',  or
  //                 'penalized-dvstar'.")
  //      ->transform(CLI::CheckedTransformer(dissimilarity_map,
  //      CLI::ignore_case));
}

void add_distance_options(CLI::App &app, cli_arguments &arguments) {
  app.add_option(
         "--in-path", arguments.in_path,
         "Path to a saved .bintree file or directory of .bintree files.")
      ->required();

  app.add_option(
      "--to-path", arguments.to_path,
      "Path to a saved .bintree file or directory of .bintree files. "
      "If empty, the distances are "
      "computed between the files in the --in-path.");

  add_shared_distance_options(app, arguments);
  add_pseudo_count_option(app, arguments);
}

void add_distance_fasta_options(CLI::App &app, cli_arguments &arguments) {

  app.add_option("--fasta-path", arguments.fasta_path,
                 "Path to a fasta file or directory of fasta files.")
      ->required();

  app.add_option("--to-path", arguments.to_path,
                 "Path to a fasta file or directory of fasta files. "
                 "If empty, the distances are "
                 "computed between the files in the --in-path.");

  app.add_option(
      "--cache-path", arguments.cache_path,
      "Path to directory where the computed VLMCs are stored. Defaults to "
      "'vlmc'. The resulting "
      "VLMCs are not cleaned after computation has finished to speed-up "
      "further computation with the same VLMCs and parameters. If a VLMC with "
      "a matching name is found here, it will not be re-computed. Note that "
      "this is not parameter-aware and so should be cleaned or changed to a "
      "different path if the parameters are changed.");

  add_shared_distance_options(app, arguments);

  add_shared_build_options(app, arguments);
}

void add_reprune_options(CLI::App &app, cli_arguments &arguments) {
  app.add_option("--in-path", arguments.in_path,
                 "Path to a previously computed .bintree file.")
      ->required();

  app.add_option("-o,--out-path", arguments.out_path,
                 "Path to output .bintree file.");

  add_memory_model_options(app, arguments);

  app.add_option("-k,--threshold", arguments.threshold,
                 "New Kullback-Leibler threshold.");

  add_pseudo_count_option(app, arguments);
}

void add_size_options(CLI::App &app, cli_arguments &arguments) {
  app.add_option("--in-path", arguments.in_path,
                 "Path to a previously computed .bintree file.")
      ->required();
}

static std::random_device rd;

static std::mt19937 gen = std::mt19937{rd()};
std::string get_random_name(const std::string &start) {
  std::stringstream ss;
  ss << start;

  // Generate 20 characters of random lower-case letters
  std::uniform_int_distribution<> distrib(97, 122);
  for (size_t i = 0; i < 20; i++) {
    auto rand = distrib(gen);
    ss << char(rand);
  }

  return ss.str();
}

template <class Comparator>
std::shared_ptr<KmerContainer<Comparator>>
parse_kmer_container(const Core &in_or_out_of_core) {
  if (in_or_out_of_core == Core::out) {
    return std::make_shared<OutOfCoreKmerContainer<Comparator>>();
  } else if (in_or_out_of_core == Core::in) {
    return std::make_shared<InCoreKmerContainer<Comparator>>();
  } else if (in_or_out_of_core == Core::hash) {
    return std::make_shared<HashMapKmerContainer<Comparator>>();
  } else {
    throw(std::invalid_argument(
        "parameter --out-or-in-of-core not 'internal', 'hash' or 'external'"));
  }
}

} // namespace vlmc
