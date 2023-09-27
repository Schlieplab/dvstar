#include "highfive/H5Easy.hpp"

#include "distance_helper.hpp"
#include "vlmc_from_kmers/bic.hpp"
#include "vlmc_from_kmers/build_vlmc.hpp"
#include "vlmc_from_kmers/dvstar.hpp"

int build_from_kmc_db(const vlmc::cli_arguments &arguments) {
  int exit_code = vlmc::build_vlmc_from_kmc_db(
      arguments.in_path, arguments.max_depth, arguments.min_count,
      arguments.threshold, arguments.out_path, arguments.in_or_out_of_core,
      arguments.pseudo_count_amount, arguments.estimator,
      arguments.sequencing_parameters);

  return exit_code;
}
int build(const vlmc::cli_arguments &arguments, bool tmp_path_existed_before) {
  if (arguments.out_path.empty() || arguments.fasta_path.empty()) {
    std::cerr
        << "Error: Both a --fasta-path and an --out-path need to be given."
        << std::endl;
    return EXIT_FAILURE;
  }

  auto out_path = arguments.out_path;
  if (!out_path.has_extension()) {
    out_path.replace_extension(".bintree");
  }

  int exit_code = vlmc::build_vlmc(
      arguments.fasta_path, arguments.max_depth, arguments.min_count,
      arguments.threshold, out_path, arguments.tmp_path,
      arguments.in_or_out_of_core, arguments.pseudo_count_amount,
      arguments.estimator, arguments.sequencing_parameters);

  if (!tmp_path_existed_before) {
    std::filesystem::remove_all(arguments.tmp_path);
  }

  return exit_code;
}

int main(int argc, char *argv[]) {
  CLI::App app{"Construction and comparisons of variable-length Markov chains "
               "with the aid of a k-mer counter."};

  vlmc::cli_arguments arguments{};

  auto build_sub = app.add_subcommand(
      "build", "Constructs a single VLMC from the provided fasta file.");
  vlmc::add_build_options(*build_sub, arguments);

  bool tmp_path_existed_before = std::filesystem::exists(arguments.tmp_path);

  build_sub->callback([&]() {
    if (arguments.in_or_out_of_core == vlmc::Core::out) {
      std::filesystem::create_directories(arguments.tmp_path);
    }

    vlmc::configure_stxxl(arguments.tmp_path);
    return build(arguments, tmp_path_existed_before);
  });

  auto build_from_kmc_sub =
      app.add_subcommand("build-from-kmc-db",
                         "Constructs a single VLMC from the provided kmc db.");
  vlmc::add_build_from_kmc_options(*build_from_kmc_sub, arguments);

  build_from_kmc_sub->callback([&]() {
    if (arguments.in_or_out_of_core == vlmc::Core::out) {
      std::filesystem::create_directories(arguments.tmp_path);
    }

    vlmc::configure_stxxl(arguments.tmp_path);
    return build_from_kmc_db(arguments);
  });

  auto dump_sub =
      app.add_subcommand("dump", "Dumps the contents of a .bintree file to raw "
                                 "text, either to a .txt file or stdout.");
  vlmc::add_dump_options(*dump_sub, arguments);

  dump_sub->callback(
      [&]() { return vlmc::dump_path(arguments.in_path, arguments.out_path); });

  auto score_sub =
      app.add_subcommand("score", "Computes the negative log-likelihood of a "
                                  "collection of fasta files for VLMCs.");
  vlmc::add_score_options(*score_sub, arguments);

  score_sub->callback([&]() {
    int nr_cores =
        vlmc::parse_degree_of_parallelism(arguments.degree_of_parallelism);

    vlmc::negative_log_likelihood_multiple(
        arguments.fasta_path, arguments.in_path, arguments.pseudo_count_amount,
        nr_cores);
  });

  auto bic_sub = app.add_subcommand(
      "bic", "Runs BIC and AIC to give insight into parameter choice of the "
             "VLMC for the given fasta file.");
  vlmc::add_bic_options(*bic_sub, arguments);

  bic_sub->callback([&]() {
    if (arguments.in_or_out_of_core == vlmc::Core::out) {
      std::filesystem::create_directories(arguments.tmp_path);
    }
    vlmc::configure_stxxl(arguments.tmp_path);

    vlmc::find_best_parameters_bic(
        arguments.fasta_path, arguments.max_depth, arguments.min_count,
        arguments.out_path, arguments.tmp_path, arguments.in_or_out_of_core);
  });

  auto dissimilarity_sub = app.add_subcommand(
      "dissimilarity",
      "Computes the dissimilarity between a collection of VLMCs.");
  vlmc::add_distance_options(*dissimilarity_sub, arguments);

  dissimilarity_sub->callback(
      [&]() { return compute_dissimilarity(arguments); });

  auto dissimilarity_fasta_sub = app.add_subcommand(
      "dissimilarity-fasta", "Computes VLMCs from the fasta files and the "
                             "dissimilarity between the resulting VLMCs.");
  vlmc::add_distance_fasta_options(*dissimilarity_fasta_sub, arguments);

  dissimilarity_fasta_sub->callback(
      [&]() { return compute_dissimilarity_fasta(arguments); });

  auto reprune_sub = app.add_subcommand(
      "reprune", "Runs the pruning steps of the VLMC construction to make a "
                 "given VLMC more general.");

  vlmc::add_reprune_options(*reprune_sub, arguments);
  reprune_sub->callback([&]() {
    return vlmc::reprune_vlmc(arguments.in_path, arguments.out_path,
                              arguments.in_or_out_of_core, arguments.threshold,
                              arguments.pseudo_count_amount);
  });

  auto size_sub = app.add_subcommand(
      "size", "Prints the size of the VLMC. The size is determined by the sum "
              "of the length of all (terminal) k-mers in the VLMCs. A k-mer is "
              "considered terminal if any branch ends with the k-mer, meaning "
              "there is no more specific k-mer in the VLMC.");
  vlmc::add_size_options(*size_sub, arguments);
  size_sub->callback([&]() {
    auto [terminal_size, sequence_size] =
        vlmc::terminal_node_sum(arguments.in_path);

    std::cout << "Terminal node sum: " << terminal_size
              << ", Sequence size: " << sequence_size << std::endl;
  });

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  // if (!tmp_path_existed_before) {
  //   std::filesystem::remove_all(arguments.tmp_path);
  // }

  return EXIT_SUCCESS;
}
