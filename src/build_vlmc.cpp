#include "highfive/H5Easy.hpp"

#include "vlmc_from_kmers/build_vlmc.hpp"
#include "distance_helper.hpp"
#include "vlmc_from_kmers/bic.hpp"
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
  int exit_code = vlmc::build_vlmc(
      arguments.fasta_path, arguments.max_depth, arguments.min_count,
      arguments.threshold, arguments.out_path, arguments.tmp_path,
      arguments.in_or_out_of_core, arguments.pseudo_count_amount,
      arguments.estimator, arguments.sequencing_parameters);

  if (!tmp_path_existed_before) {
    std::filesystem::remove_all(arguments.tmp_path);
  }

  return exit_code;
}

int main(int argc, char *argv[]) {
  CLI::App app{"Construction and comparisons of variable-length Markov chains "
               "with the aid of "
               "a k-mer counter."};

  vlmc::cli_arguments arguments{};
  add_options(app, arguments);
  add_distance_options(app, arguments);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  bool tmp_path_existed_before = std::filesystem::exists(arguments.tmp_path);

  if (arguments.in_or_out_of_core == vlmc::Core::out) {
    std::filesystem::create_directories(arguments.tmp_path);
  }

  if (arguments.mode == vlmc::Mode::build) {
    vlmc::configure_stxxl(arguments.tmp_path);
    return build(arguments, tmp_path_existed_before);

  } else if (arguments.mode == vlmc::Mode::build_from_kmc_db) {
    vlmc::configure_stxxl(arguments.tmp_path);
    return build_from_kmc_db(arguments);

  } else if (arguments.mode == vlmc::Mode::dump) {
    return vlmc::dump_path(arguments.in_path, arguments.out_path);
  } else if (arguments.mode == vlmc::Mode::score_sequence) {
    int nr_cores =
        vlmc::parse_degree_of_parallelism(arguments.degree_of_parallelism);
    vlmc::negative_log_likelihood_multiple(
        arguments.fasta_path, arguments.in_path,
        arguments.pseudo_count_amount, nr_cores);

  } else if (arguments.mode == vlmc::Mode::bic) {
    vlmc::find_best_parameters_bic(
        arguments.fasta_path, arguments.max_depth, arguments.min_count,
        arguments.out_path, arguments.tmp_path, arguments.in_or_out_of_core);

  } else if (arguments.mode == vlmc::Mode::dissimilarity) {
    compute_dissimilarity(arguments);

  } else if (arguments.mode == vlmc::Mode::reprune) {
    return vlmc::reprune_vlmc(arguments.in_path, arguments.out_path,
                              arguments.in_or_out_of_core, arguments.threshold,
                              arguments.pseudo_count_amount);

  } else if (arguments.mode == vlmc::Mode::size) {
    auto [terminal_size, sequence_size] =
        vlmc::terminal_node_sum(arguments.in_path);
    std::cout << "Terminal node sum: " << terminal_size
              << ", Sequence size: " << sequence_size << std::endl;
  }

  //if (!tmp_path_existed_before) {
  //  std::filesystem::remove_all(arguments.tmp_path);
  //}

  return EXIT_SUCCESS;
}
