#include "vlmc_from_kmers/build_vlmc.hpp"
#include "vlmc_from_kmers/bic.hpp"
#include "vlmc_from_kmers/dvstar.hpp"

int main(int argc, char *argv[]) {
  CLI::App app{"Construction and comparisons of variable-length Markov chains "
               "with the aid of "
               "a k-mer counter."};

  vlmc::cli_arguments arguments{};
  add_options(app, arguments);

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  bool tmp_path_existed_before = std::filesystem::exists(arguments.tmp_path);

  if (arguments.in_or_out_of_core == vlmc::Core::out) {
    std::filesystem::create_directories(arguments.tmp_path);
  }
  vlmc::configure_stxxl(arguments.tmp_path);

  if (arguments.mode == vlmc::Mode::build) {
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

  } else if (arguments.mode == vlmc::Mode::build_from_kmc_db) {
    int exit_code = vlmc::build_vlmc_from_kmc_db(
        arguments.in_path, arguments.max_depth, arguments.min_count,
        arguments.threshold, arguments.out_path, arguments.in_or_out_of_core,
        arguments.pseudo_count_amount, arguments.estimator,
        arguments.sequencing_parameters);

    return exit_code;
  } else if (arguments.mode == vlmc::Mode::dump) {
    return vlmc::dump_path(arguments.in_path, arguments.out_path);
  } else if (arguments.mode == vlmc::Mode::score_sequence) {
    vlmc::negative_log_likelihood(
        arguments.fasta_path, arguments.tmp_path, arguments.in_path,
        arguments.in_or_out_of_core, arguments.max_depth);
  } else if (arguments.mode == vlmc::Mode::bic) {
    vlmc::find_best_parameters_bic(
        arguments.fasta_path, arguments.max_depth, arguments.min_count,
        arguments.out_path, arguments.tmp_path, arguments.in_or_out_of_core);
  } else if (arguments.mode == vlmc::Mode::dissimilarity) {
    double dissimilarity;
    if (arguments.dissimilarity == vlmc::Dissimilarity::dvstar_dissimliarity) {
      dissimilarity = vlmc::dvstar(arguments.in_path, arguments.to_path, 0,
                                   arguments.pseudo_count_amount);
    } else if (arguments.dissimilarity ==
               vlmc::Dissimilarity::penalized_dvstar_dissimliarity) {
      dissimilarity =
          vlmc::dvstar_missing_penalized(arguments.in_path, arguments.to_path,
                                         0, arguments.pseudo_count_amount);
    }
    std::cout << dissimilarity << std::endl;

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
