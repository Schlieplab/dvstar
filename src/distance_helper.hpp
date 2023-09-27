#include <highfive/H5Easy.hpp>

#include "vlmc_from_kmers/build_vlmc.hpp"
#include "vlmc_from_kmers/cli_helper.hpp"
#include "vlmc_from_kmers/distances/calc_dists.hpp"
#include "vlmc_from_kmers/distances/get_cluster.hpp"
#include "vlmc_from_kmers/distances/global_aliases.hpp"
#include "vlmc_from_kmers/distances/utils.hpp"

void output_square_format(std::ostream &stream, const matrix_t &distance_matrix,
                          const std::vector<std::string> &ids_from,
                          const std::vector<std::string> &ids_to) {
  for (size_t j = 0; j < distance_matrix.cols(); j++) {
    stream << ids_to[j] << "\t";
  }
  stream << std::endl;

  for (size_t i = 0; i < distance_matrix.rows(); i++) {
    stream << ids_from[i] << "\t";
    for (size_t j = 0; j < distance_matrix.cols(); j++) {
      stream << distance_matrix(i, j) << "\t";
    }
    stream << std::endl;
  }
}

void output_phylip_format(std::ostream &stream, const matrix_t &distance_matrix,
                          const std::vector<std::string> &ids_from) {
  stream << ids_from.size() << std::endl;

  for (size_t i = 0; i < distance_matrix.rows(); i++) {
    stream << ids_from[i] << "\t";
    for (size_t j = 0; j < distance_matrix.cols(); j++) {
      stream << distance_matrix(i, j) << "\t";
    }
    stream << std::endl;
  }
}

void output_distances(const vlmc::cli_arguments &arguments,
                      const matrix_t &distance_matrix,
                      const std::vector<std::string> &ids_from,
                      const std::vector<std::string> &ids_to) {
  auto format = arguments.distances_output_format;
  if (arguments.out_path.extension() == ".h5" ||
      arguments.out_path.extension() == ".hdf5") {
    format = vlmc::DistancesFormat::hdf5;
  }

  if (arguments.out_path.empty() && format == vlmc::DistancesFormat::hdf5) {
    std::cerr << "Need to provide path in '--out-path' to write hdf5 files. "
                 "Will output to stdout in square format instead.";
    format = vlmc::DistancesFormat::square;
  }

  if (format == vlmc::DistancesFormat::hdf5) {
    H5Easy::File file{arguments.out_path, HighFive::File::OpenOrCreate};

    H5Easy::dump<Eigen::MatrixXd>(file, "/distances", distance_matrix,
                                  H5Easy::DumpMode::Overwrite);

    H5Easy::dump(file, "/ids", ids_from, H5Easy::DumpMode::Overwrite);
    H5Easy::dump(file, "/ids_to", ids_to, H5Easy::DumpMode::Overwrite);

    std::clog << "Wrote distances to: " << arguments.out_path.string()
              << std::endl;
  } else if (arguments.out_path.empty()) {
    if (format == vlmc::DistancesFormat::phylip) {
      output_phylip_format(std::cout, distance_matrix, ids_from);
    } else {
      output_square_format(std::cout, distance_matrix, ids_from, ids_to);
    }
  } else {
    auto out_path = arguments.out_path;
    if (std::filesystem::is_directory(out_path)) {
      out_path = out_path / "distances.dist";
    }
    if (!out_path.has_extension()) {
      out_path.replace_extension(".dist");
    }

    std::ofstream stream{out_path};

    if (format == vlmc::DistancesFormat::phylip) {
      output_phylip_format(stream, distance_matrix, ids_from);
    } else {
      output_square_format(stream, distance_matrix, ids_from, ids_to);
    }
    std::clog << "Wrote distances to: " << out_path.string() << std::endl;
  }
}

template <typename VC>
std::tuple<matrix_t, std::vector<std::string>, std::vector<std::string>>
calculate_cluster_distance(const std::filesystem::path &in_path,
                           const double pseudo_count_amount,
                           int background_order, const size_t n_threads) {
  auto [cluster, ids] = vlmc::get_cluster<VC>(
      in_path, n_threads, background_order, pseudo_count_amount);

  std::clog << "Calculating distances matrix of size " << cluster.size() << "x"
            << cluster.size() << std::endl;
  return {vlmc::calc_dist::calculate_distances<VC>(cluster, n_threads), ids,
          ids};
}

template <typename VC>
std::tuple<matrix_t, std::vector<std::string>, std::vector<std::string>>
calculate_cluster_distance(const std::filesystem::path &in_path,
                           const std::filesystem::path &to_path,
                           const double pseudo_count_amount,
                           int background_order, const size_t n_threads) {
  auto [cluster, ids] = vlmc::get_cluster<VC>(
      in_path, n_threads, background_order, pseudo_count_amount);

  auto [cluster_to, ids_to] = vlmc::get_cluster<VC>(
      to_path, n_threads, background_order, pseudo_count_amount);

  std::clog << "Calculating distances matrix of size " << cluster.size() << "x"
            << cluster_to.size() << std::endl;

  return {
      vlmc::calc_dist::calculate_distances<VC>(cluster, cluster_to, n_threads),
      ids, ids_to};
}

template <typename VC>
std::tuple<matrix_t, std::vector<std::string>, std::vector<std::string>>
calculate_cluster_distance(const vlmc::cli_arguments &arguments,
                           const size_t n_threads) {
  if (arguments.to_path.empty()) {
    return calculate_cluster_distance<VC>(arguments.in_path,
                                          arguments.pseudo_count_amount,
                                          arguments.background_order, n_threads);

  } else {
    return calculate_cluster_distance<VC>(arguments.in_path, arguments.to_path,
                                          arguments.pseudo_count_amount,
                                          arguments.background_order, n_threads);
  }
}

std::tuple<matrix_t, std::vector<std::string>, std::vector<std::string>>
apply_container(const vlmc::cli_arguments &arguments,
                vlmc::VLMCRepresentation vlmc_container,
                const size_t n_threads) {
  if (vlmc_container == vlmc::VLMCRepresentation::vlmc_sorted_vector) {
    return calculate_cluster_distance<vlmc::container::SortedVector>(arguments,
                                                                     n_threads);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_b_tree) {
    return calculate_cluster_distance<vlmc::container::BTree>(arguments,
                                                              n_threads);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_hashmap) {
    return calculate_cluster_distance<vlmc::container::HashMap>(arguments,
                                                                n_threads);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_veb) {
    return calculate_cluster_distance<vlmc::container::VanEmdeBoasTree>(
        arguments, n_threads);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_ey) {
    return calculate_cluster_distance<vlmc::container::EytzingerTree>(arguments,
                                                                      n_threads);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_sorted_search) {
    return calculate_cluster_distance<vlmc::container::SortedSearch>(arguments,
                                                                     n_threads);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_kmer_major) {
    throw std::logic_error{"K-mer major currently not implemented."};
    //    return calculate_kmer_major(arguments, n_threads);
  }
  return {};
}

int compute_dissimilarity(vlmc::cli_arguments &arguments) {
  if (arguments.in_path.empty()) {
    std::cerr << "Error: A input path to .bintree files has to be given for "
                 "comparison operation."
              << std::endl;
    return EXIT_FAILURE;
  }

  size_t n_threads =
      vlmc::parse_degree_of_parallelism(arguments.degree_of_parallelism);

  auto [distance_matrix, ids_from, ids_to] =
      apply_container(arguments, arguments.vlmc_representation, n_threads);

  output_distances(arguments, distance_matrix, ids_from, ids_to);

  return EXIT_SUCCESS;
}

int compute_dissimilarity_fasta(vlmc::cli_arguments &arguments) {
  if (arguments.fasta_path.empty()) {
    std::cerr << "Error: A --fasta-path needs to be provided." << std::endl;
    return EXIT_FAILURE;
  }
  if (arguments.cache_path.empty()) {
    std::cerr << "Error: An --cache-path needs to be provided." << std::endl;
    return EXIT_FAILURE;
  }

  auto fasta_paths = vlmc::get_recursive_paths(
      arguments.fasta_path, {".fasta", ".fa", ".fna", ".gz"});

  auto bintree_path = arguments.cache_path;

  bool tmp_path_existed_before = std::filesystem::exists(arguments.tmp_path);

  vlmc::configure_stxxl(arguments.tmp_path);

  auto start_building = std::chrono::steady_clock::now();

  auto n_threads =
      vlmc::parse_degree_of_parallelism(arguments.degree_of_parallelism);

  for (auto &fasta_path : fasta_paths) {
    auto out_path = bintree_path / fasta_path.stem();
    out_path.replace_extension(".bintree");

    if (!std::filesystem::exists(out_path)) {
      int exit_code = vlmc::build_vlmc(
          fasta_path, arguments.max_depth, arguments.min_count,
          arguments.threshold, out_path, arguments.tmp_path,
          arguments.in_or_out_of_core, arguments.pseudo_count_amount,
          arguments.estimator, arguments.sequencing_parameters, true);

      if (exit_code != EXIT_SUCCESS) {
        return exit_code;
      }
    }
  }
  auto done_building = std::chrono::steady_clock::now();
  std::chrono::duration<double> building_duration =
      done_building - start_building;

  std::clog << "VLMC construction time: " << building_duration.count() << "s\n";

  auto start_distances = std::chrono::steady_clock::now();

  auto [distance_matrix, ids_from, ids_to] =
      calculate_cluster_distance<vlmc::container::SortedSearch>(
          bintree_path, arguments.pseudo_count_amount,
          arguments.background_order, n_threads);

  auto done_distances = std::chrono::steady_clock::now();
  std::chrono::duration<double> distances_duration =
      done_distances - start_distances;
  std::clog << "VLMC distances time: " << distances_duration.count() << "s\n";

  output_distances(arguments, distance_matrix, ids_from, ids_to);

  if (!tmp_path_existed_before) {
    std::filesystem::remove_all(arguments.tmp_path);
  }

  return EXIT_SUCCESS;
}
