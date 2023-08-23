#include <highfive/H5Easy.hpp>

#include "vlmc_from_kmers/cli_helper.hpp"
#include "vlmc_from_kmers/distances/calc_dists.hpp"
#include "vlmc_from_kmers/distances/get_cluster.hpp"
#include "vlmc_from_kmers/distances/global_aliases.hpp"
#include "vlmc_from_kmers/distances/parser.hpp"
#include "vlmc_from_kmers/distances/utils.hpp"

matrix_t calculate_kmer_major(const vlmc::cli_arguments &arguments,
                              const size_t nr_cores) {
  size_t use_cores = nr_cores;
  size_t max_cores = std::thread::hardware_concurrency();
  if (max_cores < nr_cores) {
    use_cores = max_cores;
  }

  auto cluster = vlmc::get_kmer_cluster(arguments.in_path, use_cores,
                                        arguments.background_order);
  if (arguments.to_path.empty()) {
    std::cout << "Calculating distances for single cluster." << std::endl;
    return vlmc::calc_dist::calculate_distance_major(cluster, cluster,
                                                     use_cores);
  }
  auto cluster_to = vlmc::get_kmer_cluster(arguments.to_path, use_cores,
                                           arguments.background_order);
  std::cout << "Calculating distances." << std::endl;
  return vlmc::calc_dist::calculate_distance_major(cluster, cluster_to,
                                                   use_cores);
}

template <typename VC>
std::tuple<matrix_t, std::vector<std::string>, std::vector<std::string>>
calculate_cluster_distance(const vlmc::cli_arguments &arguments,
                           const size_t nr_cores) {
  auto [cluster, ids] = vlmc::get_cluster<VC>(arguments.in_path, nr_cores,
                                              arguments.background_order);
  if (arguments.to_path.empty()) {
    std::clog << "Calculating distances matrix of size " << cluster.size()
              << "x" << cluster.size() << std::endl;
    return {vlmc::calc_dist::calculate_distances<VC>(cluster, nr_cores), ids,
            ids};
  }
  auto [cluster_to, ids_to] = vlmc::get_cluster<VC>(arguments.to_path, nr_cores,
                                                    arguments.background_order);

  std::clog << "Calculating distances matrix of size " << cluster.size() << "x"
            << cluster_to.size() << std::endl;

  return {
      vlmc::calc_dist::calculate_distances<VC>(cluster, cluster_to, nr_cores),
      ids, ids_to};
}

std::tuple<matrix_t, std::vector<std::string>, std::vector<std::string>>
apply_container(const vlmc::cli_arguments &arguments,
                vlmc::VLMCRepresentation vlmc_container,
                const size_t nr_cores) {
  if (vlmc_container == vlmc::VLMCRepresentation::vlmc_sorted_vector) {
    return calculate_cluster_distance<vlmc::container::SortedVector>(arguments,
                                                                     nr_cores);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_b_tree) {
    return calculate_cluster_distance<vlmc::container::BTree>(arguments,
                                                              nr_cores);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_hashmap) {
    return calculate_cluster_distance<vlmc::container::HashMap>(arguments,
                                                                nr_cores);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_veb) {
    return calculate_cluster_distance<vlmc::container::VanEmdeBoasTree>(
        arguments, nr_cores);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_ey) {
    return calculate_cluster_distance<vlmc::container::EytzingerTree>(arguments,
                                                                      nr_cores);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_sorted_search) {
    return calculate_cluster_distance<vlmc::container::SortedSearch>(arguments,
                                                                     nr_cores);
  } else if (vlmc_container == vlmc::VLMCRepresentation::vlmc_kmer_major) {
    throw std::logic_error{"K-mer major currently not implemented."};
    //    return calculate_kmer_major(arguments, nr_cores);
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

  size_t nr_cores =
      vlmc::parse_degree_of_parallelism(arguments.degree_of_parallelism);

  auto [distance_matrix, ids_from, ids_to] =
      apply_container(arguments, arguments.vlmc_representation, nr_cores);

  if (arguments.out_path.empty()) {
    utils::print_matrix(distance_matrix, ids_from, ids_to);
  } else if (arguments.out_path.extension() == ".h5" ||
             arguments.out_path.extension() == ".hdf5") {
    H5Easy::File file{arguments.out_path, HighFive::File::OpenOrCreate};

    H5Easy::dump<Eigen::MatrixXd>(file, "/distances", distance_matrix,
                                  H5Easy::DumpMode::Overwrite);

    H5Easy::dump(file, "/ids", ids_from, H5Easy::DumpMode::Overwrite);
    H5Easy::dump(file, "/ids_to", ids_to, H5Easy::DumpMode::Overwrite);

    std::clog << "Wrote distances to: " << arguments.out_path.string()
              << std::endl;
  }

  return EXIT_SUCCESS;
}