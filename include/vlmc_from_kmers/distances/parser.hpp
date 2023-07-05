#pragma once

#include <filesystem>

#include "CLI/App.hpp"
#include "CLI/Config.hpp"
#include "CLI/Formatter.hpp"

#include "global_aliases.hpp"
#include "vlmc_container.hpp"

namespace vlmc::parser {
enum VLMCRep {
  vlmc_sorted_search,
  vlmc_sorted_vector,
  vlmc_b_tree,
  vlmc_ey,
  vlmc_hashmap,
  vlmc_kmer_major,
  vlmc_veb
};

struct cli_arguments {
  std::filesystem::path first_VLMC_path{};
  std::filesystem::path second_VLMC_path{};
  std::filesystem::path out_path{};
  size_t dop{1};
  int set_size{-1};
  VLMCRep vlmc{VLMCRep::vlmc_sorted_search};
  size_t background_order{0};
};

size_t parse_dop(size_t requested_cores) {
  if (requested_cores < 1) {
    throw std::invalid_argument("Too low degree of parallelism, must be >= 1");
  } else {
    return requested_cores;
  }
}

void add_options(CLI::App &app, cli_arguments &arguments) {
  std::map<std::string, VLMCRep> VLMC_Rep_map{
      {"sbs", VLMCRep::vlmc_sorted_search},
      {"sorted-vector", VLMCRep::vlmc_sorted_vector},
      {"b-tree", VLMCRep::vlmc_b_tree},
      {"eytzinger", VLMCRep::vlmc_ey},
      {"HashMap", VLMCRep::vlmc_hashmap},
      {"kmer-major", VLMCRep::vlmc_kmer_major},
      {"veb", VLMCRep::vlmc_veb}};

  app.add_option("-p,--VLMC-path", arguments.first_VLMC_path,
                 "Required for distance calculation. 'Primary' path to saved "
                 "bintree directory. If '-s' is empty it will compute the "
                 "inter-distance between the trees of the directory.");

  app.add_option("-s,--snd-VLMC-path", arguments.second_VLMC_path,
                 "Optional 'Secondary' path to saved bintree directory. "
                 "Calculates distance between the trees specified in -p "
                 "(primary) and -s (secondary).");

  app.add_option("-o,--matrix-path", arguments.out_path,
                 "Path to hdf5 file where scores will be stored.");

  app.add_option("-n,--max-dop", arguments.dop,
                 "Degree of parallelism. Default 1 (sequential).");

  app.add_option("-v,--vlmc-rep", arguments.vlmc,
                 "Vlmc container representation to use.")
      ->transform(CLI::CheckedTransformer(VLMC_Rep_map, CLI::ignore_case));

  app.add_option("-b,--background-order", arguments.background_order,
                 "Background order.");

  app.add_option("-a, --set-size", arguments.set_size,
                 "Number of VLMCs to compute distance function on.");
}
} // namespace vlmc::parser
