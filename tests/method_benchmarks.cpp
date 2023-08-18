#include <benchmark/benchmark.h>
#include <random>
#include <vector>

#include <kmc_file.h>

#include "vlmc_from_kmers/build_vlmc.hpp"
#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_from_kmers/support_pruning.hpp"
#include "read_helper.hpp"

class MethodBenchmarks : public benchmark::Fixture {
  static void config_stxxl() {
    static bool have_configured = false;
    static std::string tmp_path = "test_tmp";
    if (!have_configured) {
      have_configured = true;
      std::filesystem::create_directories(tmp_path);
      configure_stxxl(tmp_path);
    }
  }

public:
  void SetUp(const ::benchmark::State &state) override { config_stxxl(); }
  std::filesystem::path fasta_path{"NC_022098.1.fasta"};
};

BENCHMARK_F(MethodBenchmarks, InCore)
(benchmark::State &state) {

  std::filesystem::path tmp_path = "test_tmp";

  std::filesystem::create_directories(tmp_path);
  std::filesystem::path run_one_path{"out-one.tree"};
  vlmc::configure_stxxl(tmp_path);

  for (auto _ : state) {
    int exit_one_code = vlmc::build_vlmc(fasta_path, 9, 2, 3.9075,
                                         run_one_path, tmp_path, Core::in);
  }
}

BENCHMARK_F(MethodBenchmarks, OutOfCore)
(benchmark::State &state) {
  std::filesystem::path tmp_path = "test_tmp";

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);
  std::filesystem::path run_path{"out-one.tree"};

  for (auto _ : state) {
    int exit_one_code = vlmc::build_vlmc(fasta_path, 9, 2, 3.9075, run_path,
                                         tmp_path, Core::out);
  }
}

BENCHMARK_F(MethodBenchmarks, HashMap)
(benchmark::State &state) {

  std::filesystem::path tmp_path = "test_tmp";

  std::filesystem::create_directories(tmp_path);
  std::filesystem::path run_path{"out-one.tree"};
  vlmc::configure_stxxl(tmp_path);

  for (auto _ : state) {
    int exit_one_code = vlmc::build_vlmc(fasta_path, 9, 2, 3.9075, run_path,
                                         tmp_path, Core::hash);
  }
}

BENCHMARK_MAIN();
