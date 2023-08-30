#include <benchmark/benchmark.h>
#include <vector>

#include "vlmc_from_kmers/dvstar.hpp"
#include "vlmc_from_kmers/distances/distances/dvstar.hpp"

class DistanceBenchmarks : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State &state) override {}

  std::filesystem::path first_bintree_path{"NC_028367.bintree"};
  std::filesystem::path second_bintree_path{"NC_045512.bintree"};
  std::filesystem::path third_bintree_path{"NC_001497.bintree"};
};

BENCHMARK_F(DistanceBenchmarks, OriginalPath)
(benchmark::State &state) {
  for (auto _ : state) {
    double path_dist = vlmc::dvstar(first_bintree_path, second_bintree_path, 0);
  }
}

BENCHMARK_F(DistanceBenchmarks, OriginalParsed)
(benchmark::State &state) {
  auto one_kmers = vlmc::get_sorted_kmers(first_bintree_path);
  auto two_kmers = vlmc::get_sorted_kmers(second_bintree_path);

  for (auto _ : state) {
    double kmers_dist = dvstar(one_kmers, two_kmers, 0);
  }
}

BENCHMARK_F(DistanceBenchmarks, SortedVector)
(benchmark::State &state) {

  auto left = vlmc::container::SortedVector{first_bintree_path, 0};
  auto right = vlmc::container::SortedVector{first_bintree_path, 0};
  for (auto _ : state) {
    vlmc::distance::dvstar<vlmc::container::SortedVector>(left, right);
  }
}

BENCHMARK_F(DistanceBenchmarks, BTree)
(benchmark::State &state) {

  auto left = vlmc::container::BTree{first_bintree_path, 0};
  auto right = vlmc::container::BTree{first_bintree_path, 0};
  for (auto _ : state) {
    vlmc::distance::dvstar<vlmc::container::BTree>(left, right);
  }
}

BENCHMARK_F(DistanceBenchmarks, HashMap)
(benchmark::State &state) {

  auto left = vlmc::container::HashMap{first_bintree_path, 0};
  auto right = vlmc::container::HashMap{first_bintree_path, 0};
  for (auto _ : state) {
    vlmc::distance::dvstar<vlmc::container::HashMap>(left, right);
  }
}

BENCHMARK_F(DistanceBenchmarks, EytzingerTree)
(benchmark::State &state) {

  auto left = vlmc::container::EytzingerTree{first_bintree_path, 0};
  auto right = vlmc::container::EytzingerTree{first_bintree_path, 0};
  for (auto _ : state) {
    vlmc::distance::dvstar<vlmc::container::EytzingerTree>(left, right);
  }
}

BENCHMARK_F(DistanceBenchmarks, SortedSearch)
(benchmark::State &state) {

  auto left = vlmc::container::SortedSearch{first_bintree_path, 0};
  auto right = vlmc::container::SortedSearch{first_bintree_path, 0};
  for (auto _ : state) {
    vlmc::distance::dvstar<vlmc::container::SortedSearch>(left, right);
  }
}

BENCHMARK_F(DistanceBenchmarks, VanEmdeBoasTree)
(benchmark::State &state) {

  auto left = vlmc::container::VanEmdeBoasTree{first_bintree_path, 0};
  auto right = vlmc::container::VanEmdeBoasTree{first_bintree_path, 0};
  for (auto _ : state) {
    vlmc::distance::dvstar<vlmc::container::VanEmdeBoasTree>(left, right);
  }
}

BENCHMARK_MAIN();
