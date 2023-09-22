#include <benchmark/benchmark.h>

#include <algorithm>

class NLLBenchmarks : public benchmark::Fixture {};

BENCHMARK_DEFINE_F(NLLBenchmarks, ShiftStringLoop)
(benchmark::State &state) {
  std::string subsequence{};
  subsequence.resize(state.range(0));
  for (auto _ : state) {
    for (int j = 0; j < state.range(0) - 1; j++) {
      subsequence[j] = subsequence[j + 1];
    }
  }
}

BENCHMARK_REGISTER_F(NLLBenchmarks, ShiftStringLoop)->DenseRange(3, 20, 2);

BENCHMARK_DEFINE_F(NLLBenchmarks, ShiftStringSubstr)
(benchmark::State &state) {
  int length = state.range(0);
  std::string subsequence{};
  subsequence.resize(length);
  for (auto _ : state) {
    subsequence = subsequence.substr(1, length - 1) + "A";
  }
}

BENCHMARK_REGISTER_F(NLLBenchmarks, ShiftStringSubstr)->DenseRange(3, 20, 2);


BENCHMARK_DEFINE_F(NLLBenchmarks, ShiftStringStdShiftLeft)
(benchmark::State &state) {
  int length = state.range(0);
  std::string subsequence{};
  subsequence.resize(length);
  for (auto _ : state) {
    std::shift_left(begin(subsequence), end(subsequence), 1);
  }
}

BENCHMARK_REGISTER_F(NLLBenchmarks, ShiftStringStdShiftLeft)->DenseRange(3, 20, 2);


BENCHMARK_MAIN();