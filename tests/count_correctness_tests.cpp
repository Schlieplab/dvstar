#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include <kmc_file.h>

#include "vlmc_from_kmers/build_vlmc.hpp"
#include "read_helper.hpp"

using namespace vlmc;

// TODO: This test checks for all k-mers in the sequence.
// This is not the same the count of all characters, as this also includes
// the last k "end"-mers, which give plus 1 to a couple of values, and plus
// k to the root node.

class CountCorrectnessTests : public ::testing::Test {
protected:
  void SetUp() override {}
};

std::map<std::string, size_t> get_kmer_counts(std::string &sequence) {
  std::map<std::string, size_t> counts{};

  for (size_t i = 0; i <= sequence.size() - 11; i++) {
    for (int j = 0; j < 12 && j + i <= sequence.size(); j++) {
      auto kmer = sequence.substr(i, j);
      if (counts.find(kmer) != counts.end()) {
        counts[kmer] += 1;
      } else {
        counts[kmer] = 1;
      }
    }
  }

  return counts;
}

TEST_F(CountCorrectnessTests, CorrectCounts) {
  std::filesystem::path fasta_path{"../python-prototype/NC_022098.1.fasta"};
  std::ifstream fasta_stream{fasta_path};

  std::string sequence{};
  for (std::string line; std::getline(fasta_stream, line);) {
    if (line.empty()) {
      continue;
    }
    if (line[0] == '>') {
      continue;
    }

    sequence += line;
  }

  auto counts = get_kmer_counts(sequence);

  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);

  std::filesystem::path run_path{"out.tree"};
  int exit_code =
      vlmc::build_vlmc(fasta_path, 10, 1, 3.9075, run_path, tmp_path, vlmc::Core::in);

  std::ifstream run_stream{run_path};

  std::ifstream file_stream(run_path, std::ios::binary);
  {
    cereal::BinaryInputArchive iarchive(file_stream);

    VLMCKmer kmer{};
    while (file_stream.peek() != EOF) {
      iarchive(kmer);
      auto kmer_string = kmer.to_string();
      EXPECT_EQ(counts[kmer_string], kmer.count) << kmer_string;
      EXPECT_EQ(counts[kmer_string + 'A'], kmer.next_symbol_counts[0])
          << kmer_string + 'A';
      EXPECT_EQ(counts[kmer_string + 'C'], kmer.next_symbol_counts[1])
          << kmer_string + 'C';
      EXPECT_EQ(counts[kmer_string + 'G'], kmer.next_symbol_counts[2])
          << kmer_string + 'G';
      ASSERT_EQ(counts[kmer_string + 'T'], kmer.next_symbol_counts[3])
          << kmer_string + 'T';
    }
  }

  file_stream.close();
}

// Hash-based count did not turn out to be faster here. The test currently does not work.
//TEST_F(CountCorrectnessTests, CorrectCountsHash) {
//  std::filesystem::path fasta_path{"../tests/NC_022098.1.fasta"};
//  std::ifstream fasta_stream{fasta_path};
//
//  std::string sequence{};
//  for (std::string line; std::getline(fasta_stream, line);) {
//    if (line.empty()) {
//      continue;
//    }
//    if (line[0] == '>') {
//      continue;
//    }
//
//    sequence += line;
//  }
//
//  auto counts = get_kmer_counts(sequence);
//
//  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();
//
//  std::filesystem::create_directories(tmp_path);
//  vlmc::configure_stxxl(tmp_path);
//
//  std::filesystem::path run_path{"out.tree"};
//  int exit_code =
//      vlmc::build_vlmc(fasta_path, 10, 1, 3.9075, run_path, tmp_path, vlmc::Core::hash);
//
//  std::ifstream run_stream{run_path};
//
//  std::ifstream file_stream(run_path, std::ios::binary);
//  {
//    cereal::BinaryInputArchive iarchive(file_stream);
//
//    VLMCKmer kmer{};
//    while (file_stream.peek() != EOF) {
//      iarchive(kmer);
//      auto kmer_string = kmer.to_string();
//      EXPECT_EQ(counts[kmer_string], kmer.count) << kmer_string;
//      EXPECT_EQ(counts[kmer_string + 'A'], kmer.next_symbol_counts[0])
//          << kmer_string + 'A';
//      EXPECT_EQ(counts[kmer_string + 'C'], kmer.next_symbol_counts[1])
//          << kmer_string + 'C';
//      EXPECT_EQ(counts[kmer_string + 'G'], kmer.next_symbol_counts[2])
//          << kmer_string + 'G';
//      ASSERT_EQ(counts[kmer_string + 'T'], kmer.next_symbol_counts[3])
//          << kmer_string + 'T';
//    }
//  }
//
//  file_stream.close();
//}