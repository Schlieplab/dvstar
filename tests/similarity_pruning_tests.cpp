#include <gtest/gtest.h>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_from_kmers/similarity_pruning.hpp"

#include <cereal/archives/binary.hpp>
#include <kmc_file.h>

#include "read_helper.hpp"

using namespace vlmc;

class SimilarityPruningTests : public ::testing::Test {
protected:
  void SetUp() override {}

  KMersPerLevel<PstKmer> kmers_per_level{4, 8};
};

TEST_F(SimilarityPruningTests, KullbackLiebler) {
  auto parent = create_kmer("TTTTT");
  parent.count = 6;
  parent.next_symbol_counts = {3, 1, 1, 1};

  auto child = create_kmer("TTTTTA");
  child.count = 3;
  child.next_symbol_counts = {0, 1, 1, 1};

  std::vector<double> parent_probs{4.0 / 10, 2.0 / 10, 2.0 / 10, 2.0 / 10};
  std::vector<double> child_probs{1.0 / 7, 2.0 / 7, 2.0 / 7, 2.0 / 7};

  double expected_kl =
      (child_probs[0] * std::log(child_probs[0] / parent_probs[0])) +
      (child_probs[1] * std::log(child_probs[1] / parent_probs[1])) +
      (child_probs[2] * std::log(child_probs[2] / parent_probs[2])) +
      (child_probs[3] * std::log(child_probs[3] / parent_probs[3]));
  expected_kl *= 3.0;

  double real_kl = kl_divergence(child, parent, 1.0);

  EXPECT_EQ(expected_kl, real_kl);
}

TEST_F(SimilarityPruningTests, SimilarityPruneSameLevel) {
  auto prev_kmer = create_kmer("TTTTTT");
  auto kmer = create_kmer("ATTTTT");

  std::ofstream file_stream("test_tmp.bin", std::ios::binary);
  cereal::BinaryOutputArchive oarchive(file_stream);

  estimator_f keep_node = kl_estimator(3.9075);
  similarity_prune(prev_kmer, kmer, kmers_per_level, oarchive, keep_node);

  ASSERT_TRUE(kmers_per_level[6][3].real_child);
  EXPECT_EQ(kmers_per_level[6][3].kmer.to_string(), "ATTTTT");
}

TEST_F(SimilarityPruningTests, SimilarityPruneParent) {
  auto a_kmer = create_kmer("ATTTTT");
  a_kmer.count = 7;
  a_kmer.next_symbol_counts = {4, 1, 1, 1};

  auto c_kmer = create_kmer("CTTTTT");
  c_kmer.count = 3;
  c_kmer.next_symbol_counts = {1, 2, 0, 0};

  auto g_kmer = create_kmer("GTTTTT");
  g_kmer.count = 4;
  g_kmer.next_symbol_counts = {2, 1, 0, 1};

  auto t_kmer = create_kmer("TTTTTT");
  t_kmer.count = 1;
  t_kmer.next_symbol_counts = {0, 1, 0, 0};
  auto kmer = create_kmer("TTTTT");

  kmer.count = 15;
  kmer.next_symbol_counts = {7, 3, 4, 1};

  kmers_per_level[6][0] = {a_kmer, false, true};
  kmers_per_level[6][1] = {c_kmer, false, true};
  kmers_per_level[6][2] = {g_kmer, false, true};
  kmers_per_level[6][3] = {t_kmer, true, true};

  std::filesystem::path tmp_path{"test_tmp.bintree"};
  std::stringstream out;

  {
    std::ofstream file_stream(tmp_path, std::ios::binary);
    cereal::BinaryOutputArchive oarchive(file_stream);

    estimator_f keep_node = kl_estimator(3.9075);
    auto [has_children, is_terminal] =
        process_parent(kmer, kmers_per_level, oarchive, keep_node);

    EXPECT_TRUE(has_children);
  }

  EXPECT_EQ(kmers_per_level[6][0].real_child, false);
  EXPECT_EQ(kmers_per_level[6][1].real_child, false);
  EXPECT_EQ(kmers_per_level[6][2].real_child, false);
  EXPECT_EQ(kmers_per_level[6][3].real_child, false);

  {

    std::ifstream file_stream(tmp_path, std::ios::binary);
    cereal::BinaryInputArchive iarchive(file_stream);

    VLMCKmer kmer{};

    while (file_stream.peek() != EOF) {
      try {
        iarchive(kmer);

        EXPECT_EQ(kmer.to_string(), "TTTTTT");
        EXPECT_EQ(kmer.count, t_kmer.count);
        EXPECT_EQ(kmer.next_symbol_counts[0], 0);
        EXPECT_EQ(kmer.next_symbol_counts[1], 1);
        EXPECT_EQ(kmer.next_symbol_counts[2], 0);
        EXPECT_EQ(kmer.next_symbol_counts[3], 0);
      } catch (const cereal::Exception &e) {
        EXPECT_TRUE(file_stream.peek() == EOF);
      }
    }
  }
}
