#include <gtest/gtest.h>

#include "vlmc_from_kmers/build_vlmc.hpp"
#include "vlmc_from_kmers/negative_log_likelihood.hpp"

#include "read_helper.hpp"

using namespace vlmc;

class NegativeLogLikelihoodTests : public ::testing::Test {
protected:
  void SetUp() override {}
  std::string binary_seq{"AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAATAAAAAAA"
                         "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAATAAAAAAA"
                         "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAATAAAAAAA"
                         "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAATAAAAAAA"
                         "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAATAAAAAAA"
                         "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAATAAAAAAA"
                         "AAAATAAAAATAAAAATAAAAAAATAAAAAAATAAAAAAAAAAAAAAAAAAAA"
                         "AAAAAAATAAAAAAA"};
};

void output_fasta(const std::filesystem::path &path,
                  const std::string &contents) {
  std::ofstream out{path};
  out << ">TEST\n";
  for (int i = 0; i < contents.length() / 80 + 1; i++) {
    auto subseq = contents.substr(i * 80, 80);
    out << subseq << "\n";
  }
  out.close();
}

void generate_fasta(const std::filesystem::path &path) {
  std::ofstream file_stream{path};
  file_stream << ">Test" << std::endl;
  file_stream << "ACGTACGATCGTACGATACGTACGATCGTACGATACGTACGATCGT" << std::endl;
  file_stream << "CCGTACGATCGTACGAGACGTACGATCGTACGATCCGTACTATCGT" << std::endl;
  file_stream << "TTCAGCTAGTCAGATCTAGCTAGCTAGCTACGTACTAGCTAGCGAT" << std::endl;
  file_stream << "ACGATCGTAGTCGATCAGTCATCTACTGACACGAGCGCTCGATGAC" << std::endl;
  file_stream << "ATCATCATCATAATTCAGTATATAGCATTATCGATTATCGATTACG" << std::endl;
  file_stream << "ATCGATCGGATCGATCGATCGATCGATCGATCGATCGATCGATCGA" << std::endl;
  file_stream.close();
}

TEST_F(NegativeLogLikelihoodTests, Correct) {
  std::filesystem::path temp_fasta_1{"tmp1.fasta"};
  std::filesystem::path temp_fasta_2{"tmp2.fasta"};
  std::filesystem::path temp{"tmp"};
  std::filesystem::path out{"_test.tree"};

  generate_fasta(temp_fasta_1);
  generate_fasta(temp_fasta_2);

  build_vlmc(temp_fasta_1, 2, 1, 3.9075, out, temp, Core::in);

  double nll = vlmc::negative_log_likelihood_kmc(temp_fasta_2, temp, out,
                                                 vlmc::Core::in, 2);

  EXPECT_NEAR(0.9023370264265379, nll, 0.00001);
}

TEST_F(NegativeLogLikelihoodTests, GetClosest) {
  ankerl::unordered_dense::map<std::string, std::array<double, 4>> kmers{};
  kmers["GT"] = {1.0, 1.0, 1.0, 1.0};
  kmers[""] = {0.7, 0.5, 0.0, 0.0};
  auto closest_prob = vlmc::details::get_closest_state(kmers, "ACGT");

  ASSERT_FLOAT_EQ(closest_prob[0], 1.0);
  ASSERT_FLOAT_EQ(closest_prob[1], 1.0);
  ASSERT_FLOAT_EQ(closest_prob[2], 1.0);
  ASSERT_FLOAT_EQ(closest_prob[3], 1.0);
}

TEST_F(NegativeLogLikelihoodTests, LogLikelihoodHandCrafted0Order) {
  ankerl::unordered_dense::map<std::string, std::array<double, 4>> kmers{};
  kmers[""] = {0.8, 0.0, 0.0, 0.2};

  auto test_path = std::filesystem::temp_directory_path() / "test.fasta";
  output_fasta(test_path, binary_seq);

  auto [identifiers, nlls] = vlmc::details::score_file(test_path, kmers, 0);
  std::filesystem::remove(test_path);

  double negative_log_likelihood_manual =
      -(std::log(0.8) * 434 + std::log(0.2) * 42);
  negative_log_likelihood_manual /= binary_seq.length();

  EXPECT_FLOAT_EQ(nlls[0], negative_log_likelihood_manual);
}
TEST_F(NegativeLogLikelihoodTests, LogLikelighoodHandCrafted1Order) {
  ankerl::unordered_dense::map<std::string, std::array<double, 4>> kmers{};
  kmers[""] = {0.8, 0.0, 0.0, 0.2};
  kmers["A"] = {0.9, 0.0, 0.0, 0.1};
  kmers["T"] = {1.0, 0.0, 0.0, 0.0};

  auto test_path = std::filesystem::temp_directory_path() / "test.fasta";
  output_fasta(test_path, binary_seq);

  auto [identifiers, nlls] = vlmc::details::score_file(test_path, kmers, 1);
  std::filesystem::remove(test_path);

  double negative_log_likelihood_manual = -(391 * std::log(0.9) +
                                          42 * std::log(0.1) +
                                          42 * std::log(1) + std::log(0.8));
  negative_log_likelihood_manual /= binary_seq.length();

  EXPECT_FLOAT_EQ(nlls[0], negative_log_likelihood_manual);
}