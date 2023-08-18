#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include <kmc_file.h>

#include "vlmc_from_kmers/build_vlmc.hpp"
#include "vlmc_from_kmers/dvstar.hpp"

using namespace vlmc;

class DvstarTests : public ::testing::Test {
protected:
  void SetUp() override {}

  std::filesystem::path first_fasta{"NC_028367.1.fa"};
  std::filesystem::path second_fasta{"NC_045512.2.fa"};
  std::filesystem::path third_fasta{"NC_001497.2.fa"};
};


TEST_F(DvstarTests, Identity) {
  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);

  std::filesystem::path run_one_path{"out-one.bintree"};
  int exit_one_code = vlmc::build_vlmc(first_fasta, 6, 2, 3.9075,
                                       run_one_path, tmp_path, Core::in);

  std::filesystem::path run_two_path{"out-two.bintree"};
  int exit_two_code = vlmc::build_vlmc(first_fasta, 6, 2, 3.9075,
                                       run_two_path, tmp_path, Core::in);

  double dist = dvstar(run_one_path, run_two_path, 0);

  EXPECT_FLOAT_EQ(dist, 0.0);
}

TEST_F(DvstarTests, Similar) {
  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);

  std::filesystem::path run_one_path{"NC_028367.bintree"};
  int exit_one_code = vlmc::build_vlmc(first_fasta, 6, 2, 3.9075,
                                       run_one_path, tmp_path, Core::in);

  std::filesystem::path run_two_path{"NC_045512.bintree"};
  int exit_two_code = vlmc::build_vlmc(second_fasta, 6, 2, 3.9075,
                                       run_two_path, tmp_path, Core::in);

  double dist = dvstar(run_one_path, run_two_path, 0);

  EXPECT_FLOAT_EQ(dist, 0.15452552);
}

TEST_F(DvstarTests, SimilarPseudo) {
  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);

  std::filesystem::path run_one_path{"NC_028367.bintree"};
  int exit_one_code = vlmc::build_vlmc(first_fasta, 6, 2, 3.9075,
                                       run_one_path, tmp_path, Core::in);

  std::filesystem::path run_two_path{"NC_045512.bintree"};
  int exit_two_code = vlmc::build_vlmc(second_fasta, 6, 2, 3.9075,
                                       run_two_path, tmp_path, Core::in);

  double dist = dvstar(run_one_path, run_two_path, 0);

  EXPECT_FLOAT_EQ(dist, 0.15452552);
}

TEST_F(DvstarTests, Different) {
  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);
  std::cout << std::filesystem::current_path() << std::endl;

  std::filesystem::path run_one_path{"NC_028367.bintree"};
  int exit_one_code = vlmc::build_vlmc(first_fasta, 6, 2, 1.2,
                                       run_one_path, tmp_path, Core::in);

  std::filesystem::path run_two_path{"NC_001497.bintree"};
  int exit_two_code = vlmc::build_vlmc(third_fasta, 6, 2, 1.2,
                                       run_two_path, tmp_path, Core::in);

  double dist = dvstar(run_one_path, run_two_path, 0);

  EXPECT_FLOAT_EQ(dist, 0.27162704);
}

TEST_F(DvstarTests, RegressionVersions) {
  std::filesystem::path tmp_path = std::filesystem::temp_directory_path();

  std::filesystem::create_directories(tmp_path);
  vlmc::configure_stxxl(tmp_path);
  std::cout << std::filesystem::current_path() << std::endl;

  std::filesystem::path run_one_path{"NC_028367.bintree"};
  int exit_one_code = vlmc::build_vlmc(first_fasta, 6, 2, 1.2,
                                       run_one_path, tmp_path, Core::in);

  std::filesystem::path run_two_path{"NC_001497.bintree"};
  int exit_two_code = vlmc::build_vlmc(third_fasta, 6, 2, 1.2,
                                       run_two_path, tmp_path, Core::in);

  auto one_kmers = vlmc::get_sorted_kmers(run_one_path);
  auto two_kmers = vlmc::get_sorted_kmers(run_two_path);

  double path_dist = dvstar(run_one_path, run_two_path, 0);
  double kmers_dist = dvstar(one_kmers, two_kmers, 0);

  EXPECT_FLOAT_EQ(path_dist, kmers_dist);
}