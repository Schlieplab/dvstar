#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include "vlmc_from_kmers/cli_helper.hpp"
#include "vlmc_from_kmers/distances/cluster_container.hpp"
#include "vlmc_from_kmers/distances/vlmc_container.hpp"

using vlmc_c = vlmc::container::SortedVector;

class ParserTests : public ::testing::Test {
protected:
  void SetUp() override {}
};

// Tests when comparing two directories
TEST_F(ParserTests, paths) {
  CLI::App app{"Distance comparison of either one directory or between two "
               "different directories."};
  vlmc::cli_arguments arguments{};
  add_options(app, arguments);

  int argc = 5;
  char *args1[] = {"./dist", "-p", "./test1", "--to-path", "./test2"};
  char **argv1 = args1;

  try {
    app.parse(argc, argv1);
  } catch (const CLI::ParseError &e) {
    FAIL();
  }

  EXPECT_TRUE(arguments.fasta_path == std::filesystem::path("./test1"));
  EXPECT_TRUE(arguments.to_path == std::filesystem::path("./test2"));

  char *args2[] = {"./dist", "--in-path", "./test3", "--to-path",
                   "./test4"};
  char **argv2 = args2;

  try {
    app.parse(argc, argv2);
  } catch (const CLI::ParseError &e) {

    FAIL();
  }

  EXPECT_TRUE(arguments.in_path == std::filesystem::path("./test3"));
  EXPECT_TRUE(arguments.to_path == std::filesystem::path("./test4"));
}
