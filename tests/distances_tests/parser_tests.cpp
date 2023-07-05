#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>

#include "vlmc_from_kmers/distances/vlmc_container.hpp"
#include "vlmc_from_kmers/distances/cluster_container.hpp"
#include "vlmc_from_kmers/distances/parser.hpp"

using vlmc_c = vlmc::container::SortedVector;

class ParserTests : public ::testing::Test {
protected:
  void SetUp() override {}
};

// Tests when comparing two directories 
TEST_F(ParserTests, paths) {
  CLI::App app{"Distance comparison of either one directory or between two different directories."};
  vlmc::parser::cli_arguments arguments{};
  add_options(app, arguments);

  int argc = 5; 
  char * args1[] = {"./dist", "-p", "./test1", "-s", "./test2"};
  char ** argv1 = args1; 
  
  try {
    app.parse(argc, argv1);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  EXPECT_TRUE(arguments.first_VLMC_path == std::filesystem::path("./test1"));
  EXPECT_TRUE(arguments.second_VLMC_path == std::filesystem::path("./test2")); 

  char * args2[] = {"./dist", "--VLMC-path", "./test3", "--snd-VLMC-path", "./test4"};
  char ** argv2 = args2; 
  
  try {
    app.parse(argc, argv2);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  EXPECT_TRUE(arguments.first_VLMC_path == std::filesystem::path("./test3"));
  EXPECT_TRUE(arguments.second_VLMC_path == std::filesystem::path("./test4")); 
}
