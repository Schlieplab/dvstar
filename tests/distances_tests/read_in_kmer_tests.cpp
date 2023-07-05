#include <gtest/gtest.h>

#include <functional>
#include <cstdlib>
#include <fstream>

#include "vlmc_from_kmers/kmer.hpp"
#include "vlmc_from_kmers/distances/read_in_kmer.hpp"
#include "../read_helper.hpp"
#include "vlmc_from_kmers/distances/global_aliases.hpp"

const out_t error_tolerance = 1E-7;

class RIKmerTest : public ::testing::Test {
protected:
  void SetUp() override {}
};

TEST_F(RIKmerTest, KmerConstructorProb) {
  vlmc::VLMCKmer old_kmer{5, 10, {5, 10, 7, 8}};
  double sum_children = 30.0+4.0;
  
  vlmc::ReadInKmer kmer{old_kmer};
  for (size_t i = 0; i < 4; i++){
    EXPECT_NEAR(kmer.next_char_prob[i], (old_kmer.next_symbol_counts[i] + 1)/sum_children, error_tolerance);
  }
}

TEST_F(RIKmerTest, KmerConstructorIntRep) {
  std::string kmer_string{"A"};
  auto old_kmer = create_kmer(kmer_string);
  
  vlmc::ReadInKmer kmer{old_kmer};
  EXPECT_EQ(kmer.integer_rep, kmer.get_index_rep(old_kmer));
}

TEST_F(RIKmerTest, KmerBackgroundRep1) {
  std::string kmer_string{"A"};
  auto old_kmer = create_kmer(kmer_string);
  vlmc::ReadInKmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 0), 0);
}
TEST_F(RIKmerTest, KmerBackgroundRep2) {
  std::string kmer_string{"A"};
  auto old_kmer = create_kmer(kmer_string);
  vlmc::ReadInKmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 1), 1);
}
TEST_F(RIKmerTest, KmerBackgroundRep3) {
  std::string kmer_string{"ATT"};
  auto old_kmer = create_kmer(kmer_string);
  vlmc::ReadInKmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 1), 4);
}
TEST_F(RIKmerTest, KmerBackgroundRep4) {
  std::string kmer_string{"ATT"};
  auto old_kmer = create_kmer(kmer_string);
  vlmc::ReadInKmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 2), 20);
}
TEST_F(RIKmerTest, KmerBackgroundRep5) {
  std::string kmer_string{"AAG"};
  auto old_kmer = create_kmer(kmer_string);
  vlmc::ReadInKmer kmer{old_kmer};
  EXPECT_EQ(kmer.background_order_index(kmer.integer_rep, 2), 7);
}
void createStrComb(std::vector<std::string> &strings, std::string str, int size){
  if (size <= 0){
    return; 
  }
  std::string str_A = str; 
  std::string str_C = str; 
  std::string str_G = str;
  std::string str_T = str; 

  str_A.append("A");
  strings.push_back(str_A);
  str_C.append("C");
  strings.push_back(str_C);
  str_G.append("G");
  strings.push_back(str_G);
  str_T.append("T");
  strings.push_back(str_T);

  createStrComb(strings, str_A, size - 1);
  createStrComb(strings, str_C, size - 1);
  createStrComb(strings, str_G, size - 1);
  createStrComb(strings, str_T, size - 1);
}

bool compare (std::pair <int, std::string> &a, std::pair <int, std::string> &b ){
   return a.first < b.first;
}

std::string get_background_context(const std::string &state,
                                   const size_t background_order) {
  if (state.size() <= background_order) {
    return state; // <- This will never happen 
  } else {
    size_t background = state.size() - background_order;
    return state.substr(background);
  }
}

TEST_F(RIKmerTest, showIntRep) {
  vlmc::ReadInKmer ri_kmer{};
  std::vector<std::string> strings{""};
  std::string current_string{""};  
  createStrComb(strings, current_string, 5);

  std::vector<std::pair<int, std::string>> cmb{};

  for (auto str : strings){
    std::string kmer_string{str};
    auto old_kmer = create_kmer(kmer_string);
    vlmc::ReadInKmer kmer{old_kmer};
    cmb.push_back({kmer.integer_rep, str}); 
  }

  sort(cmb.begin(), cmb.end(), compare);

  for (int background_order = 0; background_order <= 3; background_order++){
    for (auto item : cmb) {
      auto idx = item.first;
      auto context = get_background_context(item.second, background_order);
      std::string kmer_string{context};
      auto old_kmer = create_kmer(kmer_string);
      vlmc::ReadInKmer kmer{old_kmer};
      auto context_idx = kmer.integer_rep;
      auto our_context_idx = ri_kmer.background_order_index(idx, background_order);
      EXPECT_EQ(context_idx, our_context_idx);
    }
  }

}
/*
TEST_F(RIKmerTest, RI_KmerTester) {
  std::array<double, 4> input_counts = {3,6,8,9};
  int sum_children = std::accumulate(input_counts.begin(), input_counts.end(), 1 * 4);
  container::ReadInKmer tester = container::ReadInKmer(input_counts, 1);
  EXPECT_EQ(30, sum_children);
}
*/