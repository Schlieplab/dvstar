#pragma once

#include <filesystem>
#include <functional>
#include <unordered_map>

#include "vlmc_container.hpp"

namespace vlmc::container{
template <typename VC> class ClusterContainer {

private:
  std::vector<VC> container{};

public:
  ClusterContainer() = default;
  ~ClusterContainer() = default;

  ClusterContainer(const size_t i) : container(i) {}

  size_t size() const { return container.size(); }

  void push(const VC vlmc) { container.push_back(vlmc); }

  VC &get(const int i) { return container[i]; }

  VC &operator[](size_t index) { return container[index]; }

  const VC &operator[](size_t index) const { return container[index]; }
};

struct KmerPair {
  ReadInKmer kmer;
  size_t id;

  KmerPair(ReadInKmer kmer, size_t id) {
    this->kmer = kmer;
    this->id = id;
  }

  KmerPair() = default;
  ~KmerPair() = default;
};

class KmerCluster {

private:
  std::unordered_map<int, std::vector<KmerPair>> container{};

  size_t vlmc_count = 0;

public:
  KmerCluster() = default;
  ~KmerCluster() = default;

  int size() const { return vlmc_count; }

  void set_size(size_t count) { this->vlmc_count = count; }

  void push(const KmerPair kmer_pair) {
    container[kmer_pair.kmer.integer_rep].push_back(kmer_pair);
  }

  void push_all(KmerCluster cluster) {
    auto begin_it = cluster.get_begin();
    auto end_it = cluster.get_end();
    while (begin_it != end_it) {
      container[begin_it->first].insert(container[begin_it->first].end(),
                                        begin_it->second.begin(),
                                        begin_it->second.end());
      begin_it++;
    }
  }

  std::vector<KmerPair> &get(int bucket_num) { return container[bucket_num]; }

  int experimental_bucket_count() {
    auto count = 0;
    for (auto it = container.begin(); it != container.end(); it++) {
      count++;
    }
    return count;
  }

  std::unordered_map<int, std::vector<KmerPair>>::iterator get_begin() {
    return container.begin();
  }

  std::unordered_map<int, std::vector<KmerPair>>::iterator get_end() {
    return container.end();
  }

  std::unordered_map<int, std::vector<KmerPair>>::iterator
  find(int bucket_num) {
    return container.find(bucket_num);
  }
};
} // namespace vlmc::cluster_container