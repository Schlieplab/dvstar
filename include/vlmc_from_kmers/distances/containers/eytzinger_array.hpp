#pragma once

#include <bits/stdc++.h>
#include <iostream>

#include "../read_in_kmer.hpp"

namespace vlmc::array {
struct EytzingerArray {
  ReadInKmer null_kmer = ReadInKmer(-1);
  int size{};
  static const int block_size = 2; // = 64 / sizeof(ReadInKmer)
  ReadInKmer *kmer_from{};
  alignas(64) std::vector<ReadInKmer> ey_sorted_kmers;

  EytzingerArray() = default;
  EytzingerArray(std::vector<ReadInKmer> &from_container) {
    size = from_container.size();
    kmer_from = from_container.data();
    ey_sorted_kmers.reserve(size + 1);
    ey_sorted_kmers[0] = null_kmer;
    construct();
  }
  ~EytzingerArray() = default;

  int construct(int i = 0, int k = 1) {
    if (k <= size) {
      i = EytzingerArray::construct(i, 2 * k);
      ey_sorted_kmers[k] = kmer_from[i++];
      i = EytzingerArray::construct(i, 2 * k + 1);
    }
    return i;
  }

  int search(int x) {
    int k = 1;
    while (k <= size) {
      __builtin_prefetch(ey_sorted_kmers.data() + k * block_size);
      k = 2 * k + (ey_sorted_kmers[k] < x);
    }
    k >>= __builtin_ffs(~k);
    return k;
  }

  ReadInKmer &get_from_array(const int i_rep) {
    return ey_sorted_kmers[search(i_rep)];
  }
};
} // namespace vlmc::array