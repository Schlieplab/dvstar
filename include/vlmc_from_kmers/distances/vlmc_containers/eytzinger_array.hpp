#pragma once

#include <iostream>

#include "read_in_kmer.hpp"
#include <bits/stdc++.h>

namespace array {
  struct Ey_array {
    kmers::RI_Kmer null_kmer = kmers::RI_Kmer(-1);
    int size;
    static const int block_size = 2; // = 64 / sizeof(RI_Kmer)
    kmers::RI_Kmer* kmer_from;
    alignas(64) std::vector<kmers::RI_Kmer> ey_sorted_kmers;

    Ey_array() = default;
    Ey_array(std::vector<kmers::RI_Kmer>& from_container) {
      size = from_container.size();
      kmer_from = from_container.data();
      ey_sorted_kmers.reserve(size + 1);
      ey_sorted_kmers[0] = null_kmer;
      construct();
    }
    ~Ey_array() = default;

    int construct(int i = 0, int k = 1) {
      if (k <= size) {
        i = Ey_array::construct(i, 2 * k);
        ey_sorted_kmers[k] = kmer_from[i++];
        i = Ey_array::construct(i, 2 * k + 1);
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

    kmers::RI_Kmer& get_from_array(const int i_rep) {
      return ey_sorted_kmers[search(i_rep)];
    }
  };
}