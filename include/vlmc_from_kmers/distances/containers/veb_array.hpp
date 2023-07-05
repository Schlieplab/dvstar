#pragma once

#include <functional>
#include <iostream>
#include <limits.h>
#include <memory>

#include "../read_in_kmer.hpp"

namespace vlmc::array {
struct VanEmdeBoasArray {
  alignas(64) std::vector<ReadInKmer> a;
  static const unsigned MAX_H = 32;
  int height;
  int n;
  typedef unsigned char h_type;
  struct dumdum {
    h_type h0;
    h_type h1;
    h_type dummy[2];
    int m0;
    int m1;
  };
  dumdum s[MAX_H + 1];

  VanEmdeBoasArray() = default;
  ~VanEmdeBoasArray() = default;

  void sequencer(int h, dumdum *s, unsigned d) {
    if (h == 0)
      return;
    int h0 = h / 2;
    int h1 = h - h0 - 1;
    sequencer(h0, s, d);
    s[d + h0].h0 = h0;
    s[d + h0].m0 = (2 << h0) - 1;
    s[d + h0].h1 = h1;
    s[d + h0].m1 = (2 << h1) - 1;
    sequencer(h1, s, d + h0 + 1);
  }

  ReadInKmer *construct(ReadInKmer *a0, int *rtl, int path, unsigned d) {
    if (d > height || rtl[d] >= n)
      return a0;

    // visit left child
    path <<= 1;
    rtl[d + 1] = rtl[d - s[d].h0] + s[d].m0 + (path & s[d].m0) * (s[d].m1);
    a0 = construct(a0, rtl, path, d + 1);

    a[rtl[d]] = *a0++;

    // visit right child
    path += 1;
    rtl[d + 1] = rtl[d - s[d].h0] + s[d].m0 + (path & s[d].m0) * (s[d].m1);
    a0 = construct(a0, rtl, path, d + 1);

    return a0;
  }

  VanEmdeBoasArray(std::vector<ReadInKmer> &from_container) {
    n = from_container.size();
    // find smallest h such that sum_i=0^h 2^h >= n
    int m = 1;
    for (height = 0; m < n; height++, m += 1 << height)
      ;

    dumdum q = {(h_type)height, 0, {0, 0}, (2 << height) - 1, 1};
    std::fill_n(s, MAX_H + 1, q);
    sequencer(height, s, 0);

    a.resize(n);
    int rtl[MAX_H + 1];
    rtl[0] = 0;
    construct(from_container.data(), rtl, 0, 0);
  }

  int search(int x) {
    int rtl[MAX_H + 1];
    int j = n;
    int i = 0;
    int p = 0;
    for (int d = 0; i < n; d++) {
      rtl[d] = i;
      if (x < a[i].integer_rep) {
        p <<= 1;
        j = i;
      } else if (x > a[i].integer_rep) {
        p = (p << 1) + 1;
      } else {
        return i;
      }
      i = rtl[d - s[d].h0] + s[d].m0 + (p & s[d].m0) * (s[d].m1);
    }
    return j;
  }

  ReadInKmer &get_from_array(const int i_rep) { return a[search(i_rep)]; }
};
} // namespace vlmc::array