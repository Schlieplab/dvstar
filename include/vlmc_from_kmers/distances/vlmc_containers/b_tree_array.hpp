#pragma once

#include <iostream>

#include "read_in_kmer.hpp"
#include <bits/stdc++.h>

namespace array {
	struct B_Tree {
		int size;
		const int block_size = 2; // 64 / sizeof(RI_Kmer)
		static const int B = 2;
		alignas(64) std::vector<kmers::RI_Kmer> a;

		B_Tree() = default;
		~B_Tree() = default;

		B_Tree(std::vector<kmers::RI_Kmer>& from_container) {
			size = from_container.size();
			a.reserve(size + 1);
			construct(from_container.begin(), 0);
		}
		static int child(unsigned c, int i) {
			return (B + 1) * i + (c + 1) * B;
		}

		template <unsigned int C>
		static const kmers::RI_Kmer* branchfree_inner_search(const kmers::RI_Kmer* base, const kmers::RI_Kmer x) {
			if (C <= 1)
				return base;
			const unsigned int half = C / 2;
			const kmers::RI_Kmer* current = &base[half];
			return branchfree_inner_search<C - half>((*current < x) ? current : base, x);
		}

		template <unsigned C>
		static int branchy_inner_search(const kmers::RI_Kmer* a, int i, kmers::RI_Kmer x) {
			if (C == 0)
				return i;
			if (x <= a[i + C / 2].integer_rep)
				return branchy_inner_search<C / 2>(a, i, x);
			return branchy_inner_search<C - C / 2 - 1>(a, i + C / 2 + 1, x);
		}

		std::vector<kmers::RI_Kmer>::iterator construct(std::vector<kmers::RI_Kmer>::iterator a0, int i) {
			if (i >= size)
				return a0;

			for (unsigned c = 0; c <= B; c++) {
				// visit c'th child
				a0 = construct(a0, child(c, i));
				if (c < B && i + c < size) {
					a[i + c] = *a0++;
				}
			}
			return a0;
		}

		int search(int x) {
			int j = size;
			int i = 0;
			while (i < size) {
				int lo = i;
				int hi = std::min(i + B, size);
				while (lo < hi) {
					int m = (lo + hi) / 2;
					if (x < a[m].integer_rep) {
						hi = m;
						j = hi;
					}
					else if (x > a[m].integer_rep) {
						lo = m + 1;
					}
					else {
						return m;
					}
				}
				i = child((unsigned)(hi - i), i);
			}
			return j;
		}

		// unrolled branchy inner serach
		int unrolled_branchy_search(int x) const {
			int j = size;
			int i = 0;
			while (i + B <= size) {
				int t = branchy_inner_search<B>(a.data(), i, x);
				j = t < i + B ? t : j;
				i = child((unsigned)(t - i), i);
			}
			if (__builtin_expect(i <= size, 0)) {
				// Now we're in the last block
				int lo = i;
				int hi = size;
				while (lo < hi) {
					int m = (lo + hi) / 2;
					if (x < a[m].integer_rep) {
						hi = m;
						j = m;
					}
					else if (x > a[m].integer_rep) {
						lo = m + 1;
					}
					else {
						return m;
					}
				}
			}
			return j;
		}

		// branch-free search (with or without prefetching)
		int unrolled_branchfree_search(int x) const {
			int j = size;
			int i = 0;
			while (i + B <= size) {
				__builtin_prefetch(a.data() + child(i, B / 2), 0, 0);
				const kmers::RI_Kmer* base = &a[i];
				const kmers::RI_Kmer* pred = branchfree_inner_search<B>(base, x);
				unsigned int nth = (*pred < x) + pred - base;
				{
					/* nth == B iff x > all values in block. */
					const kmers::RI_Kmer current = base[nth % B];
					int next = i + nth;
					j = (current >= x) ? next : j;
				}
				i = child(nth, i);
			}
			if (__builtin_expect(i < size, 0)) {
				// last (partial) block
				const kmers::RI_Kmer* base = &a[i];
				int m = size - i;
				while (m > 1) {
					int half = m / 2;
					const kmers::RI_Kmer* current = &base[half];

					base = (*current < x) ? current : base;
					m -= half;
				}

				int ret = (*base < x) + base - a.data();
				return (ret == size) ? j : ret;
			}
			return j;
		}

		kmers::RI_Kmer& get_from_array(const int i_rep) {
			return a[unrolled_branchfree_search(i_rep)];
		}
	};
}