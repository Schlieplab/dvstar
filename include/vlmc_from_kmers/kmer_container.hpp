#pragma once

#include <execution>
#include <functional>
#include <memory>
#include <unordered_map>
#include <robin_hood.h>

#include <stxxl/sort>
#include <stxxl/sorter>

#include "kmer.hpp"

namespace vlmc {

enum Core { out, in, hash };
enum class Iteration { parallel, sequential };

constexpr int max_k = 255;

template <class Comparator = ReverseKMerComparator<max_k>> class KmerContainer {
  static_assert(
      std::is_base_of<VirtualKMerComparator<max_k>, Comparator>::value,
      "Invalid comparator template");

public:
  KmerContainer() = default;
  ~KmerContainer() = default;

  VLMCKmer null_kmer{};
  [[nodiscard]] virtual size_t size() const { return 0; };
  virtual void push(const VLMCKmer &kmer){};
  virtual void push(VLMCKmer &kmer){};
  virtual void sort(){};
  virtual void clear(){};
  virtual void for_each(const std::function<void(VLMCKmer &kmer)> &){};
  virtual VLMCKmer &get(const VLMCKmer &kmer) { return null_kmer; };
};

template <class Comparator = ReverseKMerComparator<max_k>>
class InCoreKmerContainer : public KmerContainer<Comparator> {
  std::deque<VLMCKmer> container{};

public:
  InCoreKmerContainer() = default;
  explicit InCoreKmerContainer(Iteration iteration_) : iteration(iteration_){};
  ~InCoreKmerContainer() = default;

  Iteration iteration = Iteration::parallel;

  [[nodiscard]] size_t size() const override { return container.size(); }

  void push(const VLMCKmer &kmer) override { container.push_back(kmer); }
  void push(VLMCKmer &kmer) override { container.push_back(kmer); }

  void sort() override {
    if (iteration == Iteration::parallel) {
      std::sort(std::execution::par_unseq, container.begin(), container.end(),
                Comparator());
    } else {
      std::sort(container.begin(), container.end(), Comparator());
    }
  };

  void clear() override { container.clear(); };

  void for_each(const std::function<void(VLMCKmer &kmer)> &f) override {
    std::for_each(container.begin(), container.end(), f);
  }
};

template <class Comparator = ReverseKMerComparator<max_k>>
class HashMapKmerContainer : public KmerContainer<Comparator> {
  robin_hood::unordered_map<VLMCKmer, VLMCKmer> container{};

public:
  HashMapKmerContainer() = default;
  explicit HashMapKmerContainer(Iteration iteration_) : iteration(iteration_){};
  ~HashMapKmerContainer() = default;

  Iteration iteration = Iteration::parallel;

  [[nodiscard]] size_t size() const override { return container.size(); }

  void push(const VLMCKmer &kmer) override { container[kmer] = kmer; }
  void push(VLMCKmer &kmer) override { container[kmer] = kmer; }

  void sort() override{};

  void clear() override { container.clear(); };

  void for_each(const std::function<void(VLMCKmer &kmer)> &f) override {
    std::vector<std::string> keys{};
    for (auto &[k, v] : container) {
      f(container[k]);
    }
  }

  bool contains(const VLMCKmer &kmer) {
    return container.find(kmer) != container.end();
  }

  VLMCKmer &get(const VLMCKmer &kmer) { return container[kmer]; }
};

template <class Comparator = ReverseKMerComparator<max_k>>
class OutOfCoreKmerContainer : public KmerContainer<Comparator> {
  template <int MAX_K>
  using kmer_sorter = stxxl::sorter<VLMCKmer, Comparator, 16 * 1024 * 1024>;
  kmer_sorter<max_k> sorter{Comparator(), 128 * 1024 * 1024};

public:
  OutOfCoreKmerContainer() = default;
  ~OutOfCoreKmerContainer() = default;

  [[nodiscard]] size_t size() const override { return sorter.size(); }

  void push(const VLMCKmer &kmer) override { sorter.push(kmer); }
  void push(VLMCKmer &kmer) override { sorter.push(kmer); }

  void sort() override { sorter.sort(); };

  void clear() override { sorter.clear(); };

  void for_each(const std::function<void(VLMCKmer &kmer)> &f) override {
    sorter.rewind();

    VLMCKmer kmer{};

    while (!sorter.empty()) {
      kmer = *sorter;

      f(kmer);

      ++sorter;
    }
  }
};

template <class Comparator = ReverseKMerComparator<max_k>>
class OutOfCoreSequentialKmerContainer : public KmerContainer<Comparator> {
  typedef stxxl::VECTOR_GENERATOR<VLMCKmer>::result vector_type;
  stxxl::vector<VLMCKmer> container{};

public:
  OutOfCoreSequentialKmerContainer() = default;
  ~OutOfCoreSequentialKmerContainer() = default;

  [[nodiscard]] size_t size() const override { return container.size(); }

  void push(const VLMCKmer &kmer) override { container.push_back(kmer); }
  void push(VLMCKmer &kmer) override { container.push_back(kmer); }

  void sort() override {
    stxxl::sort(container.begin(), container.end(), Comparator(),
                512 * 1024 * 1024);
  };

  void clear() override { container.clear(); };

  void for_each(const std::function<void(VLMCKmer &kmer)> &f) override {
    std::for_each(container.begin(), container.end(), f);
  }
};
} // namespace vlmc
