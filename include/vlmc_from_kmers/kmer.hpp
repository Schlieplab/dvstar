#pragma once

#include <kmc_file.h>

#include <array>
#include <functional>
#include <unordered_map>

#include <cereal/archives/binary.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/array.hpp>

namespace vlmc {
// For stxxl, this class needs to be POD
// Otherwise, borrows implementations from KMCs CKmerApi,
// but disregards the byte alignment.

static const std::unordered_map<char, uchar> code_codes{
    {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
static constexpr char char_codes[4] = {'A', 'C', 'G', 'T'};


struct VLMCKmer {
  VLMCKmer() = default;
  VLMCKmer(uint32 length_, uint64 count_,
           std::array<uint64, 4> next_symbol_counts_)
      : length(length_), count(count_), next_symbol_counts(next_symbol_counts_),
        divergence(-1.0), kmer_data(), is_terminal(false), has_children(false),
        to_be_removed(true) {
    this->n_rows = (length_ / 32.0) + 1; // std::ceil, but no float conversion
    if (length_ == 0) {
      this->n_rows = 0;
    }
  }

  ~VLMCKmer() = default;

  std::array<uint64, 4> kmer_data;
  uint64 count;
  std::array<uint64, 4> next_symbol_counts;
  double divergence;
  uint32 length;
  uint32 n_rows;
  bool is_terminal;
  bool has_children;
  bool to_be_removed;

  static VLMCKmer create_prefix_kmer(const VLMCKmer &kmer, uint32 length,
                                     uint64 count,
                                     const std::array<uint64, 4> &child_counts,
                                     VLMCKmer &out_kmer) {
    out_kmer.length = length;
    out_kmer.count = count;
    out_kmer.next_symbol_counts = child_counts;

    for (int i = 0; i < 4; i++) {
      out_kmer.kmer_data[i] = 0;
    }

    uint32 row = length >> 5;

    for (int i = 0; i < row; i++) {
      out_kmer.kmer_data[i] = kmer.kmer_data[i];
    }

    uint32 length_in_row = length & 31;

    uint32 positions_to_remove_from_row = (31 - length_in_row) * 2;
    unsigned long long mask = 0xFFFFFFFFFFFFFFFF
                              << positions_to_remove_from_row;

    out_kmer.kmer_data[row] = (kmer.kmer_data[row] & mask);

    return out_kmer;
  }

  static VLMCKmer
  create_prefix_kmer(const VLMCKmer &kmer, uint32 length, uint64 count,
                     const std::array<uint64, 4> &child_counts) {
    VLMCKmer prefix_kmer{length, count, child_counts};

    return create_prefix_kmer(kmer, length, count, child_counts, prefix_kmer);
  }

  static VLMCKmer create_suffix_kmer(const VLMCKmer &kmer, VLMCKmer &out_kmer) {
    out_kmer.length = kmer.length - 1;
    if (out_kmer.length == 0) {
      out_kmer.n_rows = 0;
      for (int i = 0; i < 4; i++) {
        out_kmer.kmer_data[i] = 0;
      }
      return out_kmer;
    }

    uint32 row = out_kmer.length >> 5;
    out_kmer.n_rows = row + 1;

    out_kmer.kmer_data[0] = kmer.kmer_data[0] << 2; // remove first char
    for (int i = 1; i < row; i++) {
      // Move char from this row to the previous row.
      unsigned long long remove_mask = 0xC000000000000000;
      unsigned long long move_data = kmer.kmer_data[i] & remove_mask;
      out_kmer.kmer_data[i] = kmer.kmer_data[i] << 2;
      out_kmer.kmer_data[i - 1] = kmer.kmer_data[i - 1] | (move_data >> 62);
    }

    return out_kmer;
  }

  static VLMCKmer replace_first_char(const VLMCKmer &kmer, uchar replace_with,
                                     VLMCKmer &out_kmer) {
    out_kmer.length = kmer.length;

    for (int i = 0; i < 4; i++) {
      out_kmer.kmer_data[i] = kmer.kmer_data[i];
    }

    unsigned long long remove_mask = 0x3FFFFFFFFFFFFFFF;
    unsigned long long new_data = replace_with;
    new_data <<= 62;

    out_kmer.kmer_data[0] = (out_kmer.kmer_data[0] & remove_mask) | new_data;

    return out_kmer;
  }

  // This method lets cereal know which data members to serialize
  template <class Archive> void serialize(Archive &archive) {
    archive(kmer_data, length, n_rows, count, next_symbol_counts, divergence,
            is_terminal);
  }

  /**
   * Finds the first position between the two kmers that is different.
   * E.g. for two kmers "ACGT" and "ACTT", it would return 2
   *
   * @param current_kmer
   * @param prev_kmer
   * @return position the k-mers differ at.
   */
  static int32 get_first_differing_position(
      const VLMCKmer &current_kmer, const VLMCKmer &prev_kmer,
      const uint32 current_kmer_prefix = 0, const uint32 prev_kmer_prefix = 0) {
    auto n_rows = current_kmer.n_rows;
    int32 offset = 0;

    for (int32 row_counter = 0; row_counter < n_rows; row_counter++) {
      unsigned long long current_data = current_kmer.kmer_data[row_counter];
      unsigned long long prev_data = prev_kmer.kmer_data[row_counter];

      if (row_counter == 0) {
        current_data = current_data << (current_kmer_prefix * 2);
        prev_data = prev_data << (prev_kmer_prefix * 2);
      }

      unsigned long long diff = (current_data) ^ (prev_data);

      // GCC instruction, counts leading zeros.
      int32 n_leading_zeros = __builtin_clzll(diff);

      if (diff != 0 && n_leading_zeros != 63) {
        int32 diff_pos = n_leading_zeros / 2;
        int32 final_diff_pos = diff_pos + offset;

        if (final_diff_pos >= current_kmer.length - current_kmer_prefix) {
          return -1;
        } else {
          return final_diff_pos;
        }
      }
      offset += 32;
    }
    return -1;
  }

  [[nodiscard]] uchar extract2bits(uint32 pos) const {
    // pos >> 5 == pos // 31
    // pos & 31 == remainder(pos, 31)

    uchar row = pos >> 5;
    uchar pos_in_row = pos & 31;
    uchar n_shift_pos_to_end = (62 - pos_in_row * 2);
    return (kmer_data[row] >> n_shift_pos_to_end) & 3;
  }

  [[nodiscard]] uchar char_pos(uint32 pos) const { return extract2bits(pos); }

  [[nodiscard]] static int n_kmers_with_length(uint32_t len) {
    static std::vector<int> vals = {1, 5, 21, 85, 341, 1365};
    //    return vals[len];

    if (len > vals.size()) {
      return (1 - std::pow(4, len + 1)) / (1 - 4);
    } else {
      return vals[len];
    }
  }

  /*
   * Get the value of the first prefix_length characters.
   *
   * Useful for splitting the kmers into subsets based on their first
   * `prefix_length` characters.
   *
   * WARNING: Assumes the prefix is smaller than 32, as this is our only
   * use case.
   *
   * @param prefix_length length of prefix.
   * @return unique value per prefix.
   */
  [[nodiscard]] uint64 get_prefix_index(int prefix_length) const {
    uint64 key = 0;

    int n_shift_pos_to_end = (32 - prefix_length) * 2;
    return kmer_data[0] >> n_shift_pos_to_end;

    for (int i = 0; i < prefix_length; i++) {
      int bits = extract2bits(i);

      // Length is reversed - suffixes should count as having only the most
      // significant digits.
      int length_from_end = prefix_length - i - 1;

      // For every position we add the number of k-mers with a specific
      // length that have to come before this one, multiplied by the bit values.
      // This can be viewed as dividing the space of all prefix_length k-mers
      // into divisions for every suffix, where every such division at a
      // specific length is n_kmers_with_length(i) large.  Thus, we need
      // to find which such division this kmer belongs to.
      // Plus one so e.g. AT and T don't get same values.

      key += bits * VLMCKmer::n_kmers_with_length(length_from_end) + 1;
    }

    return key;
  }

  [[nodiscard]] bool reverse_less_than(const VLMCKmer &kmer) const {
    int this_pos = this->length - 1;
    int kmer_pos = kmer.length - 1;

    while (this_pos >= 0 && kmer_pos >= 0) {
      auto this_2_bits = this->extract2bits(this_pos--);
      auto kmer_2_bits = kmer.extract2bits(kmer_pos--);

      if (this_2_bits != kmer_2_bits) {
        return this_2_bits > kmer_2_bits;
      }
    }

    if (kmer.length != this->length) {
      return kmer.length < this->length;
    }
    return kmer.count > this->count;
  }

  inline bool operator<(const VLMCKmer &kmer) const {
    auto min_length = std::min(kmer.length, length);

    for (int i = 0; i < min_length; i++) {
      auto this_2_bits = this->extract2bits(i);
      auto kmer_2_bits = kmer.extract2bits(i);
      if (this_2_bits != kmer_2_bits) {
        return this_2_bits < kmer_2_bits;
      }
    }

    if (kmer.length != this->length) {
      return kmer.length > this->length;
    }
    return kmer.count < this->count;
  };

  inline bool operator==(const VLMCKmer &kmer) const {
    if (kmer.length != length) {
      return false;
    }

    for (int i = 0; i < length; i++) {
      auto this_2_bits = this->extract2bits(i);
      auto kmer_2_bits = kmer.extract2bits(i);
      if (this_2_bits != kmer_2_bits) {
        return false;
      }
    }

    return true;
  };

  [[nodiscard]] std::string to_string() const {
    if (this->length <= 0) {
      return "";
    }
    std::string out_string(this->length, ' ');

    uchar *byte_ptr;
    uchar c;
    uint32 cur_string_size = 0;

    for (uint32 row_counter = 0; row_counter < this->n_rows; row_counter++) {
      byte_ptr = reinterpret_cast<uchar *>(
          &(const_cast<VLMCKmer *>(this)->kmer_data[row_counter]));

      byte_ptr += 7; // Read from the left

      for (uint32 i = 0; (i < length) && (i < 32); i += 4) {

        c = 0xc0 & *byte_ptr; // 11000000
        c = c >> 6;
        out_string[cur_string_size++] = char_codes[c];

        if (cur_string_size == length)
          break;

        c = 0x30 & *byte_ptr; // 00110000
        c = c >> 4;
        out_string[cur_string_size++] = char_codes[c];
        if (cur_string_size == length)
          break;

        c = 0x0c & *byte_ptr; // 00001100
        c = c >> 2;
        out_string[cur_string_size++] = char_codes[c];
        if (cur_string_size == length)
          break;

        c = 0x03 & *byte_ptr; // 00000011
        out_string[cur_string_size++] = char_codes[c];
        if (cur_string_size == length)
          break;

        byte_ptr--;
      }
    }

    return out_string;
  }

  void output(std::ostream &stream) const {
    stream << this->to_string() << " " << this->count << " ";
    for (auto &c : this->next_symbol_counts) {
      stream << c << " ";
    }
    stream << divergence << " ";
    stream << this->is_terminal << " ";
    stream << this->to_be_removed << std::endl;
  }
};

class VLMCTranslator : public CKmerAPI {
public:
  explicit VLMCTranslator(int length) : CKmerAPI(length) {}
  ~VLMCTranslator() = default;

  VLMCKmer construct_vlmc_kmer() {
    // To make sure the VLMCKmer class stays a POD, we're only allowing k-mers
    // up to 128 or so.  Should be fine.
    VLMCKmer new_kmer{this->kmer_length, 0, {}};

    return construct_vlmc_kmer(new_kmer);
  }

  VLMCKmer construct_vlmc_kmer(VLMCKmer &new_kmer) {
    // To make sure the VLMCKmer class stays a POD, we're only allowing k-mers
    // up to 128 or so.  Should be fine.

    for (int i = 0; i < 4; i++) {
      new_kmer.kmer_data[i] = 0;
    }
    new_kmer.length = this->kmer_length;

    if (this->kmer_length == 0) {
      return new_kmer;
    }

    if (this->kmer_length < 31 - this->byte_alignment) {
      new_kmer.kmer_data[0] = this->kmer_data[0] << (this->byte_alignment * 2);

      return new_kmer;
    } else if (this->kmer_length > 31 - this->byte_alignment) {
      for (int pos = 0; pos < this->kmer_length; pos++) {
        uint32 row = pos >> 5;
        uint64 val = this->extract2bits(pos);
        new_kmer.kmer_data[row] += val << (31 - (pos & 31)) * 2;
      }

      return new_kmer;

    } else {
      new_kmer.kmer_data[0] = 0;

      return new_kmer;
    }
  }
};

template <int MAX_K> struct VirtualKMerComparator {
  virtual bool operator()(const VLMCKmer &a, const VLMCKmer &b) const = 0;
  [[nodiscard]] virtual VLMCKmer min_value() const {
    return VLMCKmer(0, 0, {});
  }
  [[nodiscard]] virtual VLMCKmer max_value() const {
    return VLMCKmer(MAX_K, 0, {});
  }
};

template <int MAX_K>
struct KMerComparator : public VirtualKMerComparator<MAX_K> {
  bool operator()(const VLMCKmer &a, const VLMCKmer &b) const override {
    return a < b;
  }
  [[nodiscard]] VLMCKmer min_value() const override {
    std::array<uint64, 4> vec{};
    return {0, 0, vec};
  }
  [[nodiscard]] VLMCKmer max_value() const override {
    VLMCKmer max_kmer(MAX_K, (size_t)-1, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
};

template <int MAX_K>
struct ReverseKMerComparator : public VirtualKMerComparator<MAX_K> {
  bool operator()(const VLMCKmer &a, const VLMCKmer &b) const override {
    return a.reverse_less_than(b);
  }
  [[nodiscard]] VLMCKmer min_value() const override {
    VLMCKmer max_kmer(MAX_K, (size_t)-1, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
  [[nodiscard]] VLMCKmer max_value() const override {
    std::array<uint64, 4> vec{};
    return {0, 0, vec};
  }
};

template <int MAX_K> struct KMerReverseKeyExtractor {
  typedef uint64 key_type;
  uint64 upper_bound = (1 - std::pow(4, MAX_K + 1)) / (1 - 4);

  [[nodiscard]] uint64 n_kmers_with_length(uint32_t len) const {
    return (1 - std::pow(4, len + 1)) / (1 - 4);
  }

  key_type operator()(const VLMCKmer &a) const {
    uint64 key = 0;

    for (int i = 0; i < a.length; i++) {
      int bits = a.extract2bits(a.length - i - 1);

      // Length is reversed - suffixes should count as having only the most
      // significant digits.
      int length_from_end = MAX_K - i - 1;

      // For every position we add the number of k-mers with a specific
      // length that have to come before this one, multiplied by the bit values.
      // This can be viewed as dividing the space of all MAX_K k-mers
      // into divisions for every suffix, where every such division at a
      // specific length is n_kmers_with_length(i) large.  Thus, we need
      // to find which such division this kmer belongs to.
      // Plus one so e.g. AT and T don't get same values.

      key += bits * this->n_kmers_with_length(length_from_end) + 1;
    }

    return upper_bound - key;
  }

  [[nodiscard]] VLMCKmer min_value() const {
    VLMCKmer max_kmer(MAX_K, (size_t)-1, {});
    for (int row = 0; row < max_kmer.n_rows; row++) {
      max_kmer.kmer_data[row] = 0xFFFFFFFFFFFFFFFF;
    }
    return max_kmer;
  }
  [[nodiscard]] VLMCKmer max_value() const {
    std::array<uint64, 4> vec{};
    return VLMCKmer(0, 0, vec);
  }
};
} // namespace vlmc

namespace std {
template <typename T> void hash_combine(size_t &seed, T const &v) {
  seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <> struct hash<vlmc::VLMCKmer> {
  std::size_t operator()(vlmc::VLMCKmer const &kmer) const noexcept {
    std::size_t seed = 0.0;
    for (int i = 0; i < kmer.n_rows; i++) {
      if (i == kmer.n_rows - 1) {
        uint32 length_in_row = (kmer.length - 1) & 31;

        uint32 positions_to_remove_from_row = (31 - length_in_row) * 2;
        unsigned long long mask = 0xFFFFFFFFFFFFFFFF
                                  << positions_to_remove_from_row;

        hash_combine(seed, kmer.kmer_data[i] & mask);
        hash_combine(seed, kmer.length);
      } else {
        hash_combine(seed, kmer.kmer_data[i]);
      }
    }

    return seed;
  }
};
} // namespace std