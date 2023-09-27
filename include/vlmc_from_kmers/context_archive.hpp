#pragma once

#include <cereal/archives/binary.hpp>
#include <fstream>

#include "kmer.hpp"
namespace vlmc::details {
bool has_next(std::ifstream &file_stream) { return file_stream.peek() != EOF; }

VLMCKmer &next(cereal::BinaryInputArchive &iarchive, VLMCKmer &kmer) {
  iarchive(kmer);
  return kmer;
}

} // namespace vlmc::details

namespace vlmc {
void iterate_archive(const std::filesystem::path &path,
                     const std::function<void(const VLMCKmer &)> &f) {
  std::ifstream fs(path, std::ios::binary);
  cereal::BinaryInputArchive archive(fs);

  try {
    VLMCKmer kmer{};
    while (details::has_next(fs)) {
      details::next(archive, kmer);
      f(kmer);
    }
    fs.close();
  } catch (cereal::Exception &e) {
    std::cerr << path.string() << " could be be loaded." << std::endl;
    throw e;
  }


}
} // namespace vlmc
