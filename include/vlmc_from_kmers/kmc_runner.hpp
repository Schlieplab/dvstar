#pragma once

#include <filesystem>
#include <string>

#include "cli_helper.hpp"

namespace vlmc {

void remove_kmc_files(const std::filesystem::path &kmc_path) {
  std::filesystem::path path_copy{kmc_path};
  std::filesystem::remove(path_copy.replace_extension(".kmc_suf"));
  std::filesystem::remove(path_copy.replace_extension(".kmc_pre"));
}

std::string find_kmc() {
  std::filesystem::path local_path{"./kmc"};
  if (std::filesystem::exists(local_path)) {
    return local_path.string();
  } else {
    return "kmc";
  }
}

std::string find_kmc_tools() {
  std::filesystem::path local_path{"./kmc_tools"};
  if (std::filesystem::exists(local_path)) {
    return local_path.string();
  } else {
    return "kmc_tools";
  }
}

std::filesystem::path run_kmc(const std::filesystem::path &fasta_path,
                              const int kmer_size,
                              const std::filesystem::path &tmp_path,
                              const Core &in_or_out_of_core) {
  std::string random_name = get_random_name("kmc_");

  std::filesystem::path kmc_db_path = tmp_path / (random_name + "_res");
  std::filesystem::path kmc_db_sorted_path =
      tmp_path / (random_name + "_res_sorted");
  std::filesystem::path kmc_db_sorted_path_with_suffix =
      tmp_path / (random_name + "_res_sorted.kmc_suf");

  std::filesystem::path kmc_tmp = tmp_path / (random_name + "_tmp");
  std::filesystem::create_directories(kmc_tmp);

  auto kmc_run = std::chrono::steady_clock::now();

  std::ostringstream kmc_run_stream;
  kmc_run_stream << find_kmc();
  kmc_run_stream << " -b -ci1 -cs4294967295 ";

  if (in_or_out_of_core == Core::in) {
    kmc_run_stream << "-r ";
  }
  kmc_run_stream << "-k" << kmer_size << " ";
  if (fasta_path.extension() == ".fastq") {
    kmc_run_stream << "-fm " << fasta_path << " ";
  } else {
    kmc_run_stream << "-fm " << fasta_path << " ";
  }
  kmc_run_stream << kmc_db_path << " " << kmc_tmp;
  std::string kmc_run_command = kmc_run_stream.str();

  auto kmc_run_code = system(kmc_run_command.c_str());
  if (kmc_run_code != 0) {
      throw std::runtime_error("kmc command " + kmc_run_command + " failed.");
  }
  auto kmc_run_done = std::chrono::steady_clock::now();

  auto kmc_sort = std::chrono::steady_clock::now();

  std::ostringstream kmc_sort_stream;
  kmc_sort_stream << find_kmc_tools();
  kmc_sort_stream << " transform ";
  kmc_sort_stream << kmc_db_path;
  kmc_sort_stream << " sort ";
  kmc_sort_stream << kmc_db_sorted_path;
  std::string kmc_sort_command = kmc_sort_stream.str();

  auto kmc_sort_code = system(kmc_sort_command.c_str()); // TODO is it possible / permitted under GPL
                                    // v3 to compile kmc into our binary?
  if (kmc_sort_code != 0) {
      throw std::runtime_error("kmc_tools command " + kmc_sort_command + " failed.");
  }
  auto kmc_sort_done = std::chrono::steady_clock::now();

  std::chrono::duration<double> kmc_run_seconds = kmc_run_done - kmc_run;
  std::cout << "KMC count time: " << kmc_run_seconds.count() << "s\n";

  std::chrono::duration<double> kmc_sort_seconds = kmc_sort_done - kmc_sort;
  std::cout << "KMC sort time: " << kmc_sort_seconds.count() << "s\n";

  std::filesystem::remove_all(kmc_tmp);

  if (std::filesystem::exists(kmc_db_sorted_path_with_suffix)) {
    remove_kmc_files(kmc_db_path);
    return kmc_db_sorted_path;
  } else {
    return kmc_db_path;
  }
}
} // namespace vlmc
