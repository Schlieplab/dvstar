#pragma once

#include <filesystem>
#include <format>
#include <vector>

namespace vlmc {
std::vector<std::filesystem::path>
get_recursive_paths(const std::filesystem::path &path,
                    const std::set<std::string> &allowed_extensions) {
  std::vector<std::filesystem::path> paths{};
  if (is_directory(path)) {
    for (const auto &dir_entry :
         std::filesystem::recursive_directory_iterator(path)) {
      if (allowed_extensions.find(dir_entry.path().extension()) !=
          allowed_extensions.end()) {
        paths.push_back(dir_entry.path());
      }
    }
  } else if (allowed_extensions.find(path.extension()) !=
             allowed_extensions.end()) {
    paths.push_back(path);
  } else {
    std::string extension_string{};
    for (const auto &extension : allowed_extensions) {
      extension_string += extension + ", ";
    }
    std::string msg =
        std::format("Path {} is not a path to a directory or a {} file.\n",
                    path.string(), extension_string);
    throw std::invalid_argument(msg);
  }
  return paths;
}
} // namespace vlmc