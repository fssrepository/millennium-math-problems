#pragma once

#include "triad_partition.hpp"

#include <string_view>

namespace lemma {

class LocalSldTriadSelection {
public:
    [[nodiscard]] static bool supports(std::string_view name);
    [[nodiscard]] static TriadSelection parse(std::string_view name);
};

}  // namespace lemma
