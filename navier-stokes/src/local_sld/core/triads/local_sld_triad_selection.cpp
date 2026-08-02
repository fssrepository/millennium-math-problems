#include "local_sld_triad_selection.hpp"

#include <stdexcept>

namespace lemma {

bool LocalSldTriadSelection::supports(std::string_view name) {
    return name == "local" ||
        name == "doubling-family" ||
        name == "doubling-remainder" ||
        name == "remainder-without-123" ||
        name == "double-triple-family" ||
        name == "double-triple-remainder" ||
        name == "double-triple-remainder-without-123";
}

TriadSelection LocalSldTriadSelection::parse(std::string_view name) {
    if (name == "local") {
        return TriadPartition::local;
    }
    if (name == "doubling-family") {
        return TriadSelection::local_equal_low_doubling();
    }
    if (name == "doubling-remainder") {
        return TriadSelection::local_without_equal_low_doubling();
    }
    if (name == "remainder-without-123") {
        return TriadSelection::
            local_without_equal_low_doubling_and_signature(1, 2, 3);
    }
    if (name == "double-triple-family") {
        return TriadSelection::local_equal_low_double_triple();
    }
    if (name == "double-triple-remainder") {
        return TriadSelection::local_without_equal_low_double_triple();
    }
    if (name == "double-triple-remainder-without-123") {
        return TriadSelection::
            local_without_equal_low_double_triple_and_signature(1, 2, 3);
    }
    throw std::invalid_argument("unsupported local SLD triad selection");
}

}  // namespace lemma
