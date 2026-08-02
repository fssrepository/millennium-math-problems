#pragma once

#include "parallel_executor.hpp"

#include <cstddef>
#include <utility>

namespace lemma {

class ProjectiveFamily {
public:
    explicit ProjectiveFamily(int requested_threads)
        : executor_(requested_threads) {}

    template <typename Worker>
    void run_cutoffs(std::size_t cutoff_count, bool kernel_parallel,
                     Worker&& worker) const {
        if (kernel_parallel) {
            for (std::size_t index = 0; index < cutoff_count; ++index) {
                worker(index);
            }
        } else {
            executor_.for_each(cutoff_count, std::forward<Worker>(worker));
        }
    }

    [[nodiscard]] int threads() const { return executor_.threads(); }

private:
    ParallelExecutor executor_;
};

}  // namespace lemma
