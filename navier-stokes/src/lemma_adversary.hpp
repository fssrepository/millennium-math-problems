#pragma once

#include "parallel_executor.hpp"

#include <cstddef>
#include <utility>

namespace lemma {

class LemmaAdversary {
public:
    explicit LemmaAdversary(int requested_threads)
        : executor_(requested_threads) {}

    template <typename Worker>
    void run_restarts(std::size_t restart_count, Worker&& worker) const {
        executor_.for_each(restart_count, std::forward<Worker>(worker));
    }

    [[nodiscard]] int threads() const { return executor_.threads(); }

private:
    ParallelExecutor executor_;
};

}  // namespace lemma
