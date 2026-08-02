#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace lemma {

class ParallelExecutor {
public:
    explicit ParallelExecutor(int requested_threads = 0) {
        const unsigned hardware = std::thread::hardware_concurrency();
        const int automatic = static_cast<int>(hardware == 0 ? 1U : hardware);
        threads_ = requested_threads == 0 ? std::min(12, automatic) : requested_threads;
        if (threads_ < 1 || threads_ > 256) {
            throw std::invalid_argument("thread count must be between 1 and 256");
        }
    }

    [[nodiscard]] int threads() const { return threads_; }

    template <typename Function>
    void for_each(std::size_t task_count, Function&& function) const {
        if (task_count == 0) {
            return;
        }
        const std::size_t worker_count = std::min<std::size_t>(
            task_count, static_cast<std::size_t>(threads_));
        if (worker_count == 1) {
            for (std::size_t index = 0; index < task_count; ++index) {
                function(index);
            }
            return;
        }

        std::atomic<std::size_t> next{0};
        std::atomic<bool> failed{false};
        std::exception_ptr exception;
        std::mutex exception_mutex;
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        for (std::size_t worker = 0; worker < worker_count; ++worker) {
            workers.emplace_back([&] {
                while (!failed.load(std::memory_order_relaxed)) {
                    const std::size_t index =
                        next.fetch_add(1, std::memory_order_relaxed);
                    if (index >= task_count) {
                        break;
                    }
                    try {
                        function(index);
                    } catch (...) {
                        failed.store(true, std::memory_order_relaxed);
                        std::lock_guard<std::mutex> lock(exception_mutex);
                        if (!exception) {
                            exception = std::current_exception();
                        }
                    }
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }
        if (exception) {
            std::rethrow_exception(exception);
        }
    }

private:
    int threads_ = 1;
};

}  // namespace lemma
