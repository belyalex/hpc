//
// Created by alex on 30.09.2020.
//

#ifndef MULMATR_TIMER_H
#define MULMATR_TIMER_H


#include <chrono>
#include <iostream>

class timer {
    using clock_t = std::chrono::high_resolution_clock;
    using nanoseconds = std::chrono::nanoseconds;
public:
    explicit timer(int d = 1)
            : d(d), start_(clock_t::now()) {
    }

    ~timer() {
        const auto finish = clock_t::now();
        const auto duration =
                std::chrono::duration_cast<nanoseconds>
                        (finish - start_).count();
        std::cout << duration / d / 1.e9;
    }

private:
    const int d;
    const clock_t::time_point start_;
};

#endif //MULMATR_TIMER_H
