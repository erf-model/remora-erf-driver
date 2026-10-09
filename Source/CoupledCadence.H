#ifndef ERF_REMORA_COUPLED_CADENCE_H_
#define ERF_REMORA_COUPLED_CADENCE_H_
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace erf_remora {
// Integer clocks avoid accumulated drift. Independent of AMReX so invalid
// schedules can be tested without initializing MPI, devices, or either solver.
struct CoupledCadence {
    double atmosphere_dt, ocean_dt, exchange_interval;
    int atmosphere_steps, ocean_steps;
    static bool close(double a, double b) {
        return std::abs(a-b) <= std::max(1.e-12, 1.e-8*std::max(std::abs(a), std::abs(b)));
    }
    static int integral_steps(double duration, double dt, const char* label) {
        if (!std::isfinite(duration) || !std::isfinite(dt) || duration < 0 || dt <= 0)
            throw std::invalid_argument(std::string(label)+": times must be finite and dt positive");
        const double ratio = duration/dt;
        if (!std::isfinite(ratio) || ratio > std::numeric_limits<int>::max())
            throw std::invalid_argument(std::string(label)+": step count exceeds integer range");
        const auto steps = static_cast<int>(std::llround(ratio));
        if (!close(steps*dt, duration))
            throw std::invalid_argument(std::string(label)+": time must be an integer multiple of dt");
        return steps;
    }
    static CoupledCadence make(double atm_dt, double ocn_dt, double exchange) {
        const int atm = integral_steps(exchange, atm_dt, "driver.exchange_interval / erf.fixed_dt");
        const int ocn = integral_steps(exchange, ocn_dt, "driver.exchange_interval / remora.fixed_dt");
        if (atm < 1 || ocn < 1)
            throw std::invalid_argument("driver.exchange_interval must contain at least one step of each solver");
        return {atm_dt, ocn_dt, exchange, atm, ocn};
    }
    int ocean_steps_before(double stop_time) const {
        if (!std::isfinite(stop_time) || stop_time < 0 ||
            stop_time/ocean_dt > std::numeric_limits<int>::max())
            throw std::invalid_argument("stop_time must be finite, nonnegative, and within the ocean step range");
        const double raw = stop_time/exchange_interval;
        const auto nearest = static_cast<int>(std::llround(raw));
        const int intervals = close(nearest*exchange_interval, stop_time)
            ? nearest : static_cast<int>(std::floor(raw));
        return intervals*ocean_steps;
    }
    void require_boundary(int steps) const {
        if (steps < 0 || steps%ocean_steps != 0)
            throw std::invalid_argument("Coupled ocean-step target must end on a driver.exchange_interval boundary");
    }
};
} // namespace erf_remora
#endif
