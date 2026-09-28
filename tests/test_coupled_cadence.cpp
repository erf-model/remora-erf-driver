#include "CoupledCadence.H"
#include <cstdlib>
#include <iostream>
#include <limits>
using erf_remora::CoupledCadence;
void require(bool condition) { if (!condition) std::abort(); }
template<class F> void rejects(F f) {
    try { f(); } catch (const std::invalid_argument&) { return; }
    std::abort();
}
int main() {
    for (double exchange : {300., 600.}) {
        const auto c = CoupledCadence::make(10, 15, exchange);
        require(c.atmosphere_steps == exchange/10 && c.ocean_steps == exchange/15);
        require(c.ocean_steps_before(2592000) == 172800);
        require(c.ocean_steps_before(exchange-1) == 0);
        require(c.ocean_steps_before(exchange+1) == c.ocean_steps);
        c.require_boundary(172800);
        c.require_boundary(0);
        rejects([&] { c.require_boundary(c.ocean_steps-1); });
        rejects([&] { CoupledCadence::integral_steps(15, exchange, "restart"); });
        require(CoupledCadence::integral_steps(3600, exchange, "restart") == 3600/exchange);
    }
    const auto c = CoupledCadence::make(.1, .3, .3);
    require(c.atmosphere_steps == 3 && c.ocean_steps == 1 && c.ocean_steps_before(.9) == 3);
    rejects([] { CoupledCadence::make(10, 15, 15); });
    rejects([] { CoupledCadence::make(10, 15, 301); });
    rejects([] { CoupledCadence::make(0, 15, 300); });
    rejects([] { CoupledCadence::make(10, -15, 300); });
    rejects([] { CoupledCadence::make(10, 15, 0); });
    rejects([] { CoupledCadence::make(1.e-30, 15, 300); });
    rejects([] { CoupledCadence::make(std::numeric_limits<double>::quiet_NaN(), 15, 300); });
    rejects([] { CoupledCadence::make(10, 15, std::numeric_limits<double>::infinity()); });
    std::cout << "cadence: month counts, 300/600 s, decimal dts, stop/restart boundaries, invalid inputs PASS\n";
}
