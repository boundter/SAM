// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "test/harmonic_oscillator_ode.hpp"
#include "include/sam/analysis/period.hpp"
#include "include/sam/system/rk4_system.hpp"

TEST_CASE("single harmonic oscillator") {
  // omega = 1 -> T = 2*pi
  sam::PeriodParameters params;
  double omega = 1;
  unsigned int N = 1;
  unsigned int dimension = 2;
  double dt = 0.01;
  std::vector<double> initial({1., 0.});
  sam::RK4System<HarmonicOscillatorODE> system(N, dimension, omega);
  system.SetPosition(initial);
  double T = CalculatePeriod(system, dt,
                             [](std::vector<double> x) { return x[1] > 0; }, params);
  REQUIRE(T == Approx(2.*M_PI).margin(1e-5));
}

TEST_CASE("two uncoupled harmonic oscillators") {
  // omega = 1 -> T = 2*pi
  sam::PeriodParameters params;
  std::vector<double> omega({1, 2});
  double coupling = 0;
  unsigned int N = 2;
  unsigned int dimension = 2;
  double dt = 0.01;
  std::vector<double> initial({1., 0., 1., 0.});
  sam::RK4System<CoupledHarmonicOscillatorODE> system(N, dimension, omega[0],
                                                      omega[1], coupling);
  system.SetPosition(initial);

  SECTION("first oscillator") {
    double T = sam::CalculatePeriod(system, dt,
                               [](std::vector<double> x) { return x[1] > 0; }, params);
    REQUIRE(T == Approx(2.*M_PI).margin(1e-5));
  }

  SECTION("second oscillator") {
    params.n_osc = 1;
    double T = sam::CalculatePeriod(system, dt,
                               [](std::vector<double> x) { return x[3] > 0; }, params);
    REQUIRE(T == Approx(M_PI).margin(1e-5));
  }
}
