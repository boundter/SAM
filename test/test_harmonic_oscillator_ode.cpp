// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "test/harmonic_oscillator_ode.hpp"

typedef std::vector<double> state_type;

TEST_CASE("harmonic oscillator ODE works", "[test]") {
  double omega = 2.;
  HarmonicOscillatorODE ode(omega);
  REQUIRE(ode.omega_ == omega);
  state_type x = {2., 1.3};
  state_type dx(2);
  ode(x, dx, 0.);
  state_type analytic = {1.3, -8.};
  REQUIRE(dx.size() == analytic.size());
  CHECK(dx[0] == Approx(analytic[0]).margin(0.01));
  CHECK(dx[1] == Approx(analytic[1]).margin(0.01));
}

TEST_CASE("coupled harmonic oscillator ODE works", "[test]") {
  double omega_1 = 2.;
  double omega_2 = 3.;
  std::vector<double> omega({omega_1, omega_2});
  double coupling = 0.5;

  SECTION("double double constructor works") {
    CoupledHarmonicOscillatorODE ode(omega_1, omega_2, coupling);
    CHECK(ode.omega_1_ == omega_1);
    CHECK(ode.omega_2_ == omega_2);
    CHECK(ode.coupling_ == coupling);
  }
  CoupledHarmonicOscillatorODE ode(omega, coupling);
  CHECK(ode.omega_1_ == omega_1);
  CHECK(ode.omega_2_ == omega_2);
  CHECK(ode.coupling_ == coupling);
  state_type x = {2., 1.3, 1.5, 3};
  state_type dx(4);
  ode(x, dx, 0.);
  state_type analytic = {1.3, -7.15, 3, -14.35};
  REQUIRE(dx.size() == analytic.size());
  CHECK(dx[0] == Approx(analytic[0]).margin(0.01));
  CHECK(dx[1] == Approx(analytic[1]).margin(0.01));
  CHECK(dx[2] == Approx(analytic[2]).margin(0.01));
  CHECK(dx[3] == Approx(analytic[3]).margin(0.01));
}