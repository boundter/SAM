// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "test/harmonic_oscillator_ode.hpp"

typedef std::vector<double> state_type;

TEST_CASE("Harmonic oscillator ODE works", "[test]") {
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
