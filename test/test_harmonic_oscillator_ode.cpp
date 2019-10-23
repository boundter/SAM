// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "test/harmonic_oscillator_ode.hpp"

namespace mt = Catch::Matchers;

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
  CHECK_THAT(dx[0],  mt::WithinAbs(analytic[0], 0.01));
  CHECK_THAT(dx[1],  mt::WithinAbs(analytic[1], 0.01));
}
