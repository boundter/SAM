// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <stdexcept>
#include <vector>

#include "test/catch.hpp"
#include "test/harmonic_oscillator_ode.hpp"
#include "include/sam/system/generic_system.hpp"

TEST_CASE("simple system", "[generic_system]") {
  double omega = 2.;
  unsigned int N = 1;
  unsigned int dimension = 2;
  sam::GenericSystem<HarmonicOscillatorODE> system(N, dimension, omega);

  SECTION("state initializes to 0") {
    std::vector<double> position = system.GetPosition();
    REQUIRE(position.size() == dimension);
    CHECK(position[0] == Approx(0.).margin(0.01));
    CHECK(position[1] == Approx(0.).margin(0.01));
  }

  SECTION("set position changes internal position") {
    std::vector<double> new_state({0.5, 0.1});
    system.SetPosition(new_state);
    std::vector<double> position = system.GetPosition();
    REQUIRE(position.size() == dimension);
    CHECK(position[0] == Approx(0.5).margin(0.01));
    CHECK(position[1] == Approx(0.1).margin(0.01));
  }

  SECTION("cannot set position with wrong size") {
    std::vector<double> too_long = {0.3, 0.1, 6.0};
    std::vector<double> too_short = {0.1};
    CHECK_THROWS_AS(system.SetPosition(too_long), std::length_error);
    CHECK_THROWS_AS(system.SetPosition(too_short), std::length_error);
  }

  SECTION("return the derivative without integrating") {
    std::vector<double> initial {0.5, 0.1};
    system.SetPosition(initial);
    std::vector<double> derivative = system.GetDerivative();
    std::vector<double> position = system.GetPosition();
    REQUIRE(derivative.size() == dimension);
    CHECK(derivative[0] == Approx(0.1).margin(0.01));
    CHECK(derivative[1] == Approx(-2.).margin(0.01));
    REQUIRE(position.size() == initial.size());
    CHECK(position[0] == Approx(initial[0]).margin(0.01));
    CHECK(position[1] == Approx(initial[1]).margin(0.01));
  }

  SECTION("can resize the system") {
    system.Resize(2);
    // check if state has correct size
    REQUIRE(system.GetPosition().size() == 4);
    // check if new state can be set
    std::vector<double> new_state = {1., 2., 3., 4.};
    REQUIRE_NOTHROW(system.SetPosition(new_state));
    std::vector<double> position = system.GetPosition();
    REQUIRE(position.size() == new_state.size());
    CHECK(position[0] == Approx(new_state[0]).margin(0.01));
    CHECK(position[1] == Approx(new_state[1]).margin(0.01));
    CHECK(position[2] == Approx(new_state[2]).margin(0.01));
    CHECK(position[3] == Approx(new_state[3]).margin(0.01));
    // check false lengths
    std::vector<double> too_short(3);
    std::vector<double> too_long(5);
    CHECK_THROWS_AS(system.SetPosition(too_short), std::length_error);
    CHECK_THROWS_AS(system.SetPosition(too_long), std::length_error);
  }

  SECTION("change parameters") {
    double new_omega = 8.3;
    system.SetParameters(new_omega);
    // indirect check using the derivative
      std::vector<double> initial {0.5, 0.1};
    system.SetPosition(initial);
    std::vector<double> derivative = system.GetDerivative();
    std::vector<double> position = system.GetPosition();
    REQUIRE(derivative.size() == dimension);
    CHECK(derivative[0] == Approx(0.1).margin(0.01));
    CHECK(derivative[1] == Approx(-34.445).margin(0.001));
    REQUIRE(position.size() == initial.size());
    CHECK(position[0] == Approx(initial[0]).margin(0.01));
    CHECK(position[1] == Approx(initial[1]).margin(0.01));
  }
}


TEST_CASE("ODE with multiple parameters", "[generic_system]") {
  double omega_1 = 2.;
  double omega_2 = 3.;
  std::vector<double> omega({omega_1, omega_2});
  double coupling = 0.5;
  unsigned int N = 2;
  unsigned int dimension = 2;

  SECTION("multiple parameters of same type") {
    sam::GenericSystem<CoupledHarmonicOscillatorODE> system(
        N, dimension, omega_1, omega_2, coupling);
    std::vector<double> position = {1., 2., 4., 5.};
    system.SetPosition(position);
    std::vector<double> derivative = system.GetDerivative();
    REQUIRE(derivative.size() == 4);
    CHECK(derivative[0] == Approx(position[1]).margin(0.01));
    CHECK(derivative[1] == Approx(-2.5).margin(0.01));
    CHECK(derivative[2] == Approx(position[3]).margin(0.01));
    CHECK(derivative[3] == Approx(-37.5).margin(0.01));
  }

  SECTION("multiple parameters of different types") {
    sam::GenericSystem<CoupledHarmonicOscillatorODE> system(
        N, dimension, omega, coupling);
    std::vector<double> position = {1., 2., 4., 5.};
    system.SetPosition(position);
    std::vector<double> derivative = system.GetDerivative();
    REQUIRE(derivative.size() == 4);
    CHECK(derivative[0] == Approx(position[1]).margin(0.01));
    CHECK(derivative[1] == Approx(-2.5).margin(0.01));
    CHECK(derivative[2] == Approx(position[3]).margin(0.01));
    CHECK(derivative[3] == Approx(-37.5).margin(0.01));
  }

  SECTION("change parameters") {
    std::vector<double> new_omega({4., 6.});
    double new_coupling = 0.;
    sam::GenericSystem<CoupledHarmonicOscillatorODE> system(
        N, dimension, omega, coupling);
    REQUIRE_NOTHROW(system.SetParameters(new_omega, new_coupling));
    std::vector<double> position = {1., 2., 4., 5.};
    system.SetPosition(position);
    std::vector<double> derivative = system.GetDerivative();
    REQUIRE(derivative.size() == 4);
    CHECK(derivative[0] == Approx(position[1]).margin(0.01));
    CHECK(derivative[1] == Approx(-16).margin(0.01));
    CHECK(derivative[2] == Approx(position[3]).margin(0.01));
    CHECK(derivative[3] == Approx(-144).margin(0.01));
  }
}

TEST_CASE("position in spherical coordinates 1d", "[generic_system]") {
  // use HarmonicOscillator as dummy ODE but with dimension 1
  double omega = 1;
  unsigned int N = 3;
  unsigned int d = 1;
  sam::GenericSystem<HarmonicOscillatorODE> system(N, d, omega);
  std::vector<double> x = {1., 2., 3.};
  system.SetPosition(x);
  std::vector<double> spherical = system.GetPositionSpherical();
  REQUIRE(spherical.size() == x.size());
  CHECK(spherical[0] == Approx(x[0]).margin(0.01));
  CHECK(spherical[1] == Approx(x[1]).margin(0.01));
  CHECK(spherical[2] == Approx(x[2]).margin(0.01));
}

TEST_CASE("position in spherical coordinates 2d", "[generic_system]") {
  // use HarmonicOscillatorODE as a dummy ODE with dimensionality 2
  double omega = 1;
  unsigned int N = 3;
  unsigned int d = 2;
  sam::GenericSystem<HarmonicOscillatorODE> system(N, d, omega);
  std::vector<double> x = {1., 2., 3., 4., 5., -6.};
  std::vector<double> analytical = {2.236, 1.107, 5, 0.927, 7.81, 5.407};
  system.SetPosition(x);
  std::vector<double> spherical = system.GetPositionSpherical();
  REQUIRE(spherical.size() == analytical.size());
  CHECK(spherical[0] == Approx(analytical[0]).margin(0.1));
  CHECK(spherical[1] == Approx(analytical[1]).margin(0.1));
  CHECK(spherical[2] == Approx(analytical[2]).margin(0.1));
  CHECK(spherical[3] == Approx(analytical[3]).margin(0.1));
  CHECK(spherical[4] == Approx(analytical[4]).margin(0.1));
  CHECK(spherical[5] == Approx(analytical[5]).margin(0.1));
}
