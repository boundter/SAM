// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "test/harmonic_oscillator_ode.hpp"
#include "include/sam/system/rk4_system.hpp"

TEST_CASE("integrate simple system"){
  // Integrate correctly changes the position
  // Analytic solution is x=A*sin(omega*t + phi)
  // for x(0)=0 and \dot{x}(0)=1 and with omega=2
  // the solution is x = 0.5*sin(omega*t), \dot{x} = cos(omega*t)
  double omega = 2.;
  unsigned int N = 1;
  unsigned int dimension = 2;
  double dt = 0.01;
  unsigned int n = 10000;
  double t = dt*static_cast<double>(n);
  std::vector<double> initial_condition({0., 1.});
  sam::RK4System<HarmonicOscillatorODE> system(N, dimension, omega);
  system.SetPosition(initial_condition);

  SECTION("integrate without observer") {
    double t_0 = system.GetTime();
    system.Integrate(dt, n);
    std::vector<double> analytical = {1./omega*sin(omega*t), cos(omega*t)};
    std::vector<double> numerical = system.GetPosition();
    REQUIRE(numerical.size() == analytical.size());
    CHECK(numerical[0] == Approx(analytical[0]).margin(0.01));
    CHECK(numerical[1] == Approx(analytical[1]).margin(0.01));
    // Integrate increases time correctly
    double t_1 = system.GetTime();
    CHECK(t == Approx(t_1 - t_0).margin(0.000001));
  }
}