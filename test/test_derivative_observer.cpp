// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "harmonic_oscillator_ode.hpp"
#include "include/sam/observer/derivative_observer.hpp"
#include "include/sam/system/generic_system.hpp"

TEST_CASE("observer saves derivative") {
  std::vector<double> pos_1({1.5, 2.6});
  std::vector<double> pos_2({2.3, 3.8});
  double t_1 = 0.1;
  double t_2 = 0.2;
  std::vector<std::vector<double>> x;
  std::vector<double> t;
  double omega = 2;
  sam::GenericSystem<HarmonicOscillatorODE> system(1, 2, omega);
  sam::DerivativeObserver<sam::GenericSystem<HarmonicOscillatorODE>,
                          std::vector<double>> observer(system, x, t);
  system.SetPosition(pos_1);
  system.SetTime(t_1);
  observer(pos_1, t_1);
  system.SetPosition(pos_2);
  system.SetTime(t_2);
  observer(pos_2, t_2);
  REQUIRE(x.size() == 2);
  CHECK(x[0] == std::vector<double>({pos_1[1], -omega*omega*pos_1[0]}));
  CHECK(x[1] == std::vector<double>({pos_2[1], -omega*omega*pos_2[0]}));
  REQUIRE(t.size() == 2);
  CHECK(t[0] == t_1);
  CHECK(t[1] == t_2);
}
