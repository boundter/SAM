// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "include/sam/observer/position_observer.hpp"

TEST_CASE("PositionObserver saves position") {
  std::vector<double> pos_1({1.5, 2.6});
  std::vector<double> pos_2({2.3, 3.8});
  double t_1 = 0.1;
  double t_2 = 0.2;
  std::vector<std::vector<double>> x;
  std::vector<double> t;
  sam::PositionObserver<std::vector<double>> observer(x, t);
  observer(pos_1, t_1);
  observer(pos_2, t_2);
  REQUIRE(x.size() == 2);
  CHECK(x[0] == pos_1);
  CHECK(x[1] == pos_2);
  REQUIRE(t.size() == 2);
  CHECK(t[0] == t_1);
  CHECK(t[1] == t_2);
}
