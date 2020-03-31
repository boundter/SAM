// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <stdexcept>
#include <vector>

#include "test/catch.hpp"
#include "include/sam/analysis/savitzky_golay.hpp"

TEST_CASE("size and coefficients") {
  std::vector<double> x(5, 0);
  double step = 0.1;
  unsigned int points = 5;
  unsigned int derivatives = 4;

  SECTION("even number points throws") {
    points = 4;
    CHECK_THROWS_AS(sam::SavitzkyGolayFilter(step, x, points, derivatives),
                    std::invalid_argument);
  }

  SECTION("more derivatives than number of points throws") {
    derivatives = 6;
    CHECK_THROWS_AS(sam::SavitzkyGolayFilter(step, x, points, derivatives),
                    std::invalid_argument);
  }

  SECTION("size is correct") {
    std::vector<std::vector<double>> filter = sam::SavitzkyGolayFilter(step, x,
        points, derivatives);
    REQUIRE(filter.size() == derivatives);
    for (size_t i = 0; i < filter.size(); ++i) {
      REQUIRE(filter[i].size() == x.size() - points + 1);
    }
  }

  SECTION("coefficients for -2 for 5 points and 4 derivatives") {
    x[0] = 1;
    std::vector<double> y_m2(
        {-3./35, 1./12./step, 2./14.*2./(pow(step, 2)),
        -1./12.*6./pow(step, 3)});
    std::vector<std::vector<double>> filter = sam::SavitzkyGolayFilter(step, x,
    points, derivatives);
    for (size_t i = 0; i < derivatives; ++i) {
      CHECK(filter[i][0] == Approx(y_m2[i]).margin(1e-6));
    }
  }

  SECTION("coefficients for -1 for 5 points and 4 derivatives") {
    x[1] = 1;
    std::vector<double> y_m1(
        {12./35., -8./12./step, -1./14.*2./(pow(step, 2)),
        2./12.*6./pow(step, 3)});
    std::vector<std::vector<double>> filter = sam::SavitzkyGolayFilter(step, x,
    points, derivatives);
    for (size_t i = 0; i < derivatives; ++i) {
      CHECK(filter[i][0] == Approx(y_m1[i]).margin(1e-6));
    }
  }

  SECTION("coefficients for 0 for 5 points and 4 derivatives") {
    x[2] = 1;
    std::vector<double> y_0(
        {17./35., 0./12./step, -2./14.*2./(pow(step, 2)),
        0./12.*6./pow(step, 3)});
    std::vector<std::vector<double>> filter = sam::SavitzkyGolayFilter(step, x,
    points, derivatives);
    for (size_t i = 0; i < derivatives; ++i) {
      CHECK(filter[i][0] == Approx(y_0[i]).margin(1e-6));
    }
  }

  SECTION("coefficients for 1 for 5 points and 4 derivatives") {
    x[3] = 1;
    std::vector<double> y_1(
        {12./35., 8./12./step, -1./14.*2./(pow(step, 2)),
        -2./12.*6./pow(step, 3)});
    std::vector<std::vector<double>> filter = sam::SavitzkyGolayFilter(step, x,
    points, derivatives);
    for (size_t i = 0; i < derivatives; ++i) {
      CHECK(filter[i][0] == Approx(y_1[i]).margin(1e-6));
    }
  }

  SECTION("coefficients for 2 for 5 points and 4 derivatives") {
    x[4] = 1;
    std::vector<double> y_2(
        {-3./35., -1./12./step, 2./14.*2./(pow(step, 2)),
        1./12.*6./pow(step, 3)});
    std::vector<std::vector<double>> filter = sam::SavitzkyGolayFilter(step, x,
    points, derivatives);
    for (size_t i = 0; i < derivatives; ++i) {
      CHECK(filter[i][0] == Approx(y_2[i]).margin(1e-6));
    }
  }
}