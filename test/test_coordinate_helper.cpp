// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "include/sam/helper/coordinate_helper.hpp"

TEST_CASE("cartesian to spherical in 2d with positive phase from iterator") {
  std::vector<double> x = {2., 3.};
  std::vector<double> analytical = {3.606, 0.9828};
  std::vector<double> spherical =
      sam::CartesianToSpherical<std::vector<double>>(x.begin(), x.end());
  REQUIRE(spherical.size() == x.size());
  CHECK(spherical[0] == Approx(analytical[0]).margin(0.1));
  CHECK(spherical[1] == Approx(analytical[1]).margin(0.1));
}

TEST_CASE("cartesian to spherical in 2d with negative phase from iterator") {
  std::vector<double> x = {2., -3.};
  std::vector<double> analytical = {3.606, 5.30};
  std::vector<double> spherical =
      sam::CartesianToSpherical<std::vector<double>>(x.begin(), x.end());
  REQUIRE(spherical.size() == x.size());
  CHECK(spherical[0] == Approx(analytical[0]).margin(0.1));
  CHECK(spherical[1] == Approx(analytical[1]).margin(0.1));
}

TEST_CASE("cartesian to spherical in 3d from iterator") {
  std::vector<double> x = {6., 2., -1.};
  std::vector<double> analytical = {6.403, 0.3567, 5.819};
  std::vector<double> spherical =
      sam::CartesianToSpherical<std::vector<double>>(x.begin(), x.end());
  REQUIRE(spherical.size() == x.size());
  CHECK(spherical[0] == Approx(analytical[0]).margin(0.01));
  CHECK(spherical[1] == Approx(analytical[1]).margin(0.01));
  CHECK(spherical[2] == Approx(analytical[2]).margin(0.01));
}

TEST_CASE("cartesian to spherical in 2d with positive phase from vector") {
  std::vector<double> x = {2., 3.};
  std::vector<double> analytical = {3.606, 0.9828};
  std::vector<double> spherical = sam::CartesianToSpherical(x);
  REQUIRE(spherical.size() == x.size());
  CHECK(spherical[0] == Approx(analytical[0]).margin(0.1));
  CHECK(spherical[1] == Approx(analytical[1]).margin(0.1));
}

TEST_CASE("cartesian to spherical in 2d with negative phase from vector") {
  std::vector<double> x = {2., -3.};
  std::vector<double> analytical = {3.606, 5.30};
  std::vector<double> spherical = sam::CartesianToSpherical(x);
  REQUIRE(spherical.size() == x.size());
  CHECK(spherical[0] == Approx(analytical[0]).margin(0.1));
  CHECK(spherical[1] == Approx(analytical[1]).margin(0.1));
}

TEST_CASE("cartesian to spherical in 3d from vector") {
  std::vector<double> x = {6., 2., -1.};
  std::vector<double> analytical = {6.403, 0.3567, 5.819};
  std::vector<double> spherical = sam::CartesianToSpherical(x);
  REQUIRE(spherical.size() == x.size());
  CHECK(spherical[0] == Approx(analytical[0]).margin(0.01));
  CHECK(spherical[1] == Approx(analytical[1]).margin(0.01));
  CHECK(spherical[2] == Approx(analytical[2]).margin(0.01));
}

TEST_CASE("cartesian to spherical for subvectors with iterators") {
  // states of different dimension (2, 3); (6, 2, -1)
  std::vector<double> x = {2, 3, 6, 2, -1};
  std::vector<double> analytical_2 = {3.606, 0.9828};
  std::vector<double> analytical_3 = {6.403, 0.3567, 5.819};
  std::vector<double> spherical_2 =
      sam::CartesianToSpherical<std::vector<double>>(x.begin(), x.begin()+2);
  REQUIRE(spherical_2.size() == analytical_2.size());
  CHECK(spherical_2[0] == Approx(analytical_2[0]).margin(0.1));
  CHECK(spherical_2[1] == Approx(analytical_2[1]).margin(0.1));
  std::vector<double> spherical_3 =
    sam::CartesianToSpherical<std::vector<double>>(x.begin()+2, x.end());
  REQUIRE(spherical_3.size() == analytical_3.size());
  CHECK(spherical_3[0] == Approx(analytical_3[0]).margin(0.1));
  CHECK(spherical_3[1] == Approx(analytical_3[1]).margin(0.1));
  CHECK(spherical_3[2] == Approx(analytical_3[2]).margin(0.1));

}