// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>

#include "test/catch.hpp"
#include "test/harmonic_oscillator_ode.hpp"
#include "include/sam/system/generic_network.hpp"

TEST_CASE("getting and setting position", "[generic_network]") {
  std::vector<unsigned int> node_sizes({3, 2});
  double omega = 1;
  sam::GenericNetwork<HarmonicOscillatorODE> system(node_sizes, 2, omega);
  std::vector<double> x({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  CHECK_NOTHROW(system.SetPosition(x));

  SECTION("setting position works") {
    std::vector<double> state = system.GetPosition();
    REQUIRE(state.size() == x.size());
    CHECK(state[0] == Approx(x[0]).margin(0.01));
    CHECK(state[1] == Approx(x[1]).margin(0.01));
    CHECK(state[2] == Approx(x[2]).margin(0.01));
    CHECK(state[3] == Approx(x[3]).margin(0.01));
    CHECK(state[4] == Approx(x[4]).margin(0.01));
    CHECK(state[5] == Approx(x[5]).margin(0.01));
    CHECK(state[6] == Approx(x[6]).margin(0.01));
    CHECK(state[7] == Approx(x[7]).margin(0.01));
    CHECK(state[8] == Approx(x[8]).margin(0.01));
    CHECK(state[9] == Approx(x[9]).margin(0.01));
  }

  SECTION("node indices are calculated correctly") {
    std::vector<unsigned int> node_indices({0, 6, 10});
    std::vector<unsigned int> indices = system.GetNodeIndices();
    REQUIRE(indices.size() == node_indices.size());
    CHECK(indices[0] == node_indices[0]);
    CHECK(indices[1] == node_indices[1]);
    CHECK(indices[2] == node_indices[2]);
  }

  SECTION("nodes are returned correctly") {
    std::vector<std::vector<double>> nodes = system.GetNodes();
    REQUIRE(nodes.size() == 2);
    REQUIRE(nodes[0].size() == 6);
    REQUIRE(nodes[1].size() == 4);
    CHECK(nodes[0][0] == Approx(x[0]).margin(0.01));
    CHECK(nodes[0][1] == Approx(x[1]).margin(0.01));
    CHECK(nodes[0][2] == Approx(x[2]).margin(0.01));
    CHECK(nodes[0][3] == Approx(x[3]).margin(0.01));
    CHECK(nodes[0][4] == Approx(x[4]).margin(0.01));
    CHECK(nodes[0][5] == Approx(x[5]).margin(0.01));
    CHECK(nodes[1][0] == Approx(x[6]).margin(0.01));
    CHECK(nodes[1][1] == Approx(x[7]).margin(0.01));
    CHECK(nodes[1][2] == Approx(x[8]).margin(0.01));
    CHECK(nodes[1][3] == Approx(x[9]).margin(0.01));
  }
}

TEST_CASE("getting position in spherical in 1d", "[generic_network]") {
  // use HarmonicOscillatorODE as a dummy ODE
  double omega = 1;
  std::vector<unsigned int> node_size = {1, 2};
  unsigned int d = 1;
  sam::GenericNetwork<HarmonicOscillatorODE> system(node_size, d, omega);
  std::vector<double> x({1., 2., 3.});
  system.SetPosition(x);
  std::vector<double> spherical = system.GetPositionSpherical();
  REQUIRE(spherical.size() == x.size());
  CHECK(spherical[0] == Approx(x[0]).margin(0.01));
  CHECK(spherical[1] == Approx(x[1]).margin(0.01));
  CHECK(spherical[2] == Approx(x[2]).margin(0.01));
}

TEST_CASE("getting position in spherical in 2d", "[generic_network]") {
  // use HarmonicOscillatorODE as a dummy ODE
  double omega = 1;
  std::vector<unsigned int> node_size({1, 2});
  unsigned int d = 2;
  sam::GenericNetwork<HarmonicOscillatorODE> system(node_size, d, omega);
  std::vector<double> x({1., 2., 3., 4., 5., 6.});
  std::vector<double> analytical({2.236, 1.107, 5, 0.927, 7.81, 0.876});
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