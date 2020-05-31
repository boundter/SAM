// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#include <vector>
#include <utility>

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

  SECTION("resizing works") {
    std::vector<unsigned int> new_sizes({4, 3});
    std::vector<unsigned int> node_indices({0, 8, 14});
    system.Resize(new_sizes);
    std::vector<unsigned int> indices = system.GetNodeIndices();
    REQUIRE(indices == node_indices);
    std::vector<double> pos = system.GetPosition();
    REQUIRE(pos.size() == 14);
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

  SECTION("getting nodes spherical in 1d") {
    std::vector<std::vector<double>> nodes = system.GetNodesSpherical();
    REQUIRE(nodes.size() == node_size.size());
    REQUIRE(nodes[0].size() == 1);
    CHECK(nodes[0][0] == Approx(x[0]).margin(0.01));
    REQUIRE(nodes[1].size() == 2);
    CHECK(nodes[1][0] == Approx(x[1]).margin(0.01));
    CHECK(nodes[1][1] == Approx(x[2]).margin(0.01));
  }
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

  SECTION("get dimension") {
    std::pair<unsigned int, unsigned int> dim = system.GetDimension();
    CHECK(dim.first == 3);
    CHECK(dim.second == 2);
  }

  SECTION("getting nodes spherical in 2d") {
    std::vector<std::vector<double>> nodes = system.GetNodesSpherical();
    REQUIRE(nodes.size() == node_size.size());
    REQUIRE(nodes[0].size() == 2);
    CHECK(nodes[0][0] == Approx(analytical[0]).margin(0.1));
    CHECK(nodes[0][1] == Approx(analytical[1]).margin(0.1));
    REQUIRE(nodes[1].size() == 4);
    CHECK(nodes[1][0] == Approx(analytical[2]).margin(0.1));
    CHECK(nodes[1][1] == Approx(analytical[3]).margin(0.1));
    CHECK(nodes[1][2] == Approx(analytical[4]).margin(0.1));
    CHECK(nodes[1][3] == Approx(analytical[5]).margin(0.1));
  }
}

TEST_CASE("Derivative", "[generic_network]") {
  double omega = 2.;
  std::vector<unsigned int> node_size({1, 2});
  unsigned int d = 2;
  sam::GenericNetwork<HarmonicOscillatorODE> system(node_size, d, omega);
  double t = 5.;
  std::vector<double> x({1., 2., 3., 4., 5., 6.});
  std::vector<double> derivative({2., -4., 4., -12., 6., -20.});
  system.SetTime(t);
  system.SetPosition(x);

  SECTION("derivative in nodes") {
    std::vector<std::vector<double>> nodes = system.GetDerivativeNodes();
    REQUIRE(nodes.size() == 2);
    REQUIRE(nodes[0].size() == 2);
    CHECK(nodes[0][0] == Approx(derivative[0]).margin(0.0001));
    CHECK(nodes[0][1] == Approx(derivative[1]).margin(0.0001));
    REQUIRE(nodes[1].size() == 4);
    CHECK(nodes[1][0] == Approx(derivative[2]).margin(0.0001));
    CHECK(nodes[1][1] == Approx(derivative[3]).margin(0.0001));
    CHECK(nodes[1][2] == Approx(derivative[4]).margin(0.0001));
    CHECK(nodes[1][3] == Approx(derivative[5]).margin(0.0001));
    }
}

TEST_CASE("Copying", "[generic_network]") {
  double omega = 2.;
  std::vector<unsigned int> node_size({1, 2});
  unsigned int d = 2;
  sam::GenericNetwork<HarmonicOscillatorODE> system(node_size, d, omega);
  double t = 5.;
  std::vector<double> x({1., 2., 3., 4., 5., 6.});
  std::vector<double> derivative({2., -4., 4., -12., 6., -20.});
  system.SetTime(t);
  system.SetPosition(x);
  double system_time = system.GetTime();
  std::vector<double> system_pos = system.GetPosition();
  std::vector<double> system_derivative = system.GetDerivative();
  CHECK(system_time == Approx(t).margin(0.0001));
  REQUIRE(system_pos.size() == 6);
  CHECK(system_pos[0] == Approx(x[0]).margin(0.0001));
  CHECK(system_pos[1] == Approx(x[1]).margin(0.0001));
  CHECK(system_pos[2] == Approx(x[2]).margin(0.0001));
  CHECK(system_pos[3] == Approx(x[3]).margin(0.0001));
  CHECK(system_pos[4] == Approx(x[4]).margin(0.0001));
  CHECK(system_pos[5] == Approx(x[5]).margin(0.0001));
  REQUIRE(system_derivative.size() == 6);
  CHECK(system_derivative[0] == Approx(derivative[0]).margin(0.0001));
  CHECK(system_derivative[1] == Approx(derivative[1]).margin(0.0001));
  CHECK(system_derivative[2] == Approx(derivative[2]).margin(0.0001));
  CHECK(system_derivative[3] == Approx(derivative[3]).margin(0.0001));
  CHECK(system_derivative[4] == Approx(derivative[4]).margin(0.0001));
  CHECK(system_derivative[5] == Approx(derivative[5]).margin(0.0001));

  SECTION("copy works") {
    sam::GenericNetwork<HarmonicOscillatorODE> copy_system = system;
    double copy_time = copy_system.GetTime();
    std::vector<double> copy_pos = copy_system.GetPosition();
    std::vector<double> copy_derivative = copy_system.GetDerivative();
    CHECK(copy_time == Approx(t).margin(0.0001));
    REQUIRE(copy_pos.size() == 6);
    CHECK(copy_pos[0] == Approx(x[0]).margin(0.0001));
    CHECK(copy_pos[1] == Approx(x[1]).margin(0.0001));
    CHECK(copy_pos[2] == Approx(x[2]).margin(0.0001));
    CHECK(copy_pos[3] == Approx(x[3]).margin(0.0001));
    CHECK(copy_pos[4] == Approx(x[4]).margin(0.0001));
    CHECK(copy_pos[5] == Approx(x[5]).margin(0.0001));
    REQUIRE(copy_derivative.size() == 6);
    CHECK(copy_derivative[0] == Approx(derivative[0]).margin(0.0001));
    CHECK(copy_derivative[1] == Approx(derivative[1]).margin(0.0001));
    CHECK(copy_derivative[2] == Approx(derivative[2]).margin(0.0001));
    CHECK(copy_derivative[3] == Approx(derivative[3]).margin(0.0001));
    CHECK(copy_derivative[4] == Approx(derivative[4]).margin(0.0001));
    CHECK(copy_derivative[5] == Approx(derivative[5]).margin(0.0001));
  }

  SECTION("deep copy independent of original") {
    sam::GenericNetwork<HarmonicOscillatorODE> copy_system = system;
    double new_time = 6.;
    std::vector<double> new_pos({6., 5., 4., 3., 2., 1.});
    std::vector<double> new_derivative({5., -24., 3., -16., 1., -8.});
    copy_system.SetTime(new_time);
    copy_system.SetPosition(new_pos);

    double copy_time = copy_system.GetTime();
    std::vector<double> copy_pos = copy_system.GetPosition();
    std::vector<double> copy_derivative = copy_system.GetDerivative();
    CHECK(copy_time == Approx(new_time).margin(0.0001));
    REQUIRE(copy_pos.size() == 6);
    CHECK(copy_pos[0] == Approx(new_pos[0]).margin(0.0001));
    CHECK(copy_pos[1] == Approx(new_pos[1]).margin(0.0001));
    CHECK(copy_pos[2] == Approx(new_pos[2]).margin(0.0001));
    CHECK(copy_pos[3] == Approx(new_pos[3]).margin(0.0001));
    CHECK(copy_pos[4] == Approx(new_pos[4]).margin(0.0001));
    CHECK(copy_pos[5] == Approx(new_pos[5]).margin(0.0001));
    REQUIRE(copy_derivative.size() == 6);
    CHECK(copy_derivative[0] == Approx(new_derivative[0]).margin(0.0001));
    CHECK(copy_derivative[1] == Approx(new_derivative[1]).margin(0.0001));
    CHECK(copy_derivative[2] == Approx(new_derivative[2]).margin(0.0001));
    CHECK(copy_derivative[3] == Approx(new_derivative[3]).margin(0.0001));
    CHECK(copy_derivative[4] == Approx(new_derivative[4]).margin(0.0001));
    CHECK(copy_derivative[5] == Approx(new_derivative[5]).margin(0.0001));

    double system_time = system.GetTime();
    std::vector<double> system_pos = system.GetPosition();
    std::vector<double> system_derivative = system.GetDerivative();
    CHECK(system_time == Approx(t).margin(0.0001));
    REQUIRE(system_pos.size() == 6);
    CHECK(system_pos[0] == Approx(x[0]).margin(0.0001));
    CHECK(system_pos[1] == Approx(x[1]).margin(0.0001));
    CHECK(system_pos[2] == Approx(x[2]).margin(0.0001));
    CHECK(system_pos[3] == Approx(x[3]).margin(0.0001));
    CHECK(system_pos[4] == Approx(x[4]).margin(0.0001));
    CHECK(system_pos[5] == Approx(x[5]).margin(0.0001));
    REQUIRE(system_derivative.size() == 6);
    CHECK(system_derivative[0] == Approx(derivative[0]).margin(0.0001));
    CHECK(system_derivative[1] == Approx(derivative[1]).margin(0.0001));
    CHECK(system_derivative[2] == Approx(derivative[2]).margin(0.0001));
    CHECK(system_derivative[3] == Approx(derivative[3]).margin(0.0001));
    CHECK(system_derivative[4] == Approx(derivative[4]).margin(0.0001));
    CHECK(system_derivative[5] == Approx(derivative[5]).margin(0.0001));
  }

  SECTION("ODE independent of original system") {
    sam::GenericNetwork<HarmonicOscillatorODE> copy_system = system;
    double new_omega = 3.;
    copy_system.SetParameters(new_omega);
    std::vector<double> new_derivative({2., -9., 4., -27., 6., -45.});
    double copy_time = copy_system.GetTime();
    std::vector<double> copy_pos = copy_system.GetPosition();
    std::vector<double> copy_derivative = copy_system.GetDerivative();

    CHECK(copy_time == Approx(t).margin(0.0001));
    REQUIRE(copy_pos.size() == 6);
    CHECK(copy_pos[0] == Approx(x[0]).margin(0.0001));
    CHECK(copy_pos[1] == Approx(x[1]).margin(0.0001));
    CHECK(copy_pos[2] == Approx(x[2]).margin(0.0001));
    CHECK(copy_pos[3] == Approx(x[3]).margin(0.0001));
    CHECK(copy_pos[4] == Approx(x[4]).margin(0.0001));
    CHECK(copy_pos[5] == Approx(x[5]).margin(0.0001));
    REQUIRE(copy_derivative.size() == 6);
    CHECK(copy_derivative[0] == Approx(new_derivative[0]).margin(0.0001));
    CHECK(copy_derivative[1] == Approx(new_derivative[1]).margin(0.0001));
    CHECK(copy_derivative[2] == Approx(new_derivative[2]).margin(0.0001));
    CHECK(copy_derivative[3] == Approx(new_derivative[3]).margin(0.0001));
    CHECK(copy_derivative[4] == Approx(new_derivative[4]).margin(0.0001));
    CHECK(copy_derivative[5] == Approx(new_derivative[5]).margin(0.0001));

    CHECK(system_time == Approx(t).margin(0.0001));
    REQUIRE(system_pos.size() == 6);
    CHECK(system_pos[0] == Approx(x[0]).margin(0.0001));
    CHECK(system_pos[1] == Approx(x[1]).margin(0.0001));
    CHECK(system_pos[2] == Approx(x[2]).margin(0.0001));
    CHECK(system_pos[3] == Approx(x[3]).margin(0.0001));
    CHECK(system_pos[4] == Approx(x[4]).margin(0.0001));
    CHECK(system_pos[5] == Approx(x[5]).margin(0.0001));
    REQUIRE(system_derivative.size() == 6);
    CHECK(system_derivative[0] == Approx(derivative[0]).margin(0.0001));
    CHECK(system_derivative[1] == Approx(derivative[1]).margin(0.0001));
    CHECK(system_derivative[2] == Approx(derivative[2]).margin(0.0001));
    CHECK(system_derivative[3] == Approx(derivative[3]).margin(0.0001));
    CHECK(system_derivative[4] == Approx(derivative[4]).margin(0.0001));
    CHECK(system_derivative[5] == Approx(derivative[5]).margin(0.0001));
  }
}
