// Copyright 2020 Erik Teichman <kontakt.teichmann@gmail.com>

#include <cmath>
#include <vector>
#include <utility>

#include "test/catch.hpp"
#include "include/sam/system/rk4_system.hpp"
#include "include/sam/analysis/period.hpp"
#include "include/sam/analysis/phase.hpp"

typedef std::vector<double> state_type;

class StuartLandauODE {
 public:
  explicit StuartLandauODE(double alpha): alpha_(alpha) {}

  void operator()(const state_type& x, state_type& dx, double t) {
    double R2 = x[0]*x[0] + x[1]*x[1];
    dx[0] = x[0] - x[1] - R2*(x[0] - alpha_*x[1]);
    dx[1] = x[1] + x[0] - R2*(x[1] + alpha_*x[0]);
  }

 private:
  double alpha_;
};

class StuartLandauForcedODE {
 public:
  explicit StuartLandauForcedODE(double alpha, double eps, double nu)
      : alpha_(alpha), eps_(eps), nu_(nu) { }

  void operator()(const state_type& x, state_type& dx, double t) {
    double R2 = x[0]*x[0] + x[1]*x[1];
    dx[0] = x[0] - x[1] - R2*(x[0] - alpha_*x[1]) + eps_*std::cos(nu_*t);
    dx[1] = x[1] + x[0] - R2*(x[1] + alpha_*x[0]);
  }

 private:
  double alpha_;
  double eps_;
  double nu_;
};

class StuartLandauLinearizedODE {
 public:
  explicit StuartLandauLinearizedODE(double alpha) : alpha_(alpha) {}

  void operator()(const state_type& x, state_type& dx, double t) {
    double R2 = x[0]*x[0] + x[1]*x[1];
    dx[0] = x[0] - x[1] - R2*(x[0] - alpha_*x[1]);
    dx[1] = x[1] + x[0] - R2*(x[1] + alpha_*x[0]);
    dx[2] = x[2] - x[3] - 3*x[0]*x[0]*x[2] - x[1]*x[1]*x[2]
            - 2*x[0]*x[1]*x[3] + 2*alpha_*x[0]*x[1]*x[2]
            + alpha_*x[0]*x[0]*x[3] + 3*alpha_*x[1]*x[1]*x[3];
    dx[3] = x[2] + x[3] - 2*x[0]*x[1]*x[2] - x[0]*x[0]*x[3] - 3*x[1]*x[1]*x[3]
            - 3*alpha_*x[0]*x[0]*x[2] - alpha_*x[1]*x[1]*x[2]
            - 2*alpha_*x[0]*x[1]*x[3];
  }

 private:
  double alpha_;
};

bool CrossingCondition(state_type x) {
  return x[1] > 0;
}

double AnalyticPhase(state_type x, double alpha) {
  double phi = std::atan2(x[1], x[0])
               - alpha*std::log(std::sqrt(x[0]*x[0] + x[1]*x[1]));
  phi -= M_PI/2.;  // because we take the y-axis for the crossing
  phi = phi < 0 ? phi += 2*M_PI : phi;
  return phi;
}

double AnalyticFrequency(state_type x, state_type dx, double alpha) {
  double denominator = x[0]*x[0] + x[1]*x[1];
  double numerator = x[0]*dx[1] - x[1]*dx[0]
                     - alpha*(x[0]*dx[0] + x[1]*dx[1]);
  return numerator/denominator;
}

TEST_CASE("phase for perturbed Stuart-Landau oscillator") {
  const double alpha = 0.1;
  const double nu = 1;
  const double eps = 0.2;
  const double dt = 0.01;
  const unsigned int n_trans = 1e4;

  state_type initial({1, 1});

  sam::RK4System<StuartLandauODE> system(1, 2, alpha);
  system.SetPosition(initial);
  system.Integrate(dt, n_trans);
  double T = sam::CalculatePeriod(system, 0, 0, dt, CrossingCondition);

  SECTION("phase on limit cycle") {
    for (unsigned int i = 0; i < 10; ++i) {
      system.Integrate(dt, 70);
      double phase_analytic = AnalyticPhase(system.GetPosition(), alpha);
      double phase = sam::PhaseOnLimitCycle(system, T, dt, 0, 0,
                                            CrossingCondition);
      REQUIRE(phase == Approx(phase_analytic).margin(0.001));
    }
  }

  SECTION("phase of perturbed system") {
    sam::RK4System<StuartLandauForcedODE> forced_system(1, 2, alpha, eps, nu);
    forced_system.SetPosition(initial);
    forced_system.Integrate(dt, n_trans);
    for (unsigned int i = 0; i < 10; ++i) {
      forced_system.Integrate(dt, 70);
      state_type pos = forced_system.GetPosition();
      double phase_analytic = AnalyticPhase(pos, alpha);
      double phase = sam::FindPhase(pos, T, system, 0, 0, CrossingCondition);
      REQUIRE(phase == Approx(phase_analytic).margin(0.001));
    }
  }

  SECTION("phase and frequency of perturbed system") {
    sam::RK4System<StuartLandauForcedODE> forced_system(1, 2, alpha, eps, nu);
    sam::RK4System<StuartLandauLinearizedODE> linearized_system(1, 4, alpha);
    forced_system.SetPosition(initial);
    forced_system.Integrate(dt, n_trans);
    for (unsigned int i = 0; i < 10; ++i) {
      forced_system.Integrate(dt, 70);
      state_type pos = forced_system.GetPosition();
      state_type deriv = forced_system.GetDerivative();
      state_type y = forced_system.GetPosition();
      y.insert(y.end(), deriv.begin(), deriv.end());
      double phase_analytic = AnalyticPhase(pos, alpha);
      double freq_analytic = AnalyticFrequency(pos, deriv, alpha);
      std::pair<double, double> phase_freq = sam::FindLinearizedPhaseFrequency(
        y, T, linearized_system, 0, 0, CrossingCondition);
      REQUIRE(phase_freq.first == Approx(phase_analytic).margin(0.001));
      REQUIRE(phase_freq.second == Approx(freq_analytic).margin(0.0001));
    }
  }
}
