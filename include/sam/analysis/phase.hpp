// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_ANALYSIS_PHASE_HPP_
#define INCLUDE_SAM_ANALYSIS_PHASE_HPP_

#include <cmath>
#include <utility>
#include <vector>

namespace sam {

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
double PhaseOnLimitCycle(system_type system, double period, double dt,
                         unsigned int n_osc, unsigned int dimension,
                         condition_func&& condition);

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
double FindPhase(const state_type& position, double period,
                 system_type unperturbed_system, unsigned int n_osc,
                 unsigned int dimension, condition_func&& condition,
                 double tolerance = 1e-7,
                 unsigned int max_iter = 100,
                 unsigned int steps_per_period = 1000);


template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
std::pair<double, double> FindLinearizedPhaseFrequency(const state_type& state,
    double period, system_type linearized_system, unsigned int n_osc,
    unsigned int dimension, condition_func&& condition, double tolerance = 1e-7,
    unsigned int max_iter = 100, unsigned int steps_per_period = 1000);


// Implementation

template <typename system_type, typename condition_func, typename state_type>
double PhaseOnLimitCycle(system_type system, double period, double dt,
                         unsigned int n_osc, unsigned int dimension,
                         condition_func&& condition) {
  system.SetTime(0.);
  std::pair<double, state_type> crossing = IntegrateToCrossingConditional(
      system, n_osc, dimension, dt, condition);
  double delta_t = crossing.first;
  return 2*M_PI*(period - delta_t)/period;
}

template<typename system_type, typename state_type>
bool RelaxOntoLimitCycle(system_type& system, double period, double tolerance,
                         unsigned int max_iter, unsigned int steps_per_period,
                         std::vector<size_t> pos_indx = std::vector<size_t>()) {
  double dt = period / static_cast<double>(steps_per_period);
  if (pos_indx.size() == 0) {
    for (size_t i = 0; i < system.GetPosition().size(); ++i) {
      pos_indx.push_back(i);
    }
  }
  state_type pos_previous;
  double error;
  unsigned int iter = 0;
  do {
    error = 0;
    pos_previous = system.GetPosition();
    system.Integrate(dt, steps_per_period);
    state_type position = system.GetPosition();
    for (auto it = pos_indx.begin(); it != pos_indx.end(); ++it) {
      error += std::fabs(pos_previous[*it] - position[*it]);
    }
    ++iter;
  } while (error > tolerance && iter < max_iter);
  return iter != max_iter;
}

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
double FindPhase(const state_type& position, double period,
                 system_type unperturbed_system, unsigned int n_osc,
                 unsigned int dimension, condition_func&& condition,
                 double tolerance, unsigned int max_iter,
                 unsigned int steps_per_period) {
  double dt = period / static_cast<double>(steps_per_period);
  unperturbed_system.SetTime(0.);
  unperturbed_system.SetPosition(position);
  bool relaxed = RelaxOntoLimitCycle<system_type, state_type>(
    unperturbed_system, period, tolerance, max_iter, steps_per_period);
  if (!relaxed) {
    return -1;
  }
  return PhaseOnLimitCycle(unperturbed_system, period, dt, n_osc, dimension,
                           condition);
}

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
std::pair<double, double> FindLinearizedPhaseFrequency(const state_type& state,
    double period, system_type linearized_system, unsigned int n_osc,
    unsigned int dimension, condition_func&& condition, double tolerance,
    unsigned int max_iter, unsigned int steps_per_period) {
  double dt = period / static_cast<double>(steps_per_period);
  linearized_system.SetTime(0.);
  linearized_system.SetPosition(state);
  std::vector<size_t> pos_indx;
  for (size_t i = 0; i < state.size()/2; ++i) {
    pos_indx.push_back(i);
  }
  bool relaxed = RelaxOntoLimitCycle<system_type, state_type>(
    linearized_system, period, tolerance, max_iter, steps_per_period, pos_indx);
  if (!relaxed) {
    return std::make_pair(-1, -1);
  }

  state_type derivative = linearized_system.GetDerivative();
  state_type position = linearized_system.GetPosition();
  double derivative_norm = 0, scalar_product = 0;
  for (size_t i = 0; i < state.size()/2; ++i) {
    derivative_norm += derivative[i]*derivative[i];
    scalar_product += position[state.size()/2 + i]*derivative[i];
  }
  double frequency = 2*M_PI/period*scalar_product/derivative_norm;

  double phase = PhaseOnLimitCycle(linearized_system, period, dt, n_osc,
                                   dimension, condition);
  return std::make_pair(phase, frequency);
}

}  // namespace sam

#endif  // INCLUDE_SAM_ANALYSIS_PHASE_HPP_
