// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_ANALYSIS_PHASE_HPP_
#define INCLUDE_SAM_ANALYSIS_PHASE_HPP_

#include <cmath>
#include <utility>
#include <vector>

#include "./henon.hpp"

namespace sam {

struct PhaseParameters : CrossingParameters {
  double phase_tolerance = 1e-7;
  unsigned int phase_max_iter = 100;
  unsigned int phase_steps_per_period = 1000;
};

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
double PhaseOnLimitCycle(system_type system, double period, double dt,
                         condition_func&& condition, CrossingParameters params);

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
double FindPhase(const state_type& position, double period,
                 system_type unperturbed_system, condition_func&& condition,
                 PhaseParameters params);


template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
std::pair<double, double> FindLinearizedPhaseFrequency(const state_type& state,
    double period, system_type linearized_system, condition_func&& condition,
    PhaseParameters params);


// Implementation

template <typename system_type, typename condition_func, typename state_type>
double PhaseOnLimitCycle(system_type system, double period, double dt,
                         condition_func&& condition,
                         CrossingParameters params) {
  system.SetTime(0.);
  std::pair<double, state_type> crossing = IntegrateToCrossingConditional(
      system, dt, condition, params);
  double delta_t = crossing.first;
  return 2*M_PI*(period - delta_t)/period;
}

template<typename system_type, typename state_type>
bool RelaxOntoLimitCycle(system_type& system, double period,
                         PhaseParameters params,
                         std::vector<size_t> pos_indx = std::vector<size_t>()) {
  double dt = period / static_cast<double>(params.phase_steps_per_period);
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
    system.Integrate(dt, params.phase_steps_per_period);
    state_type position = system.GetPosition();
    for (auto it = pos_indx.begin(); it != pos_indx.end(); ++it) {
      error += std::fabs(pos_previous[*it] - position[*it]);
    }
    ++iter;
  } while (error > params.phase_tolerance && iter < params.phase_max_iter);
  return iter != params.phase_max_iter;
}

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
double FindPhase(const state_type& position, double period,
                 system_type unperturbed_system, condition_func&& condition,
                 PhaseParameters params) {
  double dt = period / static_cast<double>(params.phase_steps_per_period);
  unperturbed_system.SetTime(0.);
  unperturbed_system.SetPosition(position);
  bool relaxed = RelaxOntoLimitCycle<system_type, state_type>(
    unperturbed_system, period, params);
  if (!relaxed) {
    return -1;
  }
  return PhaseOnLimitCycle(unperturbed_system, period, dt,
                           condition, params);
}

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
std::pair<double, double> FindLinearizedPhaseFrequency(const state_type& state,
    double period, system_type linearized_system, condition_func&& condition,
    PhaseParameters params) {
  double dt = period / static_cast<double>(params.phase_steps_per_period);
  linearized_system.SetTime(0.);
  linearized_system.SetPosition(state);
  std::vector<size_t> pos_indx;
  for (size_t i = 0; i < state.size()/2; ++i) {
    pos_indx.push_back(i);
  }
  bool relaxed = RelaxOntoLimitCycle<system_type, state_type>(
    linearized_system, period, params, pos_indx);
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

  double phase = PhaseOnLimitCycle(linearized_system, period, dt, condition,
                                   params);
  return std::make_pair(phase, frequency);
}

}  // namespace sam

#endif  // INCLUDE_SAM_ANALYSIS_PHASE_HPP_
