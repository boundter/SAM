// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_ANALYSIS_HENON_HPP_
#define INCLUDE_SAM_ANALYSIS_HENON_HPP_

#include <utility>
#include <vector>

namespace sam {

template<typename system_type, typename state_type = std::vector<double>>
std::pair<double, state_type> HenonTrick(const system_type& system,
                                         unsigned int n_osc,
                                         unsigned int dimension, double target);

// Implementation

template<typename system_type, typename state_type>
std::pair<double, state_type> HenonTrick(const system_type& system,
                                         unsigned int n_osc,
                                         unsigned int dimension,
                                         double target) {
  std::pair<unsigned int, unsigned int> dimensionality = system.GetDimension();
  // dimensionality.second is the dimension of each oscillator
  size_t indx = dimensionality.second * n_osc + dimension;
  state_type derivative = system.GetDerivative();
  for (size_t i = 0; i < derivative.size(); ++i) {
    if (i == indx) {
      derivative[indx] = 1./derivative[indx];
      continue;
    } else {
      derivative[i] /= derivative[indx];
    }
  }

  double t = system.GetTime();
  state_type state = system.GetPosition();
  double dx = target - state[indx];
  t += derivative[indx]*dx;
  for (size_t i = 0; i < state.size(); ++i) {
    if (i == indx) {
      state[i] = target;
      continue
    } else {
      state[i] += derivative[i]*dx;
    }
  }
  std::pair<double, state_type> ret_values(t, state);
  return ret_values;
}

}  // namespace sam

#endif  // INCLUDE_SAM_ANALYSIS_HENON_HPP