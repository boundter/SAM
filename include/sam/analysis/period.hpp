// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_ANALYSIS_PERIOD_HPP_
#define INCLUDE_SAM_ANALYSIS_PERIOD_HPP_

#include <cmath>
#include <vector>

#include <sam/analysis/henon.hpp>

namespace sam {

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
double CalculatePeriod(system_type& system, unsigned int n_osc,
                       unsigned int dimension, double dt,
                       condition_func&& condition,
                       double target = 0, double precision = 1e-6,
                       unsigned int max_iter = 100) {
  double t_prev;
  double t_current = IntegrateToCrossingConditional(system, n_osc, dimension,
                                                    dt, condition, target).first;
  unsigned int n_iter = 0;
  do {
    t_prev = t_current;
    t_current = IntegrateToCrossingConditional(system, n_osc, dimension,
                                               dt, condition, target).first;
    ++n_iter;
  } while (std::fabs(t_current - t_prev) < precision && n_iter < max_iter);
  return n_iter == max_iter ? -1 : t_current - t_prev;
}

}  // namespace sam

#endif  // INCLUDE_SAM_ANALYSIS_PERIOD_HPP_
