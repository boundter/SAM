// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_ANALYSIS_PERIOD_HPP_
#define INCLUDE_SAM_ANALYSIS_PERIOD_HPP_

#include <cmath>
#include <vector>

#include "./henon.hpp"

namespace sam {

struct PeriodParameters : CrossingParameters {
  double period_precision = 1e-6;
  unsigned int period_max_iter = 100;
};

template <typename system_type, typename condition_func,
          typename state_type = std::vector<double>>
double CalculatePeriod(system_type& system, double dt,
                       condition_func&& condition, PeriodParameters params) {
  double t_prev;
  double t_current = IntegrateToCrossingConditional(system, dt, condition,
                                                    params).first;
  unsigned int n_iter = 0;
  do {
    t_prev = t_current;
    t_current = IntegrateToCrossingConditional(system, dt, condition,
                                               params).first;
    ++n_iter;
  } while (std::fabs(t_current - t_prev) < params.period_precision
           && n_iter < params.period_max_iter);
  return n_iter == params.period_max_iter ? -1 : t_current - t_prev;
}

}  // namespace sam

#endif  // INCLUDE_SAM_ANALYSIS_PERIOD_HPP_
