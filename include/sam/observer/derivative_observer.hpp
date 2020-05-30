// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OBSERVER_DERIVATIVE_OBSERVER_HPP_
#define INCLUDE_SAM_OBSERVER_DERIVATIVE_OBSERVER_HPP_

#include <vector>

namespace sam {

template<typename system_type, typename state_type>
struct DerivativeObserver {
  system_type& system_;
  std::vector<state_type>& derivative_;
  std::vector<double>& time_;

  explicit DerivativeObserver(system_type& system,
                              std::vector<state_type>& position,
                              std::vector<double>& time);

  void operator()(const state_type& x, double t) const;
};

// Implementation

template<typename system_type, typename state_type>
DerivativeObserver<system_type, state_type>::DerivativeObserver(
    system_type& system, std::vector<state_type>& derivative,
    std::vector<double>& time)
    : system_(system), derivative_(derivative), time_(time) {}

template<typename system_type, typename state_type>
void DerivativeObserver<system_type, state_type>::operator()(
    const state_type& x, double t) const {
  derivative_.push_back(system_.GetDerivative());
  time_.push_back(t);
}

}  // namespace sam

#endif  // INCLUDE_SAM_OBSERVER_DERIVATIVE_OBSERVER_HPP_
