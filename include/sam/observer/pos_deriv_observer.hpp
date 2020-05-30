// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OBSERVER_POS_DERIV_OBSERVER_HPP_
#define INCLUDE_SAM_OBSERVER_POS_DERIV_OBSERVER_HPP_

#include <vector>

#include <sam/observer/position_observer.hpp>
#include <sam/observer/derivative_observer.hpp>

namespace sam {

template<typename system_type, typename state_type>
struct PosDerivObserver {
  explicit PosDerivObserver(system_type& system,
                            std::vector<state_type>& position,
                            std::vector<state_type>& derivative,
                            std::vector<double>& time);

  void operator()(const state_type& x, double t) const;

 protected:
  std::shared_ptr<sam::PositionObserver<state_type>> position_observer_;
  std::shared_ptr<sam::DerivativeObserver<system_type, state_type>>
      derivative_observer_;
  std::vector<double> dummy_time_;
};

// Implementation

template<typename system_type, typename state_type>
PosDerivObserver<system_type, state_type>::PosDerivObserver(system_type& system,
    std::vector<state_type>& position, std::vector<state_type>& derivative,
    std::vector<double>& time) {
  position_observer_ = std::make_shared<sam::PositionObserver<state_type>>(
      position, time);
  derivative_observer_ = std::make_shared<sam::DerivativeObserver<system_type,
                                                                  state_type>>(
      system, derivative, dummy_time_);
}

template<typename system_type, typename state_type>
void PosDerivObserver<system_type, state_type>::operator()(const state_type& x,
                                                           double t) const {
    position_observer_->operator()(x, t);
    derivative_observer_->operator()(x, t);
}

}  // namespace sam

#endif  // INCLUDE_SAM_OBSERVER_POS_DERIV_OBSERVER_HPP_
