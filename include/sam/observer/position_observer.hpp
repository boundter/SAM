// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OBSERVER_POSITION_OBSERVER_HPP_
#define INCLUDE_SAM_OBSERVER_POSITION_OBSERVER_HPP_

#include <vector>

namespace sam {

/*!
 * \brief Observer for the position and the time.
 *
 * The observer saves the position and time of the system during an integration.
 */
template<typename state_type>
struct PositionObserver {
  std::vector<state_type>& position_;
  std::vector<double>& time_;

  explicit PositionObserver(std::vector<state_type>& position,
                            std::vector<double>& time);

  void operator()(const state_type& x, double t) const;
};

// Implementation

template<typename state_type>
PositionObserver<state_type>::PositionObserver(
    std::vector<state_type>& position, std::vector<double>& time)
    : position_(position), time_(time) {}

template<typename state_type>
void PositionObserver<state_type>::operator()(const state_type& x, double t)
    const {
  position_.push_back(x);
  time_.push_back(t);
}

}  // namespace sam

#endif  // INCLUDE_SAM_OBSERVER_POSITION_OBSERVER_HPP_
