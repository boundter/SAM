// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_OBSERVER_BASE_OBSERVER_HPP_
#define INCLUDE_SAM_OBSERVER_BASE_OBSERVER_HPP_

template<typename state_type>
struct BaseObserver {
  virtual void operator()(const state_type& x, double t) {}
};

#endif  // INCLUDE_SAM_OBSERVER_BASE_OBSERVER_HPP_
