// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_SYSTEM_RK4_SYSTEM_HPP_
#define INCLUDE_SAM_SYSTEM_RK4_SYSTEM_HPP_

#include <vector>

#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>

#include <sam/system/generic_system.hpp>

namespace sam {

template<typename ODE, typename state_type = std::vector<double>>
class RK4System: public GenericSystem<ODE, state_type> {
 public:
  template<typename... Ts>
  RK4System(unsigned int system_size, unsigned int dimension,
            Ts... parameters)
      : GenericSystem<ODE, state_type>(system_size, dimension, parameters...) {}

  template<typename observer_type = boost::numeric::odeint::null_observer>
  void Integrate(double dt, unsigned int number_steps,
                 observer_type observer
                     = boost::numeric::odeint::null_observer()) {
      this->t_ = boost::numeric::odeint::integrate_n_steps(
          stepper_, *(this->ode_), this->x_, this->t_, dt, number_steps,
          observer);
  }

 private:
  boost::numeric::odeint::runge_kutta4<state_type> stepper_;
};

}  // namespace sam
#endif  // INCLUDE_SAM_SYSTEM_RK4_SYSTEM_HPP_
