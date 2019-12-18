// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_SYSTEM_DETAIL_RK4_SYSTEM_HPP_
#define INCLUDE_SAM_SYSTEM_DETAIL_RK4_SYSTEM_HPP_

#include <sam/system/rk4_system.hpp>

namespace sam {

template<typename ODE, typename state_type>
template<typename... Ts>
RK4System<ODE, state_type>::RK4System(unsigned int system_size,
                                      unsigned int dimension,
                                      Ts... parameters)
    : GenericSystem<ODE, state_type>(system_size, dimension, parameters...) {}

template<typename ODE, typename state_type>
template<typename observer_type>
void RK4System<ODE, state_type>::Integrate(double dt, unsigned int number_steps,
                                           observer_type observer) {
      this->t_ = boost::numeric::odeint::integrate_n_steps(
          stepper_, *(this->ode_), this->x_, this->t_, dt, number_steps,
          observer);
}

}  // namespace sam

#endif  // INCLUDE_SAM_SYSTEM_DETAIL_RK4_SYSTEM_HPP_
