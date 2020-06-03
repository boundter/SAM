// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_SYSTEM_RK4_NETWORK_HPP_
#define INCLUDE_SAM_SYSTEM_RK4_NETWORK_HPP_

#include <vector>

#include <boost/numeric/odeint/integrate/check_adapter.hpp>
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>

#include "./generic_network.hpp"

namespace sam {

/*! \brief A network that is integrated with a 4th order Runge-Kutta method.
 *
 * A network that is integrated with a Runge-Kutta method of 4th order. The
 * integrator is implemented by the odeint library in boost
 * (https://www.boost.org/doc/libs/1_72_0/libs/numeric/odeint/doc/html/index.html).
 */
template<typename ODE, typename data_type = double>
class RK4Network: public GenericNetwork<ODE, data_type> {
 public:
  using typename GenericNetwork<ODE, data_type>::node_size_type;
  using typename GenericNetwork<ODE, data_type>::state_type;
  using typename GenericNetwork<ODE, data_type>::matrix_type;

  template<typename... Ts>
  explicit RK4Network(node_size_type node_sizes, unsigned int dimension,
                     Ts... parameters);

  template<typename observer_type = boost::numeric::odeint::null_observer>
  void Integrate(double dt, unsigned int number_steps,
                 observer_type observer
                     = boost::numeric::odeint::null_observer());

 private:
  boost::numeric::odeint::runge_kutta4<state_type> stepper_;
};

// Implementation

template<typename ODE, typename data_type>
template<typename... Ts>
RK4Network<ODE, data_type>::RK4Network(node_size_type node_sizes,
                                      unsigned int dimension,
                                      Ts... parameters)
    : GenericNetwork<ODE, data_type>(node_sizes, dimension, parameters...) {}

template<typename ODE, typename data_type>
template<typename observer_type>
void RK4Network<ODE, data_type>::Integrate(double dt, unsigned int number_steps,
                                           observer_type observer) {
      this->t_ = boost::numeric::odeint::integrate_n_steps(
          stepper_, *(this->ode_), this->x_, this->t_, dt, number_steps,
          observer);
}

}  // namespace sam

#endif  // INCLUDE_SAM_SYSTEM_RK4_NETWORK_HPP_
