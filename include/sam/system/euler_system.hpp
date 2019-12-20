// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_SYSTEM_EULER_SYSTEM_HPP_
#define INCLUDE_SAM_SYSTEM_EULER_SYSTEM_HPP_

#include <vector>

#include <boost/numeric/odeint/integrate/null_observer.hpp>

#include <sam/system/generic_system.hpp>

namespace sam {

template<typename ODE, typename state_type = std::vector<double>>
class EulerSystem: public GenericSystem<ODE, state_type> {
 public:
  template<typename... Ts>
  explicit EulerSystem(unsigned int system_size, unsigned int dimension,
                       Ts... parameters);

  template<typename observer_type = boost::numeric::odeint::null_observer>
  void Integrate(double dt, unsigned int number_steps,
                 observer_type observer
                     = boost::numeric::odeint::null_observer());

 private:
  template<typename system_type, typename observer_type>
  double EulerMethod(system_type system, state_type& x, double t, double dt,
                     unsigned int number_steps, observer_type observer);
};

// Implementation

template<typename ODE, typename state_type>
template<typename... Ts>
EulerSystem<ODE, state_type>::EulerSystem(unsigned int system_size,
                                          unsigned int dimension,
                                          Ts... parameters)
      : GenericSystem<ODE, state_type>(system_size, dimension, parameters...) {}

template<typename ODE, typename state_type>
template<typename observer_type>
void EulerSystem<ODE, state_type>::Integrate(double dt,
                                             unsigned int number_steps,
                                             observer_type observer) {
    this->t_ = EulerMethod(*(this->ode_), this->x_, this->t_, dt, number_steps,
                           observer);
}

template<typename ODE, typename state_type>
template<typename system_type, typename observer_type>
double EulerSystem<ODE, state_type>::EulerMethod(system_type system,
                                                 state_type& x, double t,
                                                 double dt,
                                                 unsigned int number_steps,
                                                 observer_type observer) {
  observer(x, t);
  for (unsigned int i = 0; i < number_steps; ++i) {
    // TODO(boundter): copy the value? how to best initialize?
    state_type dx = x;
    system(x, dx, t);
    // TODO(boundter): Rather use iterators
    for (size_t j = 0; j < x.size(); ++j) x[j] += dx[j]*dt;
    t += dt;
    observer(x, t);
  }
  return t;
}

}  // namespace sam

#endif  // INCLUDE_SAM_SYSTEM_EULER_SYSTEM_HPP_