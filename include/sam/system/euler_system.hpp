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

}  // namespace sam

#include <sam/system/detail/euler_system.hpp>

#endif  // INCLUDE_SAM_SYSTEM_EULER_SYSTEM_HPP_
