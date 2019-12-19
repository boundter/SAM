// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_SYSTEM_GENERIC_SYSTEM_HPP_
#define INCLUDE_SAM_SYSTEM_GENERIC_SYSTEM_HPP_

#include <stdexcept>  // length_error,
#include <vector>

#include <sam/helper/coordinate_helper.hpp>

namespace sam {

/*!
 * A wrapper for ODEs. It is especially useful for integrating coupled
 * differential equations. Be careful: Integration does not support
 * polymorphism.
 */
template<typename ODE, typename state_type = std::vector<double>>
class GenericSystem {
 public:
  /*!
   *  Initialzer for the GenericSystem class.
   *  @param system_size Number of elements in the system.
   *  @param dimension Number of ODEs per element.
   *  @param parameters All the parameters that need to be passed to the ODE.
   */
  template<typename ...Ts>
  explicit GenericSystem(unsigned int system_size, unsigned int dimension,
                         Ts... parameters);

  /*!
   *  \brief Return the position in the state space for all elements.
   */
  state_type GetPosition() const;

  /*!
   *  \brief Set the position in the state space.
   *
   *  Set the new position in the state space. If the new position has the
   *  wrong size an exception will be thrown.
   *
   *  @throws std::length_error
   */
  void SetPosition(const state_type& new_position);

  /*!
   *  \brief Get the time of the system.
   */
  double GetTime() const;

  /*!
   *  \brief Return the derivative at the current position at time.
   */
  state_type GetDerivative() const;

  /*!
   * \brief Integrate the system whith an observer.
   *
   * Integrate the system whith an observer. The observer is a user-specified
   * struct/class that receives the current time and state. If none is
   * specified the null_observer will be used which does nothing.
   *
   * @param dt Timestep for the integration.
   * @param number_steps The total number of timesteps.
   * @param observer The observer of the integration. The observer will be
   *  called after every timestep with the current position and time.
   */
  template<typename observer_type>
  void Integrate(double dt, unsigned int number_steps,
                 observer_type observer) {}

  /*!
   * \brief Change the number of oscillators.
   *
   * Change the number of oscillators in the system. No new ODE object will be
   * generated, so if ti depends on the number of oscillattors the parameters
   * will have to be changed afterwards.
   *
   * @param N The new number of oscillators.
   */
  void Resize(unsigned int N);

  /*!
   *  \brief Set the parameters for the ODE. This creates a new pointer
   *  to the ODE.
   */
  template<typename ...Ts>
  void SetParameters(Ts... parameters);

  /*!
   * \brief Transform cartesian coordinates into hyperspherical ones.
   *
   * Calculates the coordinates on a sphere of the same dimension as
   * the phase space. If the dimension is 1, the corrdinates will be wrapped
   * around the unit circle. The first coordinate is the radius and the later
   * ones are the phases. Careful: in 3-d this is not the same as spherical
   * coordinates with polar angle and azimuth!
   */
  state_type GetPositionSpherical() const;

//   /*!
//     *  \brief Returns the average position of all elements in the state
//     *  space.
//     */
//   state_type CalculateMeanField() {
//     return CalculateMeanField(x_.begin(), x_.end());
//   }


//   /*!
//     * \brief Calculates the coordinates on a sphere of the same dimension as
//     * the phase space. If the deimension is 1, the corrdinates will be wrapped
//     * around the unit circle. The first coordinate is the radius and the later
//     * ones are the phases. Careful: in 3-d this is not the same as spherical
//     * coordinates with polar angle and azimuth!
//     */
//   state_type CalculateMeanFieldSpherical() {
//     return CalculateMeanFieldSpherical(x_.begin(), x_.end());
//   }




//   // TODO(boundter): Check for NULL-Pointer
//   // TODO(boundter): remove infinte loop
//   /*!
//     *  \brief Calculate the period of the average of all elements in the
//     *  state space.
//     *
//     *  This function calculates the period of the mean field. For one
//     *  oscillator this is the same as using the position. The period is
//     *  calculated by measuring the time between crossings of the Poincare
//     *  manifold, a surface in the state space that is crossed transversally
//     *  by the trajectory.
//     *
//     *  @param CrossedPoincareManifold This is a user-defined funtion. It
//     *  recieves the previous state and the current state and should
//     *  return True, if the mean field crossed the Poincare manifold which
//     *  defines the surface where the period is measured. Otherwise it
//     *  should return False.
//     *
//     *  @param n_crossings Number of crossings of the Poincare manifold
//     *  between periods. This is useful for oscillators with higher
//     *  winding numbers.
//     *
//     *  @param ApproximateCrossingPoincareManifold This is a user-defined
//     *  function. It receives the previous and current time and state from
//     *  this it should approximate the time of crossing. If none is
//     *  specified it will approximate the time by taking the middle between
//     *  the time of the previous and current step.
//     */
//   template <typename observer_type = boost::numeric::odeint::null_observer>
//   double CalculatePeriod(unsigned int n_average, double dt,
//                          bool (*CrossedPoincareManifold)(
//                              const state_type& /*previous_state*/,
//                              const state_type& /*current_state*/),
//                          unsigned int n_crossings = 1,
//                          double (*ApproximateCrossingPoincareManifold)(
//                              const state_type& /*previous_state*/,
//                              double /*previous_t*/,
//                              const state_type& /*current_state*/,
//                              double /*current_t*/) = NULL,
//                          observer_type observer
//                              = boost::numeric::odeint::null_observer()) {
//     std::vector<double> times_of_crossing;
//     state_type previous_state = CalculateMeanField();
//     unsigned int n_observed_crossings = 0;
//     observer(x_, t_);
//     // we need one more time of crossing than periods
//     while (times_of_crossing.size() < n_average + 1) {
//       Integrate(dt, 1);
//       observer(x_, t_);
//       state_type current_state = CalculateMeanField();
//       if (CrossedPoincareManifold(previous_state, current_state)) {
//         n_observed_crossings += 1;
//         if (n_observed_crossings == n_crossings) {
//           double current_time = GetTime();
//           double t_approx;
//           if (ApproximateCrossingPoincareManifold == NULL) {
//             t_approx = BifurcationZerothOrderCrossingPoincare(
//                 previous_state, current_time - dt, current_state,
//                 current_time);
//           } else {
//             t_approx = ApproximateCrossingPoincareManifold(
//                 previous_state, current_time - dt, current_state,
//                 current_time);
//           }
//           times_of_crossing.push_back(t_approx);
//           n_observed_crossings = 0;
//         }
//       }
//       previous_state = current_state;
//     }
//     return CalculatePeriodFromCrossings(times_of_crossing);
//   }


//   /*!
//     *  \brief Integrates the system by a phase in the range (0, 2*pi)
//     *
//     *  Integrates the system by a phase \f$ \varphi \f$ defined with the period
//     *  \f$ T \f$ as
//     *  \f[ \varphi = 2\pi \frac{t}{T}. \f]
//     *
//     *  @param phase length of the phase to integrate by.
//     *  @param T the period of the system.
//     *  @param dt the timestep.
//     */
//   template <typename observer_type = boost::numeric::odeint::null_observer>
//   void IntegrateByPhase(double phase, double T, double dt,
//                         observer_type observer
//                             = boost::numeric::odeint::null_observer()) {
//     double time_difference = phase/(2*M_PI)*T;
//     double steps = time_difference/dt;
//     // casting steps to unsigned int will loose the last fraction of a
//     // timestep
//     this->Integrate(dt, static_cast<unsigned int>(steps), observer);
//     // interate the last fraction
//     double remaining_time = time_difference
//         - dt*static_cast<double>(static_cast<unsigned int>(steps));
//     this->Integrate(remaining_time, 1);
//     // call the observer separately, because it will be called at the begin of
//     // the Integrate function and at the end, so we have one call to many
//     observer(x_, t_);
//   }


//   /*!
//     *  \brief Integrates the system to a reference point
//     *
//     *  Integrates the system to a reference point/surface of the mean field, defined by
//     *  the CrossedReference function. The position is checked after every
//     *  integration, so the state after calling this function will always be
//     *  slightly after the reference point/surface. In the worst case scenario
//     *  the difference will be the timestep dt.
//     *
//     *  @param dt the timestep
//     *  @param CrossedReference this is a user defined function taking the
//     *  previous mean field and the current mean field and returns a boolean
//     *  value of the crossing of the reference point/surface.
//     */
//   template <typename observer_type = boost::numeric::odeint::null_observer>
//   void IntegrateToReference(double dt,
//                             bool (*CrossedReference)(
//                                 const state_type& /*previous_state*/,
//                                 const state_type& /*current_state*/),
//                             observer_type observer
//                                 = boost::numeric::odeint::null_observer()) {
//     observer(x_, t_);
//     state_type previous_state = this->CalculateMeanField();
//     this->Integrate(dt, 1);
//     observer(x_, t_);
//     state_type current_state = this->CalculateMeanField();
//     while (!CrossedReference(previous_state, current_state)) {
//       this->Integrate(dt, 1);
//       observer(x_, t_);
//       previous_state = current_state;
//       current_state = this->CalculateMeanField();
//     }
//   }


 protected:
  std::unique_ptr<ODE> ode_;
  unsigned int N_, d_;
  state_type x_;
  double t_;
  typedef typename state_type::iterator iterator_type;


//   /*!
//     *  \brief Initializer for the GenericSystem class. No ODE will be
//     *  initialized, it is intended for the use in inherited classes.
//     */
//   GenericSystem(unsigned int system_size, unsigned int dimension) {
//     N_ = system_size;
//     d_ = dimension;
//     x_.resize(N_*d_);
//     t_ = 0.;
//   }


//   // Calculate the mean field using iterators to allow easy calculation of the
//   // mean field in a network
//   state_type CalculateMeanField(const iterator_type start,
//                                 const iterator_type end) {
//     double N = static_cast<double>(end-start)/static_cast<double>(d_);
//     if (N != static_cast<unsigned int>(N)) {
//       throw std::length_error("Mean Field cannot be calculated, if not all "
//                               "oscillators are given.");
//     }
//     state_type mean_field(d_);
//     // TODO(boundter): Check size
//     for (iterator_type i = start; i < end; i += d_) {
//       for (iterator_type j = i; j < i + d_; ++j) {
//         mean_field[j-i] += (*j);
//       }
//     }
//     for (iterator_type i = mean_field.begin(); i != mean_field.end(); ++i) {
//       (*i) /= static_cast<double>(N);
//     }
//     return mean_field;
//   }


//   // Calculate the mean field using iterators to allow easy calculation of the
//   // mean field in a network
//   state_type CalculateMeanFieldSpherical(const iterator_type start,
//                                          const iterator_type end) {
//     double N = static_cast<double>(end-start)/static_cast<double>(d_);
//     if (N != static_cast<unsigned int>(N)) {
//       throw std::length_error("Mean Field cannot be calculated, if not all "
//                               "oscillators are given.");
//     }
//     // for d = 1 wrap around unit circle
//     state_type spherical_mean_field;
//     if (d_ == 1) {
//       double x = 0, y = 0;
//       for (iterator_type i = start; i != end; ++i) {
//         x += cos((*i));
//         y += sin((*i));
//       }
//       spherical_mean_field.resize(2);
//       spherical_mean_field[0] = 1./N*sqrt(x*x + y*y);
//       spherical_mean_field[1] = atan2(y, x);
//     } else {
//       spherical_mean_field = CartesianToSpherical(CalculateMeanField(start,
//                                                                      end));
//     }
//     return spherical_mean_field;
//   }


//  private:
//   boost::numeric::odeint::runge_kutta4<state_type> stepper_;

//   // Return the time of crossing as t_before_crossing + dt/2.
//   static double BifurcationZerothOrderCrossingPoincare(
//       const state_type& previous_state, double previous_t,
//       const state_type& current_state, double current_t) {
//     return (previous_t + current_t)/2.;
//   }


//   // Calculates the differences between the elements in a vector and
//   // averages this difference vector. This is used to calculate the
//   // period from measuring the crossing of a Poincare manifold.
//   double CalculatePeriodFromCrossings(
//       const std::vector<double>& times_of_crossing) {
//     double period = 0.;
//     for (size_t i = 1; i < times_of_crossing.size(); ++i) {
//       period += times_of_crossing[i] - times_of_crossing[i-1];
//     }
//     return period/static_cast<double>(times_of_crossing.size()-1);
//   }
};

// Implementation

template<typename ODE, typename state_type>
template<typename... Ts>
GenericSystem<ODE, state_type>::GenericSystem(unsigned int system_size,
                                              unsigned int dimension,
                                              Ts... parameters) {
  N_ = system_size;
  d_ = dimension;
  x_.resize(N_*d_);
  t_ = 0.;
  ode_ = std::make_unique<ODE>(parameters...);
}

template<typename ODE, typename state_type>
state_type GenericSystem<ODE, state_type>::GetPosition() const {
  return x_;
}

template<typename ODE, typename state_type>
void GenericSystem<ODE, state_type>::SetPosition(
      const state_type& new_position) {
  if (new_position.size() != N_*d_) {
    throw std::length_error("Trying to set new position of wrong length!");
  }
  x_ = new_position;
}

template<typename ODE, typename state_type>
double GenericSystem<ODE, state_type>::GetTime() const {
  return t_;
}

template<typename ODE, typename state_type>
state_type GenericSystem<ODE, state_type>::GetDerivative() const {
  state_type intermediate(x_.size());
  ode_->operator()(x_, intermediate, t_);
  return intermediate;
}

template<typename ODE, typename state_type>
void GenericSystem<ODE, state_type>::Resize(unsigned int N) {
  N_ = N;
  x_.resize(N_*d_);
}

template<typename ODE, typename state_type>
template<typename... Ts>
void GenericSystem<ODE, state_type>::SetParameters(Ts... parameters) {
  ode_ = std::make_unique<ODE>(parameters...);
}

template<typename ODE, typename state_type>
state_type GenericSystem<ODE, state_type>::GetPositionSpherical() const {
  if (d_ == 1) {
    return x_;
  } else {
    state_type spherical;
    for (unsigned int i = 0; i < N_; ++i) {
      state_type coord = CartesianToSpherical<state_type>(
          x_.begin() + i*d_, x_.begin() + (i+1)*d_);
      spherical.insert(spherical.end(), coord.begin(), coord.end());
    }
    return spherical;
  }
}

}  // namespace sam

#endif  // INCLUDE_SAM_SYSTEM_GENERIC_SYSTEM_HPP_
