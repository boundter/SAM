// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_SYSTEM_GENERIC_SYSTEM_HPP_
#define INCLUDE_SAM_SYSTEM_GENERIC_SYSTEM_HPP_

#include <stdexcept>  // length_error,
#include <utility>
#include <vector>

#include "../helper/coordinate_helper.hpp"

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

  GenericSystem(const GenericSystem<ODE, state_type>& other_system);

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
   *  \brief Set the current time for the system.
   */
  void SetTime(double t);

  /*!
   *  \brief Return the derivative at the current position at time.
   */
  state_type GetDerivative() const;

  /*!
   *  \brief Return the dimensionality of the system.
   *
   *  Get the dimensionality as a pair in the form
   *  (number of oscillators, dimensionality of the oscillator).
   */
  std::pair<unsigned int, unsigned int> GetDimension() const;

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

  /*!
   *  \brief Returns the average position of all elements in the state
   *  space.
   *
   * The mean field is the average position in the state space. It has the same
   * dimension as one unit of the system.
   *
   * @returns A position in state space that is the average of all units.
   */
  state_type CalculateMeanField() const;


  /*!
   * \brief Calculates the coordinates on a sphere of the same dimension as
   * the phase space. If the deimension is 1, the corrdinates will be wrapped
   * around the unit circle. The first coordinate is the radius and the later
   * ones are the phases. Careful: in 3-d this is not the same as spherical
   * coordinates with polar angle and azimuth!
   */
  state_type CalculateMeanFieldSpherical() const;

 protected:
  std::unique_ptr<ODE> ode_;
  unsigned int N_, d_;
  state_type x_;
  double t_;

  /*!
   *  \brief Initializer for the GenericSystem class. No ODE will be
   *  initialized, it is intended for the use in inherited classes.
   */
  explicit GenericSystem(unsigned int system_size, unsigned int dimension);


  // Calculate the mean field using iterators to allow easy calculation of the
  // mean field in a network
  state_type CalculateMeanField(
    const typename state_type::const_iterator start,
    const typename state_type::const_iterator end) const;


  // Calculate the mean field using iterators to allow easy calculation of the
  // mean field in a network
  state_type CalculateMeanFieldSpherical(
    const typename state_type::const_iterator start,
    const typename state_type::const_iterator end) const;
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
GenericSystem<ODE, state_type>::GenericSystem(
    const GenericSystem<ODE, state_type>& other_system) {
  N_ = other_system.N_;
  d_ = other_system.d_;
  x_ = other_system.x_;
  t_ = other_system.t_;
  ode_ = std::unique_ptr<ODE>(new ODE(*other_system.ode_));
}


template<typename ODE, typename state_type>
GenericSystem<ODE, state_type>::GenericSystem(unsigned int system_size,
                                              unsigned int dimension) {
  N_ = system_size;
  d_ = dimension;
  x_.resize(N_*d_);
  t_ = 0.;
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
void GenericSystem<ODE, state_type>::SetTime(double t) {
  t_ = t;
}

template<typename ODE, typename state_type>
state_type GenericSystem<ODE, state_type>::GetDerivative() const {
  state_type intermediate(x_.size());
  ode_->operator()(x_, intermediate, t_);
  return intermediate;
}

template<typename ODE, typename state_type>
std::pair<unsigned int, unsigned int> GenericSystem<ODE, state_type>::
    GetDimension() const {
  std::pair<unsigned int, unsigned int> dimension;
  dimension.first = N_;
  dimension.second = d_;
  return dimension;
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

template<typename ODE, typename state_type>
state_type GenericSystem<ODE, state_type>::CalculateMeanField() const {
  return CalculateMeanField(x_.begin(), x_.end());
}

template<typename ODE, typename state_type>
state_type GenericSystem<ODE, state_type>::CalculateMeanField(
    const typename state_type::const_iterator start,
    const typename state_type::const_iterator end) const {
  double N = static_cast<double>(end-start)/static_cast<double>(d_);
  if (N != static_cast<unsigned int>(N)) {
    throw std::length_error("Mean Field cannot be calculated, if not all "
                            "oscillators are given.");
  }
  state_type mean_field(d_);
  // TODO(boundter): Check size
  for (typename state_type::const_iterator i = start; i < end; i += d_) {
    for (typename state_type::const_iterator j = i; j < i + d_; ++j) {
      mean_field[j-i] += (*j);
    }
  }
  for (typename state_type::iterator i = mean_field.begin();
       i != mean_field.end(); ++i) {
    (*i) /= static_cast<double>(N);
  }
  return mean_field;
}

template<typename ODE, typename state_type>
state_type GenericSystem<ODE, state_type>::CalculateMeanFieldSpherical() const {
  return CalculateMeanFieldSpherical(x_.begin(), x_.end());
}

template<typename ODE, typename state_type>
state_type GenericSystem<ODE, state_type>::CalculateMeanFieldSpherical(
    const typename state_type::const_iterator start,
    const typename state_type::const_iterator end) const {
  double N = static_cast<double>(end-start)/static_cast<double>(d_);
  if (N != static_cast<unsigned int>(N)) {
    throw std::length_error("Mean Field cannot be calculated, if not all "
                            "oscillators are given.");
  }
  // for d = 1 wrap around unit circle
  state_type spherical_mean_field;
  if (d_ == 1) {
    double x = 0, y = 0;
    for (typename state_type::const_iterator i = start; i != end; ++i) {
      x += cos((*i));
      y += sin((*i));
    }
    spherical_mean_field.resize(2);
    spherical_mean_field[0] = 1./N*sqrt(x*x + y*y);
    spherical_mean_field[1] = atan2(y, x);
  } else {
    spherical_mean_field = CartesianToSpherical(CalculateMeanField(start,
                                                                    end));
  }
  return spherical_mean_field;
}

}  // namespace sam

#endif  // INCLUDE_SAM_SYSTEM_GENERIC_SYSTEM_HPP_
