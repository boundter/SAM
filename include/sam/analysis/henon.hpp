// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_ANALYSIS_HENON_HPP_
#define INCLUDE_SAM_ANALYSIS_HENON_HPP_

#include <cmath>
#include <utility>
#include <vector>

namespace sam {

/*!
 * \brief Parameters for the point of crossing of an axis.
 */
struct CrossingParameters {
  unsigned int n_osc = 0;  //! The oscillator which should be considered.
  unsigned int dimension = 0;  //! The dimension to consider.
  double target = 0.;  //! The target value for the dimension.
};

/*!
 * \brief Return crossing of axis with Henon trick.
 *
 * The Henon allows for the calculation of the crossing of an axis by using
 * the relevant dimension as an independent variable. This allows a very precise
 * approximation. The actual approximation is done using the Euler-method, so
 * a close start to the target is needed.
 * For more information see his paper at:
 * https://www.sciencedirect.com/science/article/abs/pii/0167278982900343
 *
 * @param system The system with a state close to the crossing.
 * @param params The parameters for the crossing.
 *
 * @returns A pair of time and state at the time of crossing.
 */
template<typename system_type, typename state_type = std::vector<double>>
std::pair<double, state_type> HenonTrick(const system_type& system,
                                         CrossingParameters params);

/*!
 * \brief Integrate the system to the next regular point after axis crossing.
 *
 * The system is integrated to the next regular point after the crossing of
 * an axis, where the value of the desired crossing can be set. Using the
 * Henon-Trick, the point and time of crossing is returned.
 *
 * @param system The system to integrate to the crossing.
 * @param params The parameters for the crossing.
 * @param dt The integration timestep.
 *
 * @returns A pair of time and state at the time of crossing.
 *
 * TODO: Implement max iter
 * TODO: Implement condition/predicate
 */
template<typename system_type, typename state_type = std::vector<double>>
std::pair<double, state_type> IntegrateToCrossing(system_type& system,
                                                  double dt,
                                                  CrossingParameters params);

/*!
 * \brief Integrate the system to the next regular point after axis crossing.
 *
 * The system is integrated to the next regular point after the crossing of
 * an axis, where the value of the desired crossing can be set additionally some
 * given condition has to be fulfilled. Using the Henon-Trick, the point and
 * time of crossing is returned.
 *
 * @param system The system to integrate to the crossing.
 * @param n_osc The oscillator which schould be considered (starts at 0).
 * @param dimension The dimension to consider (start at 0).
 * @param dt The integration timestep.
 * @param condtion A boolean function that takes the current position.
 * @param target The target value for the dimension.
 *
 * @returns A pair of time and state at the time of crossing.
 *
 * TODO: Implement max iter
 * TODO: Implement condition/predicate
 */
template<typename system_type, typename condition_func,
         typename state_type = std::vector<double>>
std::pair<double, state_type> IntegrateToCrossingConditional(
    system_type& system, double dt, condition_func&& condition,
    CrossingParameters params);

// Implementation

template<typename system_type, typename state_type>
std::pair<double, state_type> HenonTrick(const system_type& system,
                                         CrossingParameters params) {
  std::pair<unsigned int, unsigned int> dimensionality = system.GetDimension();
  // dimensionality.second is the dimension of each oscillator
  size_t indx = dimensionality.second * params.n_osc + params.dimension;
  state_type derivative = system.GetDerivative();
  for (size_t i = 0; i < derivative.size(); ++i) {
    if (i == indx) {
      derivative[indx] = 1./derivative[indx];
      continue;
    } else {
      derivative[i] /= derivative[indx];
    }
  }

  double t = system.GetTime();
  state_type state = system.GetPosition();
  double dx = params.target - state[indx];
  // Euler Integration
  t += derivative[indx]*dx;
  for (size_t i = 0; i < state.size(); ++i) {
    if (i == indx) {
      state[i] = params.target;
      continue;
    } else {
      state[i] += derivative[i]*dx;
    }
  }
  std::pair<double, state_type> ret_values(t, state);
  return ret_values;
}

template<typename system_type, typename state_type = std::vector<double>>
std::pair<double, state_type> IntegrateToCrossing(system_type& system,
                                                  double dt,
                                                  CrossingParameters params) {
  return IntegrateToCrossingConditional(system, dt,
    [](state_type) { return true; }, params);
}

template<typename system_type, typename condition_func,
         typename state_type = std::vector<double>>
std::pair<double, state_type> IntegrateToCrossingConditional(
    system_type& system, double dt, condition_func&& condition,
    CrossingParameters params) {
  std::pair<unsigned int, unsigned int> dimensionality = system.GetDimension();
  // dimensionality.second is the dimension of each oscillator
  size_t indx = dimensionality.second * params.n_osc + params.dimension;
  state_type previous_state;
  do {
    previous_state = system.GetPosition();
    system.Integrate(dt, 1);
  } while (std::copysign(1., previous_state[indx] - params.target)
           == std::copysign(1., system.GetPosition()[indx] - params.target)
           || !condition(system.GetPosition()));
  return HenonTrick(system, params);
}

}  // namespace sam

#endif  // INCLUDE_SAM_ANALYSIS_HENON_HPP_
