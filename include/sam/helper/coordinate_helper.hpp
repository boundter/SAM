// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_HELPER_COORDINATE_HELPER_HPP_
#define INCLUDE_SAM_HELPER_COORDINATE_HELPER_HPP_

namespace sam {

/*!
 * \brief Transform cartesian coordinates into hyperspherical ones.
 *
 * Calculates the coordinates on a sphere of the same dimension as
 * the phase space. If the dimension is 1, the corrdinates will be wrapped
 * around the unit circle. The first coordinate is the radius and the later
 * ones are the phases. Careful: in 3-d this is not the same as spherical
 * coordinates with polar angle and azimuth!
 *
 * @param begin An iterator pointing to the first element of the coordinates.
 * @param end An iterator pointing to the last element of the coordinates.
 *
 * @returns The spherical coordinates in the same type as the cartesian
 *          coordinates. The first entry is the magnitude, the following
 *          entries the phases in the range [0, 2*pi].
 */
// TODO(boundter): Handle special cases
template<typename state_type>
state_type CartesianToSpherical(typename state_type::const_iterator begin,
                                typename state_type::const_iterator end);

/*!
 * \brief Transform cartesian coordinates into hyperspherical ones.
 *
 * Calculates the coordinates on a sphere of the same dimension as
 * the phase space. If the dimension is 1, the corrdinates will be wrapped
 * around the unit circle. The first coordinate is the radius and the later
 * ones are the phases. Careful: in 3-d this is not the same as spherical
 * coordinates with polar angle and azimuth!
 *
 * @param cartesian The cartesian coordinates.
 *
 * @returns The spherical coordinates in the same type as the cartesian
 *          coordinates. The first entry is the magnitude, the following
 *          entries the phases in the range [0, 2*pi].
 */
template<typename state_type>
state_type CartesianToSpherical(const state_type& cartesian);

// Implementation

template<typename state_type>
state_type CartesianToSpherical(typename state_type::const_iterator begin,
                                typename state_type::const_iterator end) {
  unsigned int dimension = end - begin;
  state_type spherical;
  // sum_squared = [sum(x_i**2), sum(x_i**2) - x_0**2,
  //                sum(x_i**2) - x_0**2 - x_1**2, ...]
  state_type sum_squared(dimension);
  for (typename state_type::const_iterator i = begin; i != end; ++i) {
    sum_squared[0] += (*i)*(*i);
  }
  for (size_t i = 1; i < dimension; ++i) {
    sum_squared[i] = sum_squared[i-1] - (*(begin + i - 1))*(*(begin + i - 1));
  }
  // radius or distance to the origin
  spherical.push_back(sqrt(sum_squared[0]));
  // phases
  for (size_t i = 0; i < dimension - 2; ++i) {
    spherical.push_back(acos((*(begin + i))/sqrt(sum_squared[i])));
  }
  if ((*(end-1)) < 0) {
    spherical.push_back(2*M_PI-acos((*(begin + dimension-2))
                                    /sqrt(sum_squared[dimension-2])));
  } else {
    spherical.push_back(acos((*(begin + dimension-2))
                             /sqrt(sum_squared[dimension-2])));
  }
  return spherical;
}

template<typename state_type>
state_type CartesianToSpherical(const state_type& cartesian) {
  return CartesianToSpherical<state_type>(cartesian.begin(), cartesian.end());
}

}  // namespace sam

#endif  // INCLUDE_SAM_HELPER_COORDINATE_HELPER_HPP_
