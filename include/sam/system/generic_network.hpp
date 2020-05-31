// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_SYSTEM_GENERIC_NETWORK_HPP_
#define INCLUDE_SAM_SYSTEM_GENERIC_NETWORK_HPP_

// TODO(boundter): Possible without dependence on vector?
// TODO(boundter): Increase test coverage
// TODO(boundter): Derivative and Spherical in matrix form

#include <utility>
#include <vector>

#include <sam/system/generic_system.hpp>

namespace sam {

/*!
 *  This class is a wrapper for the GenericSystem for the use with a network.
 *  A network consists of linked nodes, where each node can consists of
 *  multiple oscillators. There is no limitation as to the type of oscillator,
 *  it can be a phase oscillator or a general limit cycle one. For now it is
 *  limited to a use of oscillators of the same dimensionality.
 */
template<typename ODE, typename data_type = double>
class GenericNetwork:
    protected GenericSystem<ODE, std::vector<data_type>> {
 public:
  typedef std::vector<data_type> state_type;
  typedef std::vector<unsigned int> node_size_type;
  typedef std::vector<state_type> matrix_type;

  /*!
   * The network is initialized to a zero state.
   *
   * @param node_sizes a vector containing the size of every single node
   * @param dimension the dimensionality of the oscillators
   * @param parameters pointer to parameters for the ODE
   */
  template<typename ...Ts>
  explicit GenericNetwork(node_size_type node_sizes, unsigned int dimension,
                 Ts... parameters);

  GenericNetwork(const GenericNetwork<ODE, data_type>& other_network);

  /*!
   * Sets the state of the system using a flattened representation of the form
   * state = {node_1x_1, node_1x_2, ..., node_2x_1, ....}.
   */
  void SetPosition(const state_type& new_state);

  /*!
   * Gets the state in a flattened representation of the from
   * state = {node_1x_1, node_1x_2, ..., node_2x_1, ....}.
   */
  state_type GetPosition() const;

  /*!
   * Gets the indices of the beginning of every new node + (the last index + 1)
   * of the flattened representation.
   */
  node_size_type GetNodeIndices() const;

  /*!
   *  Gets the state as a vector of vector representation
   *  state = {{node_1x_1, node_1x_2, ....}, {node_2x_1, ...}, ...}.
   */
  matrix_type GetNodes() const;

  /*!
   *  Gets the state as a vector of vector representation
   *  state = {{node_1x_1, node_1x_2, ....}, {node_2x_1, ...}, ...},
   * where the coordinates are spherical like in GetPositionSpherical().
   */
  matrix_type GetNodesSpherical() const;

  /*!
   *  Gets the derivative in a flattened representation.
   */
  state_type GetDerivative() const;

  /*!
   *  Gets the derivative in a node representation.
   */
  matrix_type GetDerivativeNodes() const;

  /*!
   * \brief Return the position in the state space in phases for all elements.
   *
   * Calculates the coordinates on a sphere of the same dimension as
   * the phase space. If the dimension is 1, the corrdinates will be wrapped
   * around the unit circle as phases, otherise the first coordinate of every
   * element is the radius and the later ones are the phases.
   * Careful: in 3-d this is not the same as spherical coordinates with polar
   * angle and azimuth!
   */
  state_type GetPositionSpherical() const;

  /*!
   * Gets the time of the system.
   */
  double GetTime() const;

  /*!
   *  \brief Set the current time for the system.
   */
  void SetTime(double t);

  /*!
   *  \brief Set the parameters for the ODE. This creates a new pointer
   *  to the ODE.
   */
  template<typename ...Ts>
  void SetParameters(Ts... parameters);

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
   * @param node_sizes The new number of oscillators in each node.
   */
  void Resize(node_size_type node_sizes);

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

 protected:
  node_size_type node_indices_;
  node_size_type node_sizes_;

  node_size_type CalculateNodeIndices(const node_size_type& node_sizes) const;

  unsigned int CalculateNumberOscillators(const node_size_type& node_sizes)
      const;
};

// Implementation

template<typename ODE, typename data_type>
template<typename ...Ts>
GenericNetwork<ODE, data_type>::GenericNetwork(
    node_size_type node_sizes, unsigned int dimension,
    Ts... parameters)
    : GenericSystem<ODE, state_type>(0, dimension, parameters...) {
  node_indices_ = CalculateNodeIndices(node_sizes);
  node_sizes_ = node_sizes;
  Resize(node_sizes_);
}

template<typename ODE, typename data_type>
GenericNetwork<ODE, data_type>::GenericNetwork(
    const GenericNetwork<ODE, data_type>& other_network)
    : GenericSystem<ODE, state_type>(other_network) {
  node_indices_ = other_network.node_indices_;
  node_sizes_ = other_network.node_sizes_;
}

template<typename ODE, typename data_type>
std::vector<unsigned int> GenericNetwork<ODE, data_type>::CalculateNodeIndices(
    const node_size_type& node_sizes) const {
  node_size_type node_indices = {0};
  unsigned int offset = 0;
  for (size_t i = 0; i < node_sizes.size(); ++i) {
    node_indices.push_back(offset + node_sizes[i]*this->d_);
    offset += node_sizes[i]*this->d_;
  }
  return node_indices;
}

template<typename ODE, typename data_type>
unsigned int GenericNetwork<ODE, data_type>::CalculateNumberOscillators(
    const node_size_type& node_sizes) const {
  unsigned int N = 0;
  for (size_t i = 0; i < node_sizes.size(); ++i) {
    N += node_sizes[i];
  }
  return N;
}

template<typename ODE, typename data_type>
void GenericNetwork<ODE, data_type>::SetPosition(const state_type& new_state) {
  GenericSystem<ODE, state_type>::SetPosition(new_state);
}

template<typename ODE, typename data_type>
std::vector<data_type> GenericNetwork<ODE, data_type>::GetPosition() const {
  return GenericSystem<ODE, state_type>::GetPosition();
}

template<typename ODE, typename data_type>
std::vector<unsigned int> GenericNetwork<ODE, data_type>::GetNodeIndices()
    const {
  return node_indices_;
}

template<typename ODE, typename data_type>
std::vector<std::vector<data_type>> GenericNetwork<ODE, data_type>::GetNodes()
    const {
  matrix_type nodes;
  for (size_t i = 0; i < node_indices_.size() - 1; ++i) {
    nodes.push_back(state_type());
    for (unsigned int j = node_indices_[i]; j < node_indices_[i+1]; ++j) {
      nodes.back().push_back(this->x_[j]);
    }
  }
  return nodes;
}

template<typename ODE, typename data_type>
std::vector<data_type> GenericNetwork<ODE, data_type>::GetDerivative() const {
  return GenericSystem<ODE, state_type>::GetDerivative();
}

template<typename ODE, typename data_type>
std::vector<std::vector<data_type>> GenericNetwork<ODE, data_type>::
    GetDerivativeNodes() const {
  state_type deriv = GenericSystem<ODE, state_type>::GetDerivative();
  matrix_type nodes;
  for (size_t i = 0; i < node_indices_.size() - 1; ++i) {
    nodes.push_back(state_type(deriv.begin() + node_indices_[i],
                               deriv.begin() + node_indices_[i+1]));
  }
  return nodes;
}

template<typename ODE, typename data_type>
std::vector<data_type> GenericNetwork<ODE, data_type>::GetPositionSpherical()
    const {
  return GenericSystem<ODE, state_type>::GetPositionSpherical();
}

template<typename ODE, typename data_type>
double GenericNetwork<ODE, data_type>::GetTime() const {
  return GenericSystem<ODE, state_type>::GetTime();
}

template<typename ODE, typename data_type>
void GenericNetwork<ODE, data_type>::SetTime(double t) {
  return GenericSystem<ODE, state_type>::SetTime(t);
}

template<typename ODE, typename data_type>
template<typename... Ts>
void GenericNetwork<ODE, data_type>::SetParameters(Ts... parameters) {
  GenericSystem<ODE, state_type>::SetParameters(parameters...);
}

template<typename ODE, typename data_type>
std::pair<unsigned int, unsigned int> GenericNetwork<ODE, data_type>::
    GetDimension() const {
  return GenericSystem<ODE, state_type>::GetDimension();
}

template<typename ODE, typename data_type>
void GenericNetwork<ODE, data_type>::Resize(node_size_type node_sizes) {
  node_sizes_ = node_sizes;
  node_indices_ = CalculateNodeIndices(node_sizes_);
  GenericSystem<ODE, state_type>::
      Resize(CalculateNumberOscillators(node_sizes_));
}

template<typename ODE, typename data_type>
std::vector<std::vector<data_type>> GenericNetwork<ODE, data_type>::
    GetNodesSpherical() const {
  matrix_type nodes;
  for (size_t i = 0; i < node_indices_.size() - 1; ++i) {
    nodes.push_back(state_type());
    for (size_t j = node_indices_[i]; j < node_indices_[i+1]; j += this->d_) {
      if (this->d_ == 1) {
        nodes.back().push_back(this->x_[j]);
      } else {
        state_type coord = CartesianToSpherical<state_type>(
          this->x_.begin() + j, this->x_.begin() + j + this->d_);
          nodes.back().insert(nodes.back().end(), coord.begin(), coord.end());
      }
    }
  }
  return nodes;
}

template<typename ODE, typename data_type>
std::vector<data_type> GenericNetwork<ODE, data_type>::CalculateMeanField()
    const {
  return GenericSystem<ODE, state_type>::CalculateMeanField();
}

}  // namespace sam

#endif  // INCLUDE_SAM_SYSTEM_GENERIC_NETWORK_HPP_
