// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_ANALYSIS_SAVITZKY_GOLAY_HPP_
#define INCLUDE_SAM_ANALYSIS_SAVITZKY_GOLAY_HPP_

#include <Eigen/Dense>

#include <stdexcept>
#include <vector>

namespace sam {

template<typename state_type = std::vector<double>>
std::vector<std::vector<double>> SavitzkyGolayFilter(double stepsize,
                                                     state_type state,
                                                     unsigned int n_points,
                                                     unsigned int derivatives);


// Implementation

Eigen::MatrixXd CalculateConvolutionCoefficients(unsigned int n_points,
                                                 unsigned int derivatives) {
  unsigned int n_left = n_points/2;
  Eigen::MatrixXd J(n_points, derivatives);
  for (unsigned int i = 0; i < n_points; ++i) {
    double z = static_cast<int>(i) - static_cast<int>(n_left);
    for (unsigned int j = 0; j < derivatives; ++j) {
      J(i, j) = pow(z, j);
    }
  }
  return (J.transpose()*J).inverse()*J.transpose();
}

int factorial(unsigned int n) {
  int factorial = 1.;
  for (unsigned int i = 1; i < n+1; ++i) factorial *= i;
  return factorial;
}

template<typename state_type>
std::vector<std::vector<double>> SavitzkyGolayFilter(double stepsize,
                                                     state_type state,
                                                     unsigned int n_points,
                                                     unsigned int derivatives) {
  if (n_points%2 == 0) {
    throw std::invalid_argument("Number of points for Savitzky-Golay Filter "
        "has to be uneven.");
  }
  if (n_points < derivatives) {
    throw std::invalid_argument("Number of points has to be bigger than "
        "number of derivatives in Savitzky-Golay Filter.");
  }
  Eigen::MatrixXd coefficients = CalculateConvolutionCoefficients(n_points,
      derivatives);
  std::vector<std::vector<double>> a(derivatives);
  for (size_t point = n_points/2; point < state.size() - n_points/2; ++point) {
    for (size_t deriv = 0; deriv < derivatives; ++deriv) {
      double sum = 0.;
      for (int i = -static_cast<int>(n_points)/2;
           i <= static_cast<int>(n_points)/2; ++i) {
        sum += coefficients(deriv, i + n_points/2)*state[point + i];
      }
      sum *= factorial(deriv)/pow(stepsize, deriv);
      a[deriv].push_back(sum);
    }
  }
  return a;
}

}  // namespace sam

#endif  // INCLUDE_SAM_ANALYSIS_SAVITZKY_GOLAY_HPP_
