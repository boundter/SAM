// Copyright 2020 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef INCLUDE_SAM_ANALYSIS_SAVITZKY_GOLAY_HPP_
#define INCLUDE_SAM_ANALYSIS_SAVITZKY_GOLAY_HPP_

#include <vector>
#include <Eigen/Dense>
#include <iostream>

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
  // TODO(boundter): Check n_points is uneven and derivatives < n_points
  // TODO(boundter): Optimize sum multiplications
  Eigen::MatrixXd coefficients = CalculateConvolutionCoefficients(n_points,
    derivatives);
  std::vector<std::vector<double>> a(derivatives);
  for (size_t point = n_points/2; point < state.size() - n_points/2; ++point) {
    for (size_t deriv = 0; deriv < derivatives; ++deriv) {
      double sum = 0.;
      for (auto i = -n_points/2; i <= static_cast<int>(n_points)/2; ++i) {
        sum += coefficients(i + n_points/2, deriv)*state[point + i];
      }
      sum *= factorial(deriv)/pow(stepsize, deriv);
      a[deriv].push_back(sum);
    }
  }
  std::cout << coefficients;
  return a;
}

}  // namespace sam

#endif  // INCLUDE_SAM_ANALYSIS_SAVITZKY_GOLAY_HPP_
