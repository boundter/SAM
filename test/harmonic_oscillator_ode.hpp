// Copyright 2019 Erik Teichmann <kontakt.teichmann@gmail.com>

#ifndef TEST_HARMONIC_OSCILLATOR_ODE_HPP_
#define TEST_HARMONIC_OSCILLATOR_ODE_HPP_

#include <vector>

/*!
 *  The Harmonic oscillator is described by the ODE
 *  \f[ \ddot{x} = - \omega^2 x. \f]
 *  It has the general analytical solution
 *  \f[ x(t) = A*\sin(\omega*t + \varphi). \f]
 */
class HarmonicOscillatorODE {
 public:
  double omega_;

  /*!
   *  @omega params params will be casted to a double pointer, where the first
   *  entry will be used as the frequency.
   */
  explicit HarmonicOscillatorODE(double omega): omega_(omega) {}

  void operator()(const std::vector<double>& x, std::vector<double>& dx,
                  double t) {
    dx[0] = x[1];
    dx[1] = -omega_*omega_*x[0];
  }
};

#endif  // TEST_HARMONIC_OSCILLATOR_ODE_HPP_
