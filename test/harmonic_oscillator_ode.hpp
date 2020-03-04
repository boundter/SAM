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
    for (unsigned int i = 0; i < x.size()/2.; ++i) {
      dx[2*i] = x[2*i+1];
      dx[2*i+1] = -omega_*omega_*x[2*i];
    }
  }
};

class CoupledHarmonicOscillatorODE {
 public:
  double omega_1_;
  double omega_2_;
  double coupling_;

  explicit CoupledHarmonicOscillatorODE(double omega_1, double omega_2,
                                        double coupling)
      : omega_1_(omega_1), omega_2_(omega_2), coupling_(coupling) {}

  explicit CoupledHarmonicOscillatorODE(std::vector<double> omega,
                                        double coupling)
      : omega_1_(omega[0]), omega_2_(omega[1]), coupling_(coupling) {}

  void operator()(const std::vector<double>& x, std::vector<double>& dx,
                  double t) {
    dx[0] = x[1];
    dx[1] = -omega_1_*omega_1_*x[0] + coupling_*(x[3] - x[1]);
    dx[2] = x[3];
    dx[3] = -omega_2_*omega_2_*x[2] + coupling_*(x[1] - x[3]);
  }
};

#endif  // TEST_HARMONIC_OSCILLATOR_ODE_HPP_
