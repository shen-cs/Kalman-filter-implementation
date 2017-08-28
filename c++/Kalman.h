/**
 Kalman filter Arduino library implemented with <BasicLinearAlgebra.h> 
 Authored by Shen Chang-Shao
 08/27/2017
 state_dim: length of state vector
 measurement_dim: length of measurement vector
 x: state vector
 P: covariance matrix of state vector
 F: state_transition matrix
 u: motion vector
 z: measurement vector
 H: measurement function (projection matrix from state to measurement)
 R: covariance matrix of measurement vector
 Q: process noise covariance
*/
#ifndef KALMAN_H
#define KALMAN_H
#include <BasicLinearAlgebra.h>
template <int state_dim, int measurement_dim>
class Kalman 
{
public:
  Kalman(Matrix<state_dim, 1> x, Matrix<state_dim, state_dim> P, Matrix<state_dim, state_dim> F,
   Matrix<state_dim, 1> u, Matrix<measurement_dim, 1> z, Matrix<measurement_dim, state_dim> H,
    Matrix<measurement_dim, measurement_dim> R, Matrix<state_dim, state_dim> Q=NULL);
  ~Kalman();
  void filterOnce(Matrix<measurement_dim, 1> measurement);
  Matrix<state_dim, 1> getState();
  Matrix<state_dim, state_dim> getCovariance();
private:
  void measurement_update();
  void state_transition();
  Matrix<state_dim, 1> _x;
  Matrix<state_dim, state_dim> _P;
  Matrix<state_dim, 1> _F;
  Matrix<state_dim, 1> _u;
  Matrix<measurement_dim, 1> _z;
  Matrix<measurement_dim, state_dim> _H;
  Matrix<measurement_dim> _R;
  Matrix<state_dim, state_dim> _Q;
};

#endif