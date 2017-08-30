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
#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H
#include <BasicLinearAlgebra.h>
template <int state_dim, int measurement_dim>
class Kalman 
{
public:
  Kalman() {}
  Kalman(Matrix<state_dim, 1> x, Matrix<state_dim, state_dim> P, Matrix<state_dim, state_dim> F,
   Matrix<state_dim, 1> u, Matrix<measurement_dim, 1> z, Matrix<measurement_dim, state_dim> H,
    Matrix<measurement_dim, measurement_dim> R, Matrix<state_dim, state_dim> Q) {
    _x = x;
  _P = P;
  _F = F;
  _u = u;
  _z = z;
  _H = H;
  _R = R;
  _Q = Q;
  }
  ~Kalman(){}
  void filterOnce(Matrix<measurement_dim, 1> measurement);
  void initialize(Matrix<state_dim, 1> x, Matrix<state_dim, state_dim> P, Matrix<state_dim, state_dim> F,
   Matrix<state_dim, 1> u, Matrix<measurement_dim, 1> z, Matrix<measurement_dim, state_dim> H,
    Matrix<measurement_dim, measurement_dim> R, Matrix<state_dim, state_dim> Q);
  Matrix<state_dim, 1> getState();
  Matrix<state_dim, state_dim> getCovariance();
private:
  void measurement_update();
  void state_transition();
  Matrix<state_dim, 1> _x;
  Matrix<state_dim, state_dim> _P;
  Matrix<state_dim, state_dim> _F;
  Matrix<state_dim, 1> _u;
  Matrix<measurement_dim, 1> _z;
  Matrix<measurement_dim, state_dim> _H;
  Matrix<measurement_dim> _R;
  Matrix<state_dim, state_dim> _Q;
};

template<int state_dim, int measurement_dim>
void Kalman<state_dim, measurement_dim>::initialize(Matrix<state_dim, 1> x, Matrix<state_dim, state_dim> P, Matrix<state_dim, state_dim> F,
   Matrix<state_dim, 1> u, Matrix<measurement_dim, 1> z, Matrix<measurement_dim, state_dim> H,
    Matrix<measurement_dim, measurement_dim> R, Matrix<state_dim, state_dim> Q) {
  _x = x;
  _P = P;
  _F = F;
  _u = u;
  _z = z;
  _H = H;
  _R = R;
  _Q = Q;
}

template <int state_dim, int measurement_dim>
void Kalman<state_dim, measurement_dim>::filterOnce(Matrix<measurement_dim, 1> measurement) {
  _z = measurement;
  state_transition();
  measurement_update();
}

template <int state_dim, int measurement_dim>
Matrix<state_dim, 1> Kalman<state_dim, measurement_dim>::getState() {
  return _x;
}

template <int state_dim, int measurement_dim>
Matrix<state_dim, state_dim> Kalman<state_dim, measurement_dim>::getCovariance() {
  return _P;
}

template <int state_dim, int measurement_dim>
void Kalman<state_dim, measurement_dim>::measurement_update() {
  Matrix<measurement_dim, 1> y = _z - _H * _x;
  Matrix<measurement_dim, measurement_dim> S = _H * _P * ~_H + _R;
  Matrix<state_dim, measurement_dim> K = _P * ~_H * S.Inverse();
  Matrix<state_dim, state_dim> I;
  I.Fill(0);
  for(uint8_t i = 0; i < _P.GetColCount(); i++)
    I(i, i) = 1;
  _x = _x + K * y;
  _P = (I - K * _H) * _P;
}

template <int state_dim, int measurement_dim>
void Kalman<state_dim, measurement_dim>::state_transition() {
  _x = _F * _x + _u;
  _P = _F * _P * ~_F + _Q;
}
#endif