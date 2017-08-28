#include "Kalman.h"

template <int state_dim, int measurement_dim>
Kalman<state_dim, measurement_dim>::Kalman(Matrix<state_dim, 1> x, Matrix<state_dim, state_dim> P, Matrix<state_dim, state_dim> F,
   Matrix<state_dim, 1> u, Matrix<measurement_dim, 1> z, Matrix<measurement_dim, state_dim> H,
    Matrix<measurement_dim, measurement_dim> R, Matrix<state_dim, state_dim> Q=NULL) {
  _x = x;
  _P = P;
  _F = F;
  _u = u;
  _z = z;
  _H = H;
  _R = R;
  if(Q == NULL) {
    Matrix<state_dim, state_dim> newQ;
    newQ.Fill(0);
    Q = newQ;
  }
  _Q = Q;
}
template <int state_dim, int measurement_dim>
Kalman<state_dim, measurement_dim>::~Kalman(){}

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
  for(uint8_t i = 0; i < _P.getColCount(); i++)
    I(i, i) = 1;
  _x = _x + K * y;
  _P = (I - K * _H) * _P;
}

template <int state_dim, int measurement_dim>
void Kalman<state_dim, measurement_dim>::state_transition() {
  _x = _F * _x + _u;
  _P = _F * _P * ~_F + _Q;
}