#include "kalman_filter.h"

#include <iostream>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const VectorXd &x_in, const MatrixXd &P_in) {
  x_ = x_in;
  P_ = P_in;
}

void KalmanFilter::Predict(const Eigen::MatrixXd &F_in, const Eigen::MatrixXd &Q_in) {
  x_ = F_in * x_;
  MatrixXd Ft = F_in.transpose();
  P_ = F_in * P_ * Ft + Q_in;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &H_in, const MatrixXd &R_in) {
  VectorXd y = z - H_in * x_;
  UpdateInternal(y, H_in, R_in);
}

void KalmanFilter::UpdateEKF(const VectorXd &z, std::function<VectorXd(const VectorXd &, const VectorXd &)> h, const MatrixXd &Hj_in, const MatrixXd &R_in) {

    UpdateInternal(h(z, x_), Hj_in, R_in);
}

void KalmanFilter::UpdateInternal(const VectorXd &y, const MatrixXd & H, const MatrixXd &R_in) {
    MatrixXd Ht = H.transpose();
    MatrixXd S = H * P_ * Ht + R_in;
    MatrixXd Si = S.inverse();
    MatrixXd K = (P_ * Ht )* Si;

    x_ = x_ + (K * y);
    P_ = (MatrixXd::Identity(x_.size(), x_.size()) - K * H) * P_;
}
