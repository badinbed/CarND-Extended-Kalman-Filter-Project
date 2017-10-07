#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include <functional>

#include "Eigen/Dense"

class KalmanFilter {
public:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;


  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   */
  void Init(const Eigen::VectorXd &x_in, const Eigen::MatrixXd &P_in);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param F_in Transition matrix
   * @param Q_in Process covariance matrix
   */
  void Predict(const Eigen::MatrixXd &F_in, const Eigen::MatrixXd &Q_in);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   * @param H Measurement matrix
   * @param R_in Measurement covariance matrix
   */
  void Update(const Eigen::VectorXd &z, const Eigen::MatrixXd &H_in, const Eigen::MatrixXd &R_in);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   * @param h measurement function
   * @param Hj Linearized measurement matrix
   * @param R_in Measurement covariance matrix
   */
  void UpdateEKF(const Eigen::VectorXd &z, std::function<Eigen::VectorXd(const Eigen::VectorXd &, const Eigen::VectorXd &)> h, const Eigen::MatrixXd &Hj, const Eigen::MatrixXd &R_in);

 private:

  /**
   * Helper function that updates the state internaly given H and y
   * @param H Measurement matrix
   * @param y differnze between z and the projected state x into measurement space
   * @param R_in Measurement covariance matrix
   */
  void UpdateInternal(const Eigen::VectorXd &y, const Eigen::MatrixXd & H, const Eigen::MatrixXd &R_in);

};

#endif /* KALMAN_FILTER_H_ */
