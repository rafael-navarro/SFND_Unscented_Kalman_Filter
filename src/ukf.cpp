#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //state dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  //Process initialize
  is_initialized_ = false;

  n_v_ = 2;
  MatrixXd V_ = MatrixXd(n_v_, n_v_);
  V_ << std_a_ * std_a_, 0,
       0 , std_yawdd_ * std_yawdd_;

  n_aug_ = n_x_ + n_v_;

  //UKF Parameters
  lambda_ = 3 - n_aug_;  

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // create vector for weights and set weights
  VectorXd weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i = 1; i < weights_.size(); ++i)
  {
      weights_(i) = 0.5 / (n_aug_ + lambda_);
  }


  // state covariance matrix P ??? Esta mal
  P_ << std_laspx_, 0, 0, 0, 0,
        0, std_laspy_, 0, 0, 0,
        0, 0, std_radr_, 0, 0,
        0, 0, 0, std_radphi_, 0,
        0, 0, 0, 0, std_radphi_;


  // // measurement covariance
  // kf_.R_ = MatrixXd(2, 2);
  // kf_.R_ << 0.0225, 0,
  //           0, 0.0225;

  // // measurement matrix
  // kf_.H_ = MatrixXd(2, 4);
  // kf_.H_ << 1, 0, 0, 0,
  //           0, 1, 0, 0;

  // // the initial transition matrix F_
  // kf_.F_ = MatrixXd(4, 4);
  // kf_.F_ << 1, 0, 1, 0,
  //           0, 1, 0, 1,
  //           0, 0, 1, 0,
  //           0, 0, 0, 1;


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  MatrixXd Xsig_aug;
  AugmentedSigmaPoints(&Xsig_aug);
  SigmaPointPrediction(Xsig_aug, delta_t, &Xsig_pred_);
  PredictMeanAndCovariance(Xsig_pred_, lambda_, &x_, &P_); 
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
    // VectorXd z_pred = H_ * x_;
    // VectorXd y = z - z_pred;
    // MatrixXd Ht = H_.transpose();
    // MatrixXd S = H_ * P_ * Ht + R_;
    // MatrixXd Si = S.inverse();
    // MatrixXd PHt = P_ * Ht;
    // MatrixXd K = PHt * Si;

    // //new estimate
    // x_ = x_ + (K * y);
    // long x_size = x_.size();
    // MatrixXd I = MatrixXd::Identity(x_size, x_size);
    // P_ = (I - K * H_) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}




void UKF::GenerateSigmaPoints(VectorXd x, MatrixXd P, double lambda, MatrixXd* x_sig_out) {

  int n_x = x.size();

  // define spreading parameter
  //double lambda = 3 - n_x;

  // create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

  // calculate square root of P
  MatrixXd L = P.llt().matrixL();

  //State vector as first col
  Xsig.col(0) = x;
  
  // set sigma points as columns of matrix Xsig
  double factor = sqrt(lambda + n_x);
  for(int i = 0; i < n_x; ++i)
  {
    Xsig.col(i + 1) = x_ + factor * L.col(i);
    Xsig.col(i + n_x + 1) = x - factor * L.col(i);
  }

  *x_sig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd* x_sig_aug_out) {

  // set augmented dimension
  //int n_aug = n_x_ + n_aug_;

  // // Process noise standard deviation longitudinal acceleration in m/s^2
  // double std_a = 0.2;

  // // Process noise standard deviation yaw acceleration in rad/s^2
  // double std_yawdd = 0.2;

  // define spreading parameter
  //double lambda = 3 - n_aug;

  // // set example state
  // VectorXd x = VectorXd(n_x);
  // x <<   5.7441,
  //        1.3800,
  //        2.2049,
  //        0.5015,
  //        0.3528;

  // create example covariance matrix
  // MatrixXd P = MatrixXd(n_x, n_x);
  // P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
  //         -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
  //          0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
  //         -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
  //         -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  // create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd x_sig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
  x_aug.head(n_x_) = x_; 

  // create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_x_, n_x_) = V_;

  GenerateSigmaPoints(x_aug, P_aug, lambda_, &x_sig_aug);

  *x_sig_aug_out = x_sig_aug;
  
  // // create square root matrix
  // MatrixXd L = P_aug.llt().matrixL();

  // // create augmented sigma points
    
  // //State vector as first col
  // Xsig_aug.col(0) = x_aug;
  
  // // set sigma points as columns of matrix Xsig
  // double factor = sqrt(lambda_ + n_aug_);
  // for(int i = 0; i < n_aug_; ++i)
  // {
  //   Xsig_aug.col(i + 1) = x_aug + factor * L.col(i);
  //   Xsig_aug.col(i + n_aug_ + 1) = x_aug - factor * L.col(i);
  // }
  
  // /**
  //  * Student part end
  //  */

  // // print result
  // std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  // // write result
  // *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd x_sig_aug, double delta_t, MatrixXd* x_sig_out) {

  // // set state dimension
  // int n_x = 5;

  // // set augmented dimension
  // int n_aug = 7;

  // // create example sigma point matrix
  // MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  // Xsig_aug <<
  //   5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
  //     1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
  //   2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
  //   0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
  //   0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
  //        0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
  //        0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //double delta_t = 0.1; // time diff in sec

  // predict sigma points
  for(int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double p_x = x_sig_aug(0, i);
    double p_y = x_sig_aug(1, i);
    double v = x_sig_aug(2, i);
    double yaw = x_sig_aug(3, i);
    double yawd = x_sig_aug(4, i);
    double nu_a = x_sig_aug(5, i);
    double nu_yawdd = x_sig_aug(6, i);
    
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    
    // write predicted sigma point into right column
    Xsig_pred(0, i) = px_p;
    Xsig_pred(1, i) = py_p;
    Xsig_pred(2, i) = v_p;
    Xsig_pred(3, i) = yaw_p;
    Xsig_pred(4, i) = yawd_p;
  }

  // write result
  *x_sig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(MatrixXd x_sig_pred, double lambda, VectorXd* x_out, MatrixXd* P_out) {

  // // set state dimension
  // int n_x = 5;

  // // set augmented dimension
  // int n_aug = 7;

  // // define spreading parameter
  // double lambda = 3 - n_aug;

  // // create example matrix with predicted sigma points
  // MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  // Xsig_pred <<
  //        5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  //          1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  //         2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  //        0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  //         0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  // create vector for weights
  // VectorXd weights = VectorXd(2*n_aug_+1);
  
  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // // set weights
  // weights(0) = lambda / (lambda + n_aug_);
  // for(int i = 1; i < 2 * n_aug_ + 1; ++i)
  // {
  //     weights(i) = 0.5/(n_aug_+lambda);
  // }

  // predict state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
  {  
    x = x + weights_(i) * x_sig_pred.col(i);
  }
  // predict state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) { 
    // state difference
    VectorXd x_diff = x_sig_pred.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose();
  }

  // write result
  *x_out = x;
  *P_out = P;
}

