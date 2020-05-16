#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.5;
  
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
  V_ = MatrixXd(n_v_, n_v_);
  V_ << std_a_ * std_a_, 0,
       0 , std_yawdd_ * std_yawdd_;

  // set augmented dimension
  n_aug_ = n_x_ + n_v_;

  // Inicial state vector
  x_ << 0,
        0,
        0,
        0,
        0;


std_laspx_ = 0.15;

  // state uncertainty covariance matrix
  P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
        0, std_laspy_*std_laspy_, 0, 0, 0,
        0, 0, std_radr_*std_radr_, 0, 0,
        0, 0, 0, std_radphi_*std_radphi_, 0,
        0, 0, 0, 0, std_radrd_*std_radrd_;

  //Radar init
  n_radar_ = 3;
  V_rad_ = MatrixXd(n_radar_,n_radar_);
  V_rad_ <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;

  //Lidar init
  n_lidar_ = 2;
  V_lid_ = MatrixXd(n_lidar_, n_lidar_);
  V_lid_ <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;


  //UKF Parameters
  // define spreading parameter
  lambda_ = 3 - n_aug_;  

  // create  matrix with predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // create vector for weights and set weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i = 1; i < weights_.size(); ++i)
  {
      weights_(i) = 0.5 / (n_aug_ + lambda_);
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        double r = meas_package.raw_measurements_[0]; 
        double phi = meas_package.raw_measurements_[1]; 
        
        x_(0) = r * cos(phi);
        x_(1) = r * sin(phi);
    }
    else {
        x_(0) = meas_package.raw_measurements_[0];
        x_(1) = meas_package.raw_measurements_[1];
    }

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    return;
  }

  double dt = static_cast<double>(meas_package.timestamp_ - time_us_) / 1e6;
  if (meas_package.sensor_type_ == meas_package.LASER && use_laser_)
  {
    time_us_ = meas_package.timestamp_;
    Prediction(dt);
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == meas_package.RADAR && use_radar_)
  {
    time_us_ = meas_package.timestamp_;
    Prediction(dt);
    UpdateRadar(meas_package);
  }
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
  PredictMeanAndCovariance(Xsig_pred_, n_x_, lambda_, &x_, &P_); 
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  MatrixXd S;
  VectorXd z_pred;
  MatrixXd Zsig;
  PredictLidarMeasurement(Xsig_pred_, &z_pred, &S, &Zsig);
  
  VectorXd z = VectorXd(n_lidar_);
  z << meas_package.raw_measurements_;

  UpdateState(z, S, Zsig, z_pred); 
  
  double NIS = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
  //std::cout << "Lidar NIS: " << NIS << std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  MatrixXd S;
  VectorXd z_pred;
  MatrixXd Zsig;

  PredictRadarMeasurement(Xsig_pred_, &z_pred, &S, &Zsig);
  
  VectorXd z = VectorXd(n_radar_);
  z << meas_package.raw_measurements_;

  UpdateState(z, S, Zsig, z_pred);
  
  double NIS = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
  //std::cout << "Radar NIS: " << NIS << std::endl;
}


void UKF::GenerateSigmaPoints(VectorXd x, MatrixXd P, double lambda, MatrixXd* x_sig_out) {

  int n_x = x.size();

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
    Xsig.col(i + 1) = x + factor * L.col(i);
    Xsig.col(i + n_x + 1) = x - factor * L.col(i);
  }

  *x_sig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd* x_sig_aug_out) {

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
  P_aug(n_x_, n_x_)= V_(0,0);
  P_aug(n_x_+1, n_x_+1)= V_(1,1);

  GenerateSigmaPoints(x_aug, P_aug, lambda_, &x_sig_aug);

  *x_sig_aug_out = x_sig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd x_sig_aug, double delta_t, MatrixXd* x_sig_out) {

  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

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

void UKF::PredictMeanAndCovariance(MatrixXd x_sig_pred, int n_x, double lambda, VectorXd* x_out, MatrixXd* P_out) {

  // create vector for predicted state
  VectorXd x = VectorXd(n_x);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);

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

void UKF::PredictMeasureMeanAndCovariance(MatrixXd x_sig_pred, int n_x, double lambda, VectorXd* x_out, MatrixXd* P_out) {

  // create vector for predicted state
  VectorXd x = VectorXd(n_x);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);

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
    while (x_diff(1)> M_PI) x_diff(1)-=2.*M_PI;
    while (x_diff(1)<-M_PI) x_diff(1)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose();
  }

  // write result
  *x_out = x;
  *P_out = P;
}



void UKF::PredictRadarMeasurement(MatrixXd Xsig_pred, VectorXd* z_out, MatrixXd* S_out, MatrixXd *Zsig_out)  {

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_radar_, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_radar_);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_radar_,n_radar_);

 // transform sigma points into measurement space
  for(int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);
    
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / Zsig(0,i);                  // r_dot
  }
    
  PredictMeasureMeanAndCovariance(Zsig, n_radar_, lambda_, &z_pred, &S); 
  
  // add measurement noise covariance matrix
  S = S + V_rad_;

  // write result
  *Zsig_out = Zsig;
  *z_out = z_pred;
  *S_out = S;
}

void UKF::PredictLidarMeasurement(MatrixXd Xsig_pred, VectorXd* z_out, MatrixXd* S_out, MatrixXd *Zsig_out)  {

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_lidar_, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_lidar_);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_lidar_,n_lidar_);

 // transform sigma points into measurement space
  for(int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    
    Zsig(0,i) = p_x; 
    Zsig(1,i) = p_y; 
  }
    
  PredictMeasureMeanAndCovariance(Zsig, n_lidar_, lambda_, &z_pred, &S); 

  // add measurement noise covariance matrix
  S = S + V_lid_;

  // write result
  *Zsig_out = Zsig;
  *z_out = z_pred;
  *S_out = S;
}


void UKF::UpdateState(VectorXd z, MatrixXd S, MatrixXd Zsig, VectorXd z_pred) { //}, VectorXd* x_out, MatrixXd* P_out) {

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, z.size());

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
    
  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}
