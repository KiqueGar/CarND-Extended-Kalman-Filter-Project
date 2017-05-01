#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;
	if (estimations.size()!=ground_truth.size()){
		//cout << "CalculateRMSE(): Invalid vector sizes!" << endl;
		rmse << 9999,9999,9999,9999;
	}
	//Acumulate Squared residuals
	for (int i=0; i<estimations.size(); i++){
		VectorXd residual = estimations[i] - ground_truth[i];
		residual= residual.array()*residual.array();
		rmse+=residual;
	}
	//Median:
	rmse/=estimations.size();
	//Square Root:
	rmse=rmse.array().pow(0.5);
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj (3,4);
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float denom_ = pow(pow(px, 2) + pow(py, 2),0.5);
	float denom2_ = pow(pow(px, 2) + pow(py, 2),1.5);
	float nume_ = vx*py - vy*px;
	float numen_ =   vy*px-vx*py;

	Hj << px/denom_, py/denom_, 0, 0,
	      -py/pow(denom_, 2), px/pow(denom_,2), 0, 0,
	      py*nume_/denom2_, px*numen_/denom2_, px/denom_, py/denom_;
	return Hj;
}