#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <math.h>


#define EIGEN_NO_DEBUG
#define PI  3.14159265358979323846

using namespace std;
using namespace Eigen;
//using Eigen::MatrixXd;
//using Eigen::VectorXd;

// ~~~~~~~~~~~~~~~~~~~~~~~~
// Declaration of functions 
// ~~~~~~~~~~~~~~~~~~~~~~~~

vector<MatrixXd> read(string filename);
VectorXd latlong(VectorXd& coords);
MatrixXd designMat(MatrixXd& satPos, VectorXd& x0);
VectorXd misclosure(VectorXd& p, MatrixXd& satPos, VectorXd& x0);
VectorXd misclosure(VectorXd& p, MatrixXd& satPos, VectorXd& x0, VectorXd& base);
double elev(double x, double y, double z, VectorXd& latlon);
MatrixXd rotation(double lat, double lon);
MatrixXd weight(MatrixXd& temp);

//#endif
