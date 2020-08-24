#include "functions.h"

// reads the data from a file to a vector of matrices
vector<MatrixXd> read(string filename) {

	// open file
	ifstream fin(filename);
	if (!fin) {
		cout << "Error opening file. Please try again." << endl;
		exit(1);
	}

	// epoch
	double time;

	// num of satellites 
	int sats;

	// vector of matrices to store satellite data [PRN, x, y, z, pseudorange]
	vector <MatrixXd> sat;

	// place holder matrix
	MatrixXd mat;

	while (!fin.eof())
	{
		// read values
		if (fin >> time) {
			fin >> sats;
			mat = MatrixXd::Zero(sats + 1, 5);

			mat(0, 0) = time;

			for (int i = 1; i <= sats; i++) {
				for (int j = 0; j < 5; j++)
					fin >> mat(i, j);
			}
			sat.push_back(mat);
		}
		else
			continue;
	}

	fin.close();
	return sat;
};

// calculates latitude and longitude from given coordinates 
VectorXd latlong(VectorXd& coords) {

	VectorXd latlon(2);

	//WGS84 ellipsoid parameters
	double a = 6378137.0;
	double b = 6356752.3;
	double e = (pow(a, 2) - pow(b, 2)) / (pow(a, 2));

	latlon(1) = atan2(coords(1, 0), coords(0, 0));

	double p = sqrt(pow(coords(0, 0), 2) + pow(coords(1, 0), 2));
	double h0 = 0;
	double h = 1;
	double lat = atan((coords(2, 0) / p) * (1 / (1 - e)));
	double W, N;
	while ((h - h0) > 0.001)
	{
		h0 = h;
		W = sqrt(1 - e * pow(sin(lat), 2));
		N = a / W;
		h = (p / cos(lat)) - N;
		lat = atan((coords(2, 0) / p) * (1 / (1 - e * (N / (N + h)))));
	}
	latlon(0) = lat;
	return latlon;
}

// builds design matrix A
MatrixXd designMat(MatrixXd& sat, VectorXd& x0) {

	MatrixXd A(sat.rows() - 1, x0.size());
	double p;

	for (int i = 0; i < A.rows(); i++) {

		p = sqrt(pow((x0(0) - sat(i + 1, 1)), 2) + pow((x0(1) - sat(i + 1, 2)), 2) + pow((x0(2) - sat(i + 1, 3)), 2));

		for (int j = 0; j < A.cols() - 1; j++)
			A(i, j) = (x0(j) - sat(i + 1, j + 1)) / p;

		A(i, 3) = -1;
	}

	return A;
}

// builds the misclosure vector w
VectorXd misclosure(VectorXd& p, MatrixXd& sat, VectorXd& x0)
{
	VectorXd w(sat.rows() - 1);
	for (int i = 0; i < w.rows(); i++) {
		w(i) = sqrt(pow((x0(0) - sat(i + 1, 1)), 2) + pow((x0(1) -
			sat(i + 1, 2)), 2) + pow((x0(2) - sat(i + 1, 3)), 2)) - p(i) - x0(3);
	}

	return w;
};

// overloaded function that builds the misclosure vector w for single differencing
VectorXd misclosure(VectorXd& p, MatrixXd& sat, VectorXd& x0, VectorXd& base)
{
	VectorXd w(sat.rows() - 1);

	for (int i = 0; i < w.rows(); i++) { 
		w(i) = (sqrt(pow((x0(0) - sat(i + 1, 1)), 2) + pow((x0(1) - sat(i + 1, 2)), 2)
			+ pow((x0(2) - sat(i + 1, 3)), 2)) - (sqrt(pow((base(0) - sat(i + 1, 6)), 2)
				+ pow((base(1) - sat(i + 1, 7)), 2) + pow((base(2) - sat(i + 1, 8)), 2))) - p(i) - x0(3));
	}
	return w;
};

// calculates the elevation angle
double elev(double x, double y, double z, VectorXd& latlon) {
	VectorXd xyz(3);
	xyz << x, y, z;
	MatrixXd R = rotation(latlon(0), latlon(1));
	VectorXd enu = R * xyz;
	double elev = asin(enu(2) / sqrt(pow(enu(0), 2) + pow(enu(1), 2) + pow(enu(2), 2)));
	return elev;
};

// calculates the rotation matrix for xyz to enu 
MatrixXd rotation(double lat, double lon) {

	// Populate transformation matrix, R
	MatrixXd R;
	R.resize(3, 3);

	R << -sin(lon), cos(lon), 0,
		-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat),
		cos(lat)* cos(lon), cos(lat)* sin(lon), sin(lat);

	return R;
};

// calculates the weight matrix P as a function of elevation angle
MatrixXd weight(MatrixXd& temp) {
	MatrixXd P = MatrixXd::Zero(temp.rows() - 1, temp.rows() - 1);
	for (int i = 1; i < temp.rows(); i++)
		P(i - 1, i - 1) = pow((1 / pow(sin(temp(i, 5)), 2)) + (1 / pow(sin(temp(i, 10)), 2)), -1);

	return P;
};

