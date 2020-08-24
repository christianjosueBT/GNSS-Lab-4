#include "functions.h"

int main()
{
	// opening outut files to store our results
	ofstream out1("Part1.txt");
	ofstream out2("Part2.txt");
	ofstream out3("Part3.txt"); 
	
	 if (out1.fail() || out2.fail() || out3.fail()) {
		cout << "Could not open output file";
		exit(1);
	}

	// declare variables to be used
	MatrixXd A;
	VectorXd w;
	MatrixXd P;
	VectorXd delta;
	VectorXd xhat;
	MatrixXd N;
	VectorXd u;
	VectorXd p;
	MatrixXd temp;
	double threshold = 0.000001;

	// given
	VectorXd rover_xyz(3);
	rover_xyz << -1641891.178, -3664882.194, 4939966.860;
	VectorXd base_xyz(3);
	base_xyz << -1641945.119, -3664803.806, 4940009.246;

	// conversion XYZ to (lat, lon)
	VectorXd rover_latlong = latlong(rover_xyz);
	VectorXd base_latlong = latlong(base_xyz);

	// initial rover values
	VectorXd x0(4);
	x0 << rover_xyz(0), rover_xyz(1), rover_xyz(2), 0;

	// PARRT 1

	// read rover data in
	vector<MatrixXd> rover = read("rover.data");
	
	// for all epochs in the data file...
	for (int i = 0; i < rover.size(); i++) {

		temp = rover[i];

		// pseudorange measurements
		p = VectorXd::Zero(temp.rows() - 1);
		for (int j = 0; j < temp.rows() - 1; j++)
			p(j) = temp(j + 1, 4);

		// resize
		A.resize(p.size(), 4);
		w.resize(p.size());
		P = MatrixXd::Identity(p.size(), p.size());

		// build
		A = designMat(temp, x0);
		w = misclosure(p, temp, x0);
		N = A.transpose() * P * A;
		u = A.transpose() * P * w;
		delta = -N.inverse() * u;
		xhat = x0 + delta;

		double max = abs(delta(0));
		int iterations = 0;

		// iterations
		while (max >= threshold) {

			VectorXd x0_LOOP = xhat;
			A = designMat(temp, x0_LOOP);
			w = misclosure(p, temp, x0_LOOP);
			N = (A.transpose()) * P * A;
			u = A.transpose() * P * w;
			delta = -(N.inverse()) * u;
			xhat = x0_LOOP + delta;

			max = abs(delta(0));
			for (int i = 0; i < delta.size(); i++) {
				if (max < abs(delta(i)))
					max = abs(delta(i));
			}
			iterations++;
		}

		// outputting the results to a text file for later use
		out1 << fixed << setprecision(6) << xhat(0) << " " << xhat(1) << " " << xhat(2) << " " << xhat(3) << endl;
	}

	// PARRT 2

	// reading in base data and declaring variables to be used
	vector<MatrixXd> base = read("Base.data");
	vector<MatrixXd> epochs;
	int count;


	// for all epochs in the rover data file
	// matching epochs of rover and base 
	for (int i = 0; i < rover.size(); i++) {

		int idx_base = 0;
		bool f_base = false;

		// for all epochs in base data file
		for (int j = 0; j < base.size(); j++) {
			// if epochs match
			if (abs(base[j](0, 0) - rover[i](0, 0)) < 0.01)
			{
				f_base = true;
				idx_base = j;
				temp = MatrixXd::Zero(1, 11);
				temp(0, 0) = base[idx_base](0, 0);
				break;
			}
		}

		if (f_base == false)
			continue;

		count = 1;
		// for every sat in each epoch of rover
		for (int k = 1; k < rover[i].rows(); k++) {

			// for every of sats in each epoch of base
			for (int l = 1; l < base[idx_base].rows(); l++) {

				// storing information of rover and base into a matrix to be pushed back 
				if (rover[i](k, 0) == base[idx_base](l, 0)) {

					temp.conservativeResize(count + 1, 11);
					temp(count, 0) = rover[i](k, 0);
					temp(count, 1) = rover[i](k, 1);
					temp(count, 2) = rover[i](k, 2);
					temp(count, 3) = rover[i](k, 3);
					temp(count, 4) = rover[i](k, 4); 
					temp(count, 5) = elev(rover[i](k, 1), rover[i](k, 2),
						rover[i](k, 3), rover_latlong);
					temp(count, 6) = base[idx_base](l, 1);
					temp(count, 7) = base[idx_base](l, 2);
					temp(count, 8) = base[idx_base](l, 3);
					temp(count, 9) = base[idx_base](l, 4); 
					temp(count, 10) = elev(base[idx_base](l, 1), base[idx_base](l, 2),
						base[idx_base](l, 3), base_latlong);
					count += 1;
				}
			}
		}
		epochs.push_back(temp);
	}

	 // for all matching epochs ... 
	for (int i = 0; i < epochs.size(); i++) {

		temp = epochs[i];

		A = designMat(temp, x0);

		p = VectorXd::Zero(temp.rows() - 1);
		for (int j = 1; j < temp.rows(); j++)
			p(j - 1) = temp(j, 4) - temp(j, 9);

		w = misclosure(p, temp, x0, base_xyz);
		P = MatrixXd::Identity(p.size(), p.size());
		N = (A.transpose()) * P * A;
		u = A.transpose() * P * w;
		delta = -(N.inverse()) * u;
		xhat = x0 + delta;

		double max = abs(delta(0));
		while (max >= threshold) {

			VectorXd x0_LOOP = xhat;
			A = designMat(temp, x0_LOOP);
			w = misclosure(p, temp, x0_LOOP, base_xyz);
			N = (A.transpose()) * P * A;
			u = A.transpose() * P * w;
			delta = -(N.inverse()) * u;
			xhat = x0_LOOP + delta;

			max = abs(delta(0));
			for (int i = 0; i < delta.rows(); i++) {
				if (max < abs(delta(i)))
					max = abs(delta(i));
			}
		}

		// outputting the results to a text file for later use
		out2 << fixed << setprecision(6) << xhat(0) << " " << xhat(1) << " " << xhat(2) << " " << xhat(3) << endl;
	}

	// PARRT 3

	// for all matching epochs
	for (int i = 0; i < epochs.size(); i++) {

		temp = epochs[i];

		x0 << base_xyz(0),
			base_xyz(1),
			base_xyz(2),
			0;

		A = designMat(temp, x0);

		p = VectorXd::Zero(temp.rows() - 1);
		for (int j = 1; j < temp.rows(); j++)
			p(j - 1) = temp(j, 4) - temp(j, 9);

		w = misclosure(p, temp, x0, base_xyz);
		P = MatrixXd::Identity(p.size(), p.size());
		P = weight(temp);
		N = (A.transpose()) * P * A;
		u = A.transpose() * P * w;
		delta = -(N.inverse()) * u;
		xhat = x0 + delta;

		double max = abs(delta(0));
		while (max >= threshold) {

			VectorXd x0_LOOP = xhat;
			A = designMat(temp, x0_LOOP);
			w = misclosure(p, temp, x0_LOOP, base_xyz);
			N = (A.transpose()) * P * A;
			u = A.transpose() * P * w;
			delta = -(N.inverse()) * u;
			xhat = x0_LOOP + delta;

			max = abs(delta(0));
			for (int i = 0; i < delta.size(); i++) {
				if (max < abs(delta(i)))
					max = abs(delta(i));
			}
		}

		// outputting the results to a text file for later use
		out3 << fixed << setprecision(6) << xhat(0) << " " << xhat(1) << " " << xhat(2) << " " << xhat(3) << endl;
	}

	out1.close();
	out2.close();
	out3.close();


	return 0;
}