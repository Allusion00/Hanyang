#include <iostream>
#include <vector>

using namespace std;


// check start time
double start_time = clock();

// mesh number
int NC = 80;
int NR = 20;

// Diffusion coefficient [cm]
double D1C = 1.3314f;
double D2C = 0.2042f;
double D1R = 1.2291f;
double D2R = 0.1588f;

// absorption X section [cm-1]
double Xa1C = 0.006431f;
double Xa2C = 0.09992f;
double Xa1R = 0.0005903f;
double Xa2R = 0.01882f;

// fission X section [cm-1]
double Xf1C = 0.003258f;
double Xf2C = 0.1387f;
double Xf1R = 0.00f;
double Xf2R = 0.00f;

// scattering X section [cm-1]
double Xs12C = 0.05979f;
double Xs12R = 0.06163f;

// region size [cm]
double sizeC = 40.0f;
double sizeR = 10.0f;

// mesh size [cm]
double hC = sizeC / NC;
double hR = sizeR / NR;

// initial effective k
double keff_old = 1.0f;
double keff_new = 0.0f;


// flux vector
vector<double> phi1_old(NC + NR + 2, 0.0f); // group 1 flux
vector<double> phi1_new(NC + NR + 2, 0.0f);
vector<double> phi2_old(NC + NR + 2, 0.0f); // group 2 flux
vector<double> phi2_new(NC + NR + 2, 0.0f);

// Boundary condition
/*	vacuum condition : r = 0.5
	current zero (reflective) condition : r = 0; 10E-10
	flux zero condition : r = inf; 10E+10
string BL = "reflective"; // Left boundary
string BR = "flux zero"; // Right boundary */

// Boundary coefficient, r
double rL = 10E-10f;
double rR = 10E+10f;

// Convergence criteria
double crik = 1.0E-5f;
double crip = 1.0E-6f;


int main() {

	// Diffusion coefficient vector define
	vector<double> D1(NC + NR + 1, 0.0f); // group 1
	vector<double> D2(NC + NR + 1, 0.0f); // group 2
	// mesh size vector define
	vector<double> h(NC + NR + 1, 0.0f);
	// absorption X section vector define
	vector<double> Xa1(NC + NR + 1, 0.0f); // group 1
	vector<double> Xa2(NC + NR + 1, 0.0f); // group 2
	// fission X section vector define
	vector<double> Xf1(NC + NR + 1, 0.0f); // group 1
	vector<double> Xf2(NC + NR + 1, 0.0f); // group 2
	// scattering X section vector define
	vector<double> Xs12(NC + NR + 1, 0.0f);
	// removal X section vector define
	vector<double> Xr1(NC + NR + 1, 0.0f); // group 1
	vector<double> Xr2(NC + NR + 1, 0.0f); // group 2

	for (int i = 1; i <= NC; ++i) {
		D1[i] = D1C;
		D2[i] = D2C;
		h[i] = hC;
		Xa1[i] = Xa1C;
		Xa2[i] = Xa2C;
		Xf1[i] = Xf1C;
		Xf2[i] = Xf2C;
		Xs12[i] = Xs12C;
		Xr1[i] = Xa1[i] + Xs12[i];
		Xr2[i] = Xa2[i];
	}
	for (int i = NC + 1; i <= NC + NR; ++i) {
		D1[i] = D1R;
		D2[i] = D2R;
		h[i] = hR;
		Xa1[i] = Xa1R;
		Xa2[i] = Xa2R;
		Xf1[i] = Xf1R;
		Xf2[i] = Xf2R;
		Xs12[i] = Xs12R;
		Xr1[i] = Xa1[i] + Xs12[i];
		Xr2[i] = Xa2[i];
	}

	// beta vector define
	vector<double> beta1(NC + NR + 2, 0.0f);
	vector<double> beta2(NC + NR + 2, 0.0f);
	beta1[0] = rL / 2;
	beta2[0] = rL / 2;
	for (int i = 1; i <= NC + NR; ++i) {
		beta1[i] = D1[i] / h[i];
		beta2[i] = D2[i] / h[i];
	}
	beta1[NC + NR + 1] = rR / 2;
	beta2[NC + NR + 1] = rR / 2;

	// D tilder vector define
	vector<double> D_tilder1(NC + NR + 1, 0.0f);
	vector<double> D_tilder2(NC + NR + 1, 0.0f);
	for (int i = 0; i <= NC + NR; ++i) {
		D_tilder1[i] = 2 * beta1[i] * beta1[i + 1] / (beta1[i] + beta1[i + 1]);
		D_tilder2[i] = 2 * beta2[i] * beta2[i + 1] / (beta2[i] + beta2[i + 1]);
	}
	
	/////////////////////////////////////// Matrix M ///////////////////////////////////////////
	// group 1
	// diagonal vector define
	vector<double> diagonal1(NC + NR + 1, 0.0f);
	for (int i = 1; i <= NC + NR; ++i)
		diagonal1[i] = D_tilder1[i - 1] + D_tilder1[i] + Xr1[i] * h[i];

	// bidiagonal vector define
	vector<double> bidiagonal1(NC + NR + 1, 0.0f);
	for (int i = 1; i <= NC + NR; ++i)
		bidiagonal1[i] = -D_tilder1[i];

	// group 2
	// diagonal vector define
	vector<double> diagonal2(NC + NR + 1, 0.0f);
	for (int i = 1; i <= NC + NR; ++i)
		diagonal2[i] = D_tilder2[i - 1] + D_tilder2[i] + Xr2[i] * h[i];

	// bidiagonal vector define
	vector<double> bidiagonal2(NC + NR + 1, 0.0f);
	for (int i = 1; i <= NC + NR; ++i)
		bidiagonal2[i] = -D_tilder2[i];
	//////////////////////////////////////////////////////////////////////////////////////////
	
	// initial flux vector define
	for (int i = 1; i <= NC + NR; ++i) {
		phi1_old[i] = 1.0;
		phi2_old[i] = 1.0;
	}
	for (int i = 1; i <= NC + NR; ++i) {
		phi1_new[i] = 1.0;
		phi2_new[i] = 1.0;
	}

	// error
	double errk = 1.0f; // keff error
	double errp = 1.0f; // flux error
	double errp1 = 1.0f;
	double errp2 = 1.0f;

	// iteration counts
	int icounts = 1;

	while (errk > crik || errp > crip) {

		// Source vector define
		vector<double> source1(NC + NR + 1, 0.0f);
		vector<double> source2(NC + NR + 1, 0.0f);
		for (int i = 1; i <= NC + NR; ++i) {
			source1[i] = 1 / keff_old * (Xf1[i] * phi1_old[i] + Xf2[i] * phi2_old[i]);
			source2[i] = Xs12[i] * phi1_old[i];
		}

		for (int counts = 0; counts < 4; ++counts)
			for (int i = 1; i <= NC + NR; ++i) {
				phi1_new[i] = (source1[i] * h[i] - bidiagonal1[i] * phi1_old[i + 1] - bidiagonal1[i - 1] * phi1_new[i - 1]) / diagonal1[i];
				phi2_new[i] = (source2[i] * h[i] - bidiagonal2[i] * phi2_old[i + 1] - bidiagonal2[i - 1] * phi2_new[i - 1]) / diagonal2[i];
			}

		// keff calculation
		double numerator = 0.0f;
		double denominator = 0.0f;
		for (int i = 1; i <= NC; ++i) {
			numerator += Xf1[i] * phi1_new[i] + Xf2[i] * phi2_new[i]; // numerator of keff
			denominator += Xf1[i] * phi1_old[i] + Xf2[i] * phi2_old[i]; // denominator of keff
		}
		keff_new = keff_old * (numerator / denominator);

		// keff error calculation
		errk = abs((keff_new - keff_old) / keff_old);

		// flux error calculation
		for (int i = 1; i < NC + NR; i++) {
			double errp11 = abs(phi1_old[i] - phi1_new[i]) / phi1_old[i]; // [i]th error of group 1
			double errp12 = abs(phi1_old[i + 1] - phi1_new[i + 1]) / phi1_old[i + 1]; // [i+1]th error of group 1
			errp1 = max(errp11, errp12); // Max error of group 1
			double errp21 = abs(phi2_old[i] - phi2_new[i]) / phi2_old[i]; // [i]th error of group 2
			double errp22 = abs(phi2_old[i + 1] - phi2_new[i + 1]) / phi2_old[i + 1]; // [i+1]th error of group 2
			errp2 = max(errp21, errp22); // Max error of group 2
			errp = max(errp1, errp2); // Max error between group 1,2
		}

		// reset old flux
		for (int i = 1; i <= NC + NR; ++i) {
			phi1_old[i] = phi1_new[i];
			phi2_old[i] = phi2_new[i];
		}

		// reset keff
		keff_old = keff_new;

		++icounts;
	}

	// check end time
	double end_time = clock();

	// duration calculation
	double duration = (end_time - start_time) / CLOCKS_PER_SEC;

	cout << "keff : " << keff_new << ", error : " << max(errk, errp) << "\n" << "duration : " << duration << "s" << ", iteration counts : " << icounts;
}