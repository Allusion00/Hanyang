//#include <iostream>
//#include <vector>
//
//using namespace std;
//
//// check start time
//double start_time = clock();
//
//// mesh number
//int NC = 80;
//int NR = 10;
//
//// Diffusion coefficient [cm]
//double DC = 6.8f;
//double DR = 5.8f;
//
//// absorption X section [cm-1]
//double XaC = 0.300f;
//double XaR = 0.01f;
//
//// fission X section [cm-1]
//double XfC = 0.65f;
//double XfR = 0.0f;
//
//// region size [cm]
//double sizeC = 40.0f;
//double sizeR = 10.0f;
//
//// mesh size [cm]
//double hC = sizeC / NC;
//double hR = sizeR / NR;
//
//// initial effective k
//double keff_old = 1.0f;
//double keff_new = 0.0f;
//
//// for effective k calculation
//double phi_sum_old = 0.0f;
//double phi_sum_new = 0.0f;
//
//// flux vector
//vector<double> phi_old(NC + NR + 2, 0.0f);
//vector<double> phi_new(NC + NR + 2, 0.0f);
//
//// Boundary condition
///*	vacuum condition : r = 0.5
//	current zero (reflective) condition : r = 0; 10E-10
//	flux zero condition : r = inf; 10E+10
//string BL = "reflective"; // Left boundary
//string BR = "flux zero"; // Right boundary */
//
//// Boundary coefficient, r
//double rL = 10E-10f;
//double rR = 10E+10f;
//
//// Convergence criteria
//double crik = 1.0E-10f;
//double crip = 1.0E-10f;
//
//
//int main() {
//
//	// Diffusion coefficient vector define
//	vector<double> D(NC + NR + 1, 0.0f);
//	// mesh size vector define
//	vector<double> h(NC + NR + 1, 0.0f);
//	// absorption X section vector define
//	vector<double> Xa(NC + NR + 1, 0.0f);
//	// fission X section vector define
//	vector<double> Xf(NC + NR + 1, 0.0f);
//	for (int i = 1; i <= NC; ++i) {
//		D[i] = DC;
//		h[i] = hC;
//		Xa[i] = XaC;
//		Xf[i] = XfC;
//	}
//	for (int i = NC + 1; i <= NC + NR; ++i) {
//		D[i] = DR;
//		h[i] = hR;
//		Xa[i] = XaR;
//		Xf[i] = XfR;
//	}
//				
//	// beta vector define
//	vector<double> beta(NC + NR + 2, 0.0f);
//	beta[0] = rL / 2;
//	for (int i = 1; i <= NC + NR; ++i)
//		beta[i] = D[i] / h[i];
//	beta[NC + NR + 1] = rR / 2;
//
//	// D tilder vector define
//	vector<double> D_tilder(NC + NR + 1, 0.0f);
//	for (int i = 0; i <= NC + NR; ++i)
//		D_tilder[i]= 2 * beta[i] * beta[i + 1] / (beta[i] + beta[i + 1]);
//
//	// diagonal vector define
//	vector<double> diagonal(NC + NR + 1, 0.0f);
//	for (int i = 1; i <= NC + NR; ++i)
//		diagonal[i] = D_tilder[i - 1] + D_tilder[i] + Xa[i] * h[i];
//
//	// bidiagonal vector define
//	vector<double> bidiagonal(NC + NR + 1, 0.0f);
//	for (int i = 1; i <= NC + NR; ++i)
//		bidiagonal[i] = -D_tilder[i];
//
//	
//	// initial flux vector define
//	for (int i = 1; i <= NC + NR; ++i)
//		phi_old[i] = 1.0;
//	for (int i = 1; i <= NC + NR; ++i)
//		phi_new[i] = 1.0;
//
//	// error
//	double errk = 1.0f; // keff error
//	double errp = 1.0f; // flux error
//
//	// iteration counts
//	int icounts = 1;
//
//	while (errk > crik || errp > crip) {
//
//		// old flux sum
//		for (int i = 1; i <= NC; ++i)
//			phi_sum_old += phi_old[i];
//	
//		// Source vector define
//		vector<double> source(NC + NR + 1, 0.0f);
//		for (int i = 1; i <= NC + NR; ++i)
//			source[i] = 1 / keff_old * Xf[i] * phi_old[i];
//
//		for (int counts = 0; counts < 5; ++counts)
//			for (int i = 1; i <= NC + NR; ++i)
//				phi_new[i] = (source[i] * h[i] - bidiagonal[i] * phi_old[i + 1] - bidiagonal[i - 1] * phi_new[i - 1]) / diagonal[i];
//		
//		// new flux sum 
//		for (int i = 1; i <= NC; ++i)
//			phi_sum_new += phi_new[i];
//
//		// keff calculation
//		keff_new = keff_old * (phi_sum_new / phi_sum_old);
//
//		// keff error calculation
//		errk = abs((keff_new - keff_old) / keff_old);
//
//		// flux error calculation
//		for (int i = 1; i < NC + NR; i++) {
//			double err1 = abs(phi_old[i] - phi_new[i]) / phi_old[i]; // [i]th error
//			double err2 = abs(phi_old[i + 1] - phi_new[i + 1]) / phi_old[i + 1]; // [i+1]th error
//			errp = max(err1, err2); // Max error
//		}
//
//		// reset old flux
//		for (int i = 1; i <= NC + NR; ++i)
//			phi_old[i] = phi_new[i];
//
//		// reset keff
//		keff_old = keff_new;
//		
//		// reset flux sum
//		phi_sum_new = 0.0f;
//		phi_sum_old = 0.0f;
//
//		++icounts;
//		//return keff_new;
//	}
//
//	// check end time
//	double end_time = clock();
//
//	// duration calculation
//	double duration = (end_time - start_time) / CLOCKS_PER_SEC;
//
//	cout << "keff : " << keff_new << ", error : " << max(errk,errp) <<  "\n" << "duration : " << duration << "s" << ", iteration counts : " << icounts;
//}