//#include <iostream>
//#include <vector>
//
//using namespace std;
//
//// check start time
//double start_time = clock();
//
//// mesh number
//int NCi = 80;
//int NRi = 20;
//int Nj = 100;
//
//// Diffusion coefficient [cm]
//double D1C = 1.3314f;
//double D2C = 0.2042f;
//double D1R = 1.2291f;
//double D2R = 0.1588f;
//
//// absorption X section [cm-1]
//double Xa1C = 0.006431f;
//double Xa2C = 0.09992f;
//double Xa1R = 0.0005903f;
//double Xa2R = 0.01882f;
//
//// fission X section [cm-1]
//double Xf1C = 0.003258f;
//double Xf2C = 0.1387f;
//double Xf1R = 0.00f;
//double Xf2R = 0.00f;
//
//// scattering X section [cm-1]
//double Xs12C = 0.05979f;
//double Xs12R = 0.06163f;
//
//// region size [cm]
//double sizeC = 40.0f;
//double sizeR = 10.0f;
//double region = 50.0f;
//
//// mesh size [cm]
//double hCi = sizeC / NCi;
//double hRi = sizeR / NRi;
//double hmj = region / Nj; // j 방향 mesh 개수
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
//vector<vector<double>> phi1_old(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f)); // group 1 flux
//vector<vector<double>> phi1_new(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
//vector<vector<double>> phi2_old(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f)); // group 1 flux
//vector<vector<double>> phi2_new(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
//
//// Boundary coefficient, r
//double rN = 10E-10f; // reflective
//double rS = 10E-10f; // reflective
//double rW = 10E-10f; // reflective
//double rE = 10E+10f; // flux zero
//
//// Convergence criteria
//double crik = 1.0E-10f;
//double crip = 1.0E-07f;
//
//
//int main() {
//
//	// Diffusion coefficient vector define
//	vector<vector<double>> D1(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> D2(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	// mesh size vector define
//	vector<double> hi(NCi + NRi + 1, 0.0f); // i 방향 mesh size
//	vector<double> hj(Nj + 1, 0.0f); // j 방향 mesh size
//	// absorption X section vector define
//	vector<vector<double>> Xa1(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> Xa2(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	// fission X section vector define
//	vector<vector<double>> Xf1(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> Xf2(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	// scattering X section vector define
//	vector<vector<double>> Xs12(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	// removal X section vector define
//	vector<vector<double>> Xr1(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> Xr2(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	for (int i = 1; i <= NCi; ++i)
//		for (int j = 1; j <= Nj; ++j) {
//			D1[i][j] = D1C;
//			D2[i][j] = D2C;
//			hi[i] = hCi;
//			hj[j] = hmj;
//			Xa1[i][j] = Xa1C;
//			Xa2[i][j] = Xa2C;
//			Xf1[i][j] = Xf1C;
//			Xf2[i][j] = Xf2C;
//			Xs12[i][j] = Xs12C;
//			Xr1[i][j] = Xa1[i][j] + Xs12[i][j];
//			Xr2[i][j] = Xa2[i][j];
//		}
//	for (int i = NCi + 1; i <= NCi + NRi; ++i)
//		for (int j = 1; j <= Nj; ++j) {
//			D1[i][j] = D1R;
//			D2[i][j] = D2R;
//			hi[i] = hRi;
//			hj[j] = hmj;
//			Xa1[i][j] = Xa1R;
//			Xa2[i][j] = Xa2R;
//			Xf1[i][j] = Xf1R;
//			Xf2[i][j] = Xf2R;
//			Xs12[i][j] = Xs12R;
//			Xr1[i][j] = Xa1[i][j] + Xs12[i][j];
//			Xr2[i][j] = Xa2[i][j];
//		}
//
//	// beta vector define
//	vector<vector<double>> beta1i(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
//	vector<vector<double>> beta1j(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
//	vector<vector<double>> beta2i(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
//	vector<vector<double>> beta2j(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
//	for (int j = 1; j <= Nj; ++j) {
//		beta1i[0][j] = rW / 2; // i 방향 1그룹 beta의 가장 서쪽 boundary 계수
//		beta2i[0][j] = rW / 2; // i 방향 2그룹 beta의 가장 서쪽 boundary 계수
//	}
//	for (int i = 1; i <= NCi + NRi; ++i) {
//		beta1j[i][0] = rN / 2; // j 방향 1그룹 beta의 가장 북쪽 boundary 계수
//		beta2j[i][0] = rN / 2; // j 방향 2그룹 beta의 가장 북쪽 boundary 계수
//	}
//	for (int i = 1; i <= NCi + NRi; ++i)
//		for (int j = 1; j <= Nj; ++j) {
//			beta1i[i][j] = D1[i][j] / hi[i]; // 1그룹 beta 정의
//			beta1j[i][j] = D1[i][j] / hj[j];
//			beta2i[i][j] = D2[i][j] / hi[i]; // 2그룹 beta 정의
//			beta2j[i][j] = D2[i][j] / hj[j];
//		}
//	for (int j = 1; j <= Nj; ++j) {
//		beta1i[NCi + NRi + 1][j] = rE / 2; // i 방향 1그룹 beta의 가장 동쪽 boundary 계수
//		beta2i[NCi + NRi + 1][j] = rE / 2; // i 방향 2그룹 beta의 가장 동쪽 boundary 계수
//	}
//	for (int i = 1; i <= NCi + NRi; ++i) {
//		beta1j[i][Nj + 1] = rS / 2; // j 방향 1그룹 beta의 가장 남쪽 boundary 계수
//		beta2j[i][Nj + 1] = rS / 2; // j 방향 2그룹 beta의 가장 남쪽 boundary 계수
//	}
//
//	// D tilder vector define
//	vector<vector<double>> D_tilder1E(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> D_tilder1W(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> D_tilder1N(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> D_tilder1S(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> D_tilder2E(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> D_tilder2W(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> D_tilder2N(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> D_tilder2S(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	for (int i = 1; i <= NCi + NRi; ++i)
//		for (int j = 1; j <= Nj; ++j) {
//			D_tilder1E[i][j] = 2 * beta1i[i][j] * beta1i[i + 1][j] / (beta1i[i][j] + beta1i[i + 1][j]);
//			D_tilder1W[i][j] = 2 * beta1i[i][j] * beta1i[i - 1][j] / (beta1i[i][j] + beta1i[i - 1][j]);
//			D_tilder1N[i][j] = 2 * beta1j[i][j] * beta1j[i][j - 1] / (beta1j[i][j] + beta1j[i][j - 1]);
//			D_tilder1S[i][j] = 2 * beta1j[i][j] * beta1j[i][j + 1] / (beta1j[i][j] + beta1j[i][j + 1]);
//			D_tilder2E[i][j] = 2 * beta2i[i][j] * beta2i[i + 1][j] / (beta2i[i][j] + beta2i[i + 1][j]);
//			D_tilder2W[i][j] = 2 * beta2i[i][j] * beta2i[i - 1][j] / (beta2i[i][j] + beta2i[i - 1][j]);
//			D_tilder2N[i][j] = 2 * beta2j[i][j] * beta2j[i][j - 1] / (beta2j[i][j] + beta2j[i][j - 1]);
//			D_tilder2S[i][j] = 2 * beta2j[i][j] * beta2j[i][j + 1] / (beta2j[i][j] + beta2j[i][j + 1]);
//		}
//
//
//	// diagonal vector define
//	vector<vector<double>> diagonal1(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> diagonal2(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	for (int i = 1; i <= NCi + NRi; ++i)
//		for (int j = 1; j <= Nj; ++j) {
//			diagonal1[i][j] = D_tilder1E[i][j] / hi[i] + D_tilder1W[i][j] / hi[i] + D_tilder1N[i][j] / hj[j] + D_tilder1S[i][j] / hj[j] + Xr1[i][j];
//			diagonal2[i][j] = D_tilder2E[i][j] / hi[i] + D_tilder2W[i][j] / hi[i] + D_tilder2N[i][j] / hj[j] + D_tilder2S[i][j] / hj[j] + Xr2[i][j];
//		}
//	
//	// bidiagonal vector define
//	vector<vector<double>> bidiagonal1E(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> bidiagonal1W(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> bidiagonal1N(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> bidiagonal1S(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> bidiagonal2E(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> bidiagonal2W(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> bidiagonal2N(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	vector<vector<double>> bidiagonal2S(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//	for (int i = 1; i <= NCi + NRi; ++i)
//		for (int j = 1; j <= Nj; ++j) {
//			bidiagonal1E[i][j] = -D_tilder1E[i][j] / hi[i];
//			bidiagonal1W[i][j] = -D_tilder1W[i][j] / hi[i];
//			bidiagonal1N[i][j] = -D_tilder1N[i][j] / hj[j];
//			bidiagonal1S[i][j] = -D_tilder1S[i][j] / hj[j];
//			bidiagonal2E[i][j] = -D_tilder2E[i][j] / hi[i];
//			bidiagonal2W[i][j] = -D_tilder2W[i][j] / hi[i];
//			bidiagonal2N[i][j] = -D_tilder2N[i][j] / hj[j];
//			bidiagonal2S[i][j] = -D_tilder2S[i][j] / hj[j];
//		}
//
//	// initial flux vector define
//	for (int i = 1; i <= NCi + NRi; ++i)
//		for (int j = 1; j <= Nj; ++j) {
//			phi1_old[i][j] = 1.0;
//			phi1_new[i][j] = 1.0;
//			phi2_old[i][j] = 1.0;
//			phi2_new[i][j] = 1.0;
//		}
//
//	// error
//	double errk = 1.0f; // keff error
//	double errp = 1.0f; // flux error
//	double errp1 = 1.0f;
//	double errp2 = 1.0f;
//
//	// iteration counts
//	int icounts = 1;
//
//	while (errk > crik || errp > crip) {
//
//		// Source vector define
//		vector<vector<double>> source1(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//		vector<vector<double>> source2(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
//		for (int i = 1; i <= NCi + NRi; ++i)
//			for (int j = 1; j <= Nj; ++j) {
//				source1[i][j] = 1 / keff_old * (Xf1[i][j] * phi1_old[i][j] + Xf2[i][j] * phi2_old[i][j]);
//				source2[i][j] = Xs12[i][j] * phi1_old[i][j];
//			}
//		
//		// inner iteration (Gauss Seidel Method)
//		for (int counts = 0; counts < 4; ++counts)
//			for (int i = 1; i <= NCi + NRi; ++i)
//				for (int j = 1; j <= Nj; ++j) {
//					phi1_new[i][j] = (source1[i][j] - bidiagonal1E[i][j] * phi1_new[i + 1][j] - bidiagonal1W[i][j] * phi1_new[i - 1][j] - bidiagonal1N[i][j] * phi1_new[i][j - 1] - bidiagonal1S[i][j] * phi1_new[i][j + 1]) / diagonal1[i][j];
//					phi2_new[i][j] = (source2[i][j] - bidiagonal2E[i][j] * phi2_new[i + 1][j] - bidiagonal2W[i][j] * phi2_new[i - 1][j] - bidiagonal2N[i][j] * phi2_new[i][j - 1] - bidiagonal2S[i][j] * phi2_new[i][j + 1]) / diagonal2[i][j];
//				}
//
//		// keff calculation
//		double numerator = 0.0f;
//		double denominator = 0.0f;
//		for (int i = 1; i <= NCi; ++i)
//			for (int j = 1; j <= Nj; ++j) {
//				numerator += Xf1[i][j] * phi1_new[i][j] + Xf2[i][j] * phi2_new[i][j]; // numerator of keff
//				denominator += Xf1[i][j] * phi1_old[i][j] + Xf2[i][j] * phi2_old[i][j]; // denominator of keff
//			}
//		keff_new = keff_old * (numerator / denominator);
//
//		// keff error calculation
//		errk = abs((keff_new - keff_old) / keff_old);
//
//		// flux error calculation
//		for (int i = 1; i < NCi + NRi; ++i)
//			for (int j = 1; j < Nj; ++j) {
//				double errp11 = abs(phi1_old[i][j] - phi1_new[i][j]) / phi1_old[i][j]; // [i]th error of group 1
//				double errp12 = abs(phi1_old[i + 1][j] - phi1_new[i + 1][j]) / phi1_old[i + 1][j]; // [i+1]th error of group 1
//				errp1 = max(errp11, errp12); // Max error of group 1
//				double errp21 = abs(phi2_old[i][j] - phi2_new[i][j]) / phi2_old[i][j]; // [i]th error of group 2
//				double errp22 = abs(phi2_old[i + 1][j] - phi2_new[i + 1][j]) / phi2_old[i + 1][j]; // [i+1]th error of group 2
//				errp2 = max(errp21, errp22); // Max error of group 2
//				errp = max(errp1, errp2); // Max error between group 1,2
//			}
//
//		// reset old flux
//		for (int i = 1; i <= NCi + NRi; ++i)
//			for (int j = 1; j <= Nj; ++j) {
//				phi1_old[i][j] = phi1_new[i][j];
//				phi2_old[i][j] = phi2_new[i][j];
//			}
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
//	double duration = (end_time - start_time) / CLOCKS_PER_SEC * 1000;
//
//	cout << "keff : " << keff_new << ", error : " << min(errk, errp) << "\n" << "duration : " << duration << "ms" << ", iteration counts : " << icounts;
//}