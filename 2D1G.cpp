#include <iostream>
#include <vector>

using namespace std;

// check start time
double start_time = clock();

// mesh number
int NCi = 80;
int NRi = 20;
int Nj = 100;

// Diffusion coefficient [cm]
double DC = 6.8f;
double DR = 5.8f;

// absorption X section [cm-1]
double XaC = 0.300f;
double XaR = 0.01f;

// fission X section [cm-1]
double XfC = 0.65f;
double XfR = 0.0f;

// region size [cm]
double sizeC = 40.0f;
double sizeR = 10.0f;
double region = 50.0f;

// mesh size [cm]
double hCi = sizeC / NCi;
double hRi = sizeR / NRi;
double hmj = region / Nj; // j 방향 mesh 개수

// initial effective k
double keff_old = 1.0f;
double keff_new = 0.0f;

// for effective k calculation
double phi_sum_old = 0.0f;
double phi_sum_new = 0.0f;

// flux vector
vector<vector<double>> phi_old(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
vector<vector<double>> phi_new(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));

// Boundary coefficient, r
double rN = 10E-10f; // reflective
double rS = 10E-10f; // reflective
double rW = 10E-10f; // reflective
double rE = 10E+10f; // flux zero

// Convergence criteria
double crik = 1.0E-10f;
double crip = 1.0E-07f;


int main() {

	// Diffusion coefficient vector define
	vector<vector<double>> D(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
	// mesh size vector define
	vector<double> hi(NCi + NRi + 1, 0.0f); // i 방향 mesh size
	vector<double> hj(Nj + 1, 0.0f); // j 방향 mesh size
	// absorption X section vector define
	vector<vector<double>> Xa(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
	// fission X section vector define
	vector<vector<double>> Xf(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
	for (int i = 1; i <= NCi; ++i)		
		for (int j = 1; j <= Nj; ++j) {
			D[i][j] = DC;
			hi[i] = hCi;
			hj[j] = hmj;
			Xa[i][j] = XaC;
			Xf[i][j] = XfC;
		}
	for (int i = NCi + 1; i <= NCi + NRi; ++i)
		for (int j = 1; j <= Nj; ++j) {
			D[i][j] = DR; 
			hi[i] = hRi;
			hj[j] = hmj;
			Xa[i][j] = XaR;
			Xf[i][j] = XfR;
		}

	// beta vector define
	vector<vector<double>> betai(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
	vector<vector<double>> betaj(NCi + NRi + 2, vector<double>(Nj + 2, 0.0f));
	for (int j = 0; j <= Nj; ++j)
		betai[0][j] = rW / 2; // i 방향 beta의 가장 서쪽 boundary 계수
	for (int i = 0; i <= NCi + NRi; ++i)
		betaj[i][0] = rN / 2; // j 방향 beta의 가장 북쪽 boundary 계수
	for (int i = 1; i <= NCi + NRi; ++i)
		for (int j = 1; j <= Nj; ++j) {
			betai[i][j] = D[i][j] / hi[i];
			betaj[i][j] = D[i][j] / hj[j];
		}
	for (int j = 0; j <= Nj; ++j)
		betai[NCi + NRi + 1][j] = rE / 2;
	for (int i = 0; i <= NCi + NRi; ++i)
		betaj[i][Nj + 1] = rS / 2;
	
	// D tilder vector define
	vector<vector<double>> D_tilderi(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
	vector<vector<double>> D_tilderj(NCi + NRi + 1, vector<double>(Nj + 1, 0.0f));
	for (int i = 0; i <= NCi + NRi; ++i)
		for (int j = 0; j <= Nj; ++j) {
			D_tilderi[i][j] = 2 * betai[i][j] * betai[i + 1][j] / (betai[i][j] + betai[i + 1][j]);
			D_tilderj[i][j] = 2 * betaj[i][j] * betaj[i][j + 1] / (betaj[i][j] + betaj[i][j + 1]);
		}
		
	// diagonal vector define
	vector<double> diagonal(NC + NR + 1, 0.0f);
	for (int i = 1; i <= NC + NR; ++i)
		diagonal[i] = D_tilder[i - 1] + D_tilder[i] + Xa[i] * h[i];

	// bidiagonal vector define
	vector<double> bidiagonal(NC + NR + 1, 0.0f);
	for (int i = 1; i <= NC + NR; ++i)
		bidiagonal[i] = -D_tilder[i];


	// initial flux vector define
	for (int i = 1; i <= NC + NR; ++i)
		phi_old[i] = 1.0;
	for (int i = 1; i <= NC + NR; ++i)
		phi_new[i] = 1.0;

	// error
	double errk = 1.0f; // keff error
	double errp = 1.0f; // flux error

	// iteration counts
	int icounts = 1;

	while (errk > crik || errp > crip) {

		// old flux sum
		for (int i = 1; i <= NC; ++i)
			phi_sum_old += phi_old[i];

		// Source vector define
		vector<double> source(NC + NR + 1, 0.0f);
		for (int i = 1; i <= NC + NR; ++i)
			source[i] = 1 / keff_old * Xf[i] * phi_old[i];

		for (int counts = 0; counts < 5; ++counts)
			for (int i = 1; i <= NC + NR; ++i)
				phi_new[i] = (source[i] * h[i] - bidiagonal[i] * phi_new[i + 1] - bidiagonal[i - 1] * phi_new[i - 1]) / diagonal[i];

		// new flux sum 
		for (int i = 1; i <= NC; ++i)
			phi_sum_new += phi_new[i];

		// keff calculation
		keff_new = keff_old * (phi_sum_new / phi_sum_old);

		// keff error calculation
		errk = abs((keff_new - keff_old) / keff_old);

		// flux error calculation
		for (int i = 1; i < NC + NR; i++) {
			double err1 = abs(phi_old[i] - phi_new[i]) / phi_old[i]; // [i]th error
			double err2 = abs(phi_old[i + 1] - phi_new[i + 1]) / phi_old[i + 1]; // [i+1]th error
			errp = max(err1, err2); // Max error
		}

		// reset old flux
		for (int i = 1; i <= NC + NR; ++i)
			phi_old[i] = phi_new[i];

		// reset keff
		keff_old = keff_new;

		// reset flux sum
		phi_sum_new = 0.0f;
		phi_sum_old = 0.0f;

		++icounts;
		//return keff_new;
	}

	// check end time
	double end_time = clock();

	// duration calculation
	double duration = (end_time - start_time) / CLOCKS_PER_SEC;

	cout << "keff : " << keff_new << ", error : " << min(errk, errp) << "\n" << "duration : " << duration << "s" << ", iteration counts : " << icounts;
}