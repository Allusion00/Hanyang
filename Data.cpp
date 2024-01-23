#include "define.h"

// 필요한 크기의 vector를 생성후 0.0으로 채운 뒤, input.file에서 정보 읽은 후 저장하는 방식

int main() {
	ifstream fileinput;
	fileinput.open("input.txt");
	string line;

	// Read CARD1 (# of Region, # of Energy group)
	cout << "[CARD1]\n";
	getline(fileinput, line);
	int region; // # of Region
	fileinput >> region;
	cout << "# of Region : " << region << "\n"; // 출력
	int g; // # of Energy group
	fileinput >> g;
	cout << "# of group : " << g << "\n\n"; // 출력
	int n = 3; // # of types of Data1

	// Read CARD2 (Boundary Condition[Left Right] : 0 = reflective, 1 = flux zero, 2 = vacuum)
	cout << "[CARD2]\n";
	getline(fileinput, line);
	getline(fileinput, line);
	int Left, Right;
	fileinput >> Left >> Right;
	double rL, rR;
	// Left Boundary coefficient
	if (Left == 0)
		rL = 10E-10f;
	else if (Left == 1)
		rL = 10E+10f;
	else
		rL = 0.5;
	// Right Boundary coefficient
	if (Right == 0)
		rR = 10E-10f;
	else if (Right == 1)
		rR = 10E+10f;
	else
		rR = 0.5;

	// 테스트용 출력
	cout << "Left Boundary Condition : " << Left << " (" << rL << ")\nRight Boundary Condition : " << Right << " (" << rR << ")\n\n";

	// Read CARD3 (Convergence criteria : # of inner iterations, relative error of keff, max relative error of flux)
	cout << "[CARD3]\n";
	getline(fileinput, line);
	getline(fileinput, line);
	int inneriter;
	double crik, crip;
	fileinput >> inneriter >> crik >> crip;

	// 테스트용 출력
	cout << "# of inner iterations : " << inneriter << "\nkeff error : " << crik << "\nsource group error : " << crip << "\n\n";

	// Read CARD4 (Size of each region)
	cout << "[CARD4]\n";
	getline(fileinput, line);
	getline(fileinput, line);
	double sizeC, sizeR;
	fileinput >> sizeC >> sizeR;
	cout << "Size of core : " << sizeC << "cm\nSize of reflector : " << sizeR << "cm\n\n";

	// Read CARD5 (# of fine meshes for each region)
	cout << "[CARD5]\n";
	getline(fileinput, line);
	getline(fileinput, line);
	int NC, NR;
	fileinput >> NC >> NR;
	cout << "# of fine meshes for core : " << NC << "\n# of fine meshes for reflector : " << NR << "\n\n";

	// Calculate Size of meshes
	vector<double> h(NC + NR + 1, 0.0f);
	for (int i = 1; i <= NC; ++i)
		h[i] = sizeC / NC;
	for (int i = NC + 1; i <= NC + NR; ++i)
		h[i] = sizeR / NR;

	// Read CARD6 (D, Xs(Fission), Xs(Absorption))
	cout << "[CARD6]\n";
	vector<vector<vector<double>>> Data1(region, vector<vector<double>>(g, vector<double>(n, 0.0f)));
	getline(fileinput, line);
	getline(fileinput, line);
	for (int i = 0; i < region; ++i) {
		for (int j = 0; j < g; ++j) {
			getline(fileinput, line);
			for (int k = 0; k < n; ++k) {
				fileinput >> Data1[i][j][k];
				cout << Data1[i][j][k] << " "; // 출력
			}
			cout << "\n"; // 출력
		}
		getline(fileinput, line);
		cout << "\n"; // 출력
	}

	// Read CARD7 (Xs(Only Down Scattering))
	cout << "[CARD7]\n";
	vector<vector<vector<double>>> Data2(region, vector<vector<double>>(g, vector<double>(g, 0.0f)));
	getline(fileinput, line);
	for (int i = 0; i < region; ++i) {
		for (int j = 0; j < g; ++j) {
			getline(fileinput, line);
			for (int k = 0; k < g; ++k) {
				fileinput >> Data2[i][j][k];
				cout << Data2[i][j][k] << " "; // 출력
			}
			cout << "\n"; // 출력
		}
		getline(fileinput, line);
		cout << "\n"; // 출력
	}

	// Read CARD8 (xg : fraction of fission neutrons)
	cout << "[CARD8]\n";
	// xg [group]
	vector<double> xg(g, 0.0f);
	getline(fileinput, line);
	for (int i = 0; i < g; ++i) {
		fileinput >> xg[i];
		cout << "group " << i + 1 << " 의 fission neutron fraction : " << xg[i] << "\n";
	}
	cout << "\n\n";

	cout << "< T E S T >\n\n"; // 테스트용 출력
	
	// Diffusion coefficient [cm] [group][mesh]
	vector<vector<double>> D(g, vector<double>(NC + NR + 1, 0.0f));
	// fission X section [cm-1] [group][mesh]
	vector<vector<double>> Xf(g, vector<double>(NC + NR + 1, 0.0f));
	// absorption X section [cm-1] [group][mesh]
	vector<vector<double>> Xa(g, vector<double>(NC + NR + 1, 0.0f));
	for (int i = 0; i < g; ++i) {
		for (int j = 1; j <= NC; ++j) {
			D[i][j] = Data1[0][i][0];
			Xf[i][j] = Data1[0][i][1];
			Xa[i][j] = Data1[0][i][2];
		}
		for (int j = NC + 1; j <= NC + NR; ++j) {
			D[i][j] = Data1[1][i][0];
			Xf[i][j] = Data1[1][i][1];
			Xa[i][j] = Data1[1][i][2];
		}
	}

	// 테스트용 출력
	for (int i = 0; i < g; ++i) {
		cout << "group " << i + 1 << " 의 D : ";
		for (int j = 0; j <= NC + NR; ++j) {
			cout << D[i][j] << " ";
		}
		cout << "\ngroup " << i + 1 << " 의 Xf : ";
		for (int j = 0; j <= NC + NR; ++j) {
			cout << Xf[i][j] << " ";
		}
		cout << "\ngroup " << i + 1 << " 의 Xa : ";
		for (int j = 0; j <= NC + NR; ++j) {
			cout << Xa[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";

	// scattering X section [cm-1] [high E group][low E group][mesh]
	vector<vector<vector<double>>> Xs(g, vector<vector<double>>(g, vector<double>(NC + NR + 1, 0.0f)));
	for (int i = 0; i < g; ++i) {
		for (int j = 0; j < g; ++j) {
			for (int k = 1; k <= NC; ++k)
				Xs[i][j][k] = Data2[0][i][j];
			for (int k = NC + 1; k <= NC + NR; ++k)
				Xs[i][j][k] = Data2[1][i][j];
		}
	}

	// 테스트용 출력
	for (int i = 0; i < g; ++i) {
		for (int j = 0; j < g; ++j) {
			for (int k = 0; k <= NC + NR; ++k) {
				cout << Xs[i][j][k] << " ";
			}
			cout << "\n";
		}
		cout << "\n";
	}

	// removal X section [cm-1] [group][mesh]
	vector<vector<double>> Xr(g, vector<double>(NC + NR + 1, 0.0f));
	vector<vector<double>> Xs_sum(g, vector<double>(NC + NR + 1, 0.0f));
	for (int i = 0; i < g; ++i)
		for (int j = 1; j <= NC + NR; ++j) {
			for (int k = 0; k < g; ++k)
				Xs_sum[i][j] += Xs[i][k][j];
			Xr[i][j] = Xa[i][j] + Xs_sum[i][j];
		}

	// 테스트용 출력
	cout << "removal X section\n";
	for (int i = 0; i < g; ++i) {
		for (int j = 0; j <= NC + NR; ++j) {
			cout << Xr[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "\n";

	// initial effective k
	double keff_old = 1.0f;
	double keff_new = 0.0f;

	// flux vector define [group][mesh]
	vector<vector<double>> phi_old(g, vector<double>(NC + NR + 2, 0.0f));
	vector<vector<double>> phi_new(g, vector<double>(NC + NR + 2, 0.0f));
	for (int i = 0; i < g; ++i) {
		for (int j = 1; j <= NC + NR; ++j) {
			phi_old[i][j] = 1.0f;
			phi_new[i][j] = 1.0f;
		}
	}

	// beta vector define [group][mesh]
	vector<vector<double>> beta(g, vector<double>(NC + NR + 2, 0.0f));
	for (int i = 0; i < g; ++i) {
		beta[i][0] = rL / 2;
		for (int j = 1; j <= NC + NR; ++j)
			beta[i][j] = D[i][j] / h[j];
		beta[i][NC + NR + 1] = rR / 2;
	}

	// D tilder vector define [group][mesh]
	vector<vector<double>> D_tilder(g, vector<double>(NC + NR + 1, 0.0f));
	for (int i = 0; i < g; ++i)
		for (int j = 0; j <= NC + NR; ++j)
			D_tilder[i][j] = 2 * beta[i][j] * beta[i][j + 1] / (beta[i][j] + beta[i][j + 1]);

	// Matrix M define diagonal, bidiagonal [group][mesh]
	vector<vector<double>> diagonal(g, vector<double>(NC + NR + 1, 0.0f));
	vector<vector<double>> bidiagonal(g, vector<double>(NC + NR + 1, 0.0f));
	for (int i = 0; i < g; ++i)
		for (int j = 1; j <= NC + NR; ++j) {
			diagonal[i][j] = D_tilder[i][j - 1] + D_tilder[i][j] + Xr[i][j] * h[j];
			bidiagonal[i][j] = -D_tilder[i][j];
		}

	// initial error
	double errk = 1.0f; // keff error
	double errp = 1.0f; // flux error
	double errp1 = 1.0f;
	double errp2 = 1.0f;

	// Source vector define [high E group][low E group][mesh] -> while문 안에 들어가야함
	vector<vector<double>> source(g, vector<double>(NC + NR + 1, 0.0f));
	vector<vector<double>> fission(g, vector<double>(NC + NR + 1, 0.0f));
	vector<vector<double>> scattering(g, vector<double>(NC + NR + 1, 0.0f));
	
	for (int i = 0; i < g; ++i) {
		for (int j = 1; j <= NC + NR; ++j) {
			for (int k = 0; k < g; ++k) {
				fission[i][j] += 1 / keff_old * Xf[k][j] * phi_old[k][j] * h[j] * xg[i];
				scattering[i][j] += Xs[i][k][j] * phi_old[k][j] * h[j];
			}
			source[i][j] = fission[i][j] + scattering[i][j];
			cout << source[i][j] << " ";
		}
		cout << "\n";
	}
}
