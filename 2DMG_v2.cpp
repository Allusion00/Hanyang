#include <iostream>		// I/O연산에 필요한 기본적인 기능
#include <fstream>		// 파일 입출력 할 수 있게 만드는 라이브러리
#include <string>
#include <vector>
using namespace std;	// 이름의 충돌을 막아주고, 프로그램의 가독성을 높여주는 코드
// Visual 2008부터는 iostream.h를 기본으로 지원하지 않기 때문에 #include <iostream>을 하거나 cout를 쓰려면 이 코드를 써야한다.

// 필요한 크기의 vector를 생성후 0.0으로 채운 뒤, input.file에서 정보 읽은 후 저장하는 방식
// check start time
double start_time = clock();

int main() {
	ifstream fileinput; // define ifstream as fileinput

		// Select group by cin
		int group;
		cout << "Select group between 1, 2, 7, 0 (C5G7), 22 (Kang) : ";
		cin >> group;
		if (group == 1)
			fileinput.open("2D1G input.txt");
		else if (group == 2)
			fileinput.open("2D2G input.txt");
		else if (group == 7)
			fileinput.open("2D7G input.txt");
		else if (group == 0)
			fileinput.open("C5G7 input.txt");
		else if (group == 22)
			fileinput.open("Kang input.txt");
		else {
			cout << "Please select group between 1, 2, 7";
			return 0;
		}

	string line; // define line as string character for getline

	// Read CARD1 (# of Region, # of Energy group)
	cout << "\n[CARD1]\n";
	getline(fileinput, line);
	int region , g; // # of Region, # of Energy group
	fileinput >> region >> g;
	cout << "# of Region : " << region << "\n" << "# of group : " << g << "\n\n"; // output for checking CARD1

	// Read CARD2 (Geometry size [row column], Geometry)
	cout << "[CARD2]\n";
	getline(fileinput, line);
	getline(fileinput, line);
	double gridsizei, gridsizej;
	fileinput >> gridsizei >> gridsizej;
	int row, column;
	fileinput >> row >> column;
	vector<vector<int>> grid(row, vector<int>(column, 0));
	for (int i = 0; i < row; ++i) {
		getline(fileinput, line);
		for (int j = 0; j < column; ++j) {
			fileinput >> grid[j][i];
			cout << grid[j][i] << " ";
		}
		cout << "\n";
	}
	
	//// Read CARD3 (# of fine meshes for each grid [horizontal \n vertical]
	//getline(fileinput, line);
	//getline(fileinput, line);
	//vector<int> horizontal(column, 0);
	//vector<int> vertical(row, 0);
	//int Ni, Nj;
	//for (int i = 0; i < column; ++i)
	//	Ni += gridsizei / horizontal[i];
	//for (int j = 0; j < row; ++j)
	//	Nj += gridsizej / vertical[j];
	
	// Read CARD3 (# of fine meshes for each grid [row column])
	cout << "\n[CARD3]\n";
	getline(fileinput, line);
	getline(fileinput, line);
	double meshnumi, meshnumj;
	fileinput >> meshnumj >> meshnumi;
	double meshsizei = gridsizei / meshnumi;
	double meshsizej = gridsizej / meshnumj;
	int Ni = meshnumi * column;
	int Nj = meshnumj * row;
	cout << "Mesh zone type [row X column] : " << Nj << " X " << Ni << "\n";
	
	// Read CARD4 (Boundary Condition[Left Right] : 0 = reflective, 1 = flux zero, 2 = vacuum)
	cout << "\n[CARD4]\n";
	getline(fileinput, line);
	getline(fileinput, line);
	vector<int> BC(4, 0);
	for (int i = 0; i < 4; ++i)
		fileinput >> BC[i];
	vector<double> r(4, 0.0f);
	for (int i = 0; i < 4; ++i) {
		if (BC[i] == 0)
			r[i] = 10E-10f; // reflective
		else if (BC[i] == 1)
			r[i] = 10E+10f; // flux zero
		else
			r[i] = 0.5; // vacuum
	}
	vector<string> cardinal{ "North", "South", "West", "East" };
	// output for checking CARD2
	for (int i = 0; i < 4; ++i)
		cout << cardinal[i] << " Boundary Condition : " << BC[i] << " (" << r[i] << ")\n\n";

	// Read CARD5 (Convergence criteria : # of inner iterations, relative error of keff, max relative error of flux)
	cout << "[CARD5]\n";
	getline(fileinput, line);
	getline(fileinput, line);
	int inneriter;
	double crik, crip;
	fileinput >> inneriter >> crik >> crip;
	// output for checking CARD3
	cout << "# of inner iterations : " << inneriter << "\nkeff error : " << crik << "\nsource group error : " << crip << "\n\n";

	// Read CARD6 (D, Xs(Fission), Xs(Absorption)) [region][group][data type]
	cout << "[CARD6]\n";
	int n = 3; // # of types of Data1
	vector<vector<vector<double>>> Data1(region, vector<vector<double>>(g, vector<double>(n, 0.0f)));
	getline(fileinput, line);
	getline(fileinput, line);
	for (int i = 0; i < region; ++i) {
		for (int j = 0; j < g; ++j) {
			getline(fileinput, line);
			for (int k = 0; k < n; ++k) {
				fileinput >> Data1[i][j][k];
				cout << Data1[i][j][k] << " "; // output for checking CARD6
			}
			cout << "\n"; // output for checking CARD6
		}
		getline(fileinput, line);
		cout << "\n"; // output for checking CARD6
	}


	// Read CARD7 (Xs(Only Down Scattering)) [region][High E group][Low E group]
	cout << "[CARD7]\n";
	vector<vector<vector<double>>> Data2(region, vector<vector<double>>(g, vector<double>(g, 0.0f)));
	getline(fileinput, line);
	for (int i = 0; i < region; ++i) {
		for (int j = 0; j < g; ++j) {
			getline(fileinput, line);
			for (int k = 0; k < g; ++k) {
				fileinput >> Data2[i][k][j];
				cout << Data2[i][j][k] << " "; // output for checking CARD7
			}
			cout << "\n"; // output for checking CARD7
		}
		getline(fileinput, line);
		cout << "\n"; // output for checking CARD7
	}

	// Read CARD8 (xg : fraction of fission neutrons)
	cout << "[CARD8]\n";
	// xg [group]
	vector<double> xg(g, 0.0f);
	getline(fileinput, line);
	for (int i = 0; i < g; ++i) {
		fileinput >> xg[i];
		cout << "group " << i + 1 << " 의 fission neutron fraction : " << xg[i] << "\n"; // output for checking CARD8
	}
	cout << "\n\n"; // output for checking CARD8

	// Diffusion coefficient [cm] [group][i mesh][j mesh]
	vector<vector<vector<double>>> D(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	// fission X section [cm-1] [group][i mesh][j mesh]
	vector<vector<vector<double>>> Xf(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	// absorption X section [cm-1] [group][i mesh][j mesh]
	vector<vector<vector<double>>> Xa(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	for (int group = 0; group < g; ++group)
		for (int i = 0; i < column; ++i)
			for (int imesh = i * meshnumi + 1; imesh <= (i + 1) * meshnumi; ++imesh)
				for (int j = 0; j < row; ++j)
					for (int jmesh = j * meshnumj + 1; jmesh <= (j + 1) * meshnumj; ++jmesh) {
						D[group][imesh][jmesh] = Data1[grid[i][j]][group][0];
						Xf[group][imesh][jmesh] = Data1[grid[i][j]][group][1];
						Xa[group][imesh][jmesh] = Data1[grid[i][j]][group][2];
					}

	// scattering X section [cm-1] [high E group][low E group][i mesh][j mesh]
	vector<vector<vector<vector<double>>>> Xs(g, vector<vector<vector<double>>>(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f))));
	for (int hgroup = 0; hgroup < g; ++hgroup)
		for (int lgroup = 0; lgroup < g; ++lgroup)
			for (int i = 0; i < column; ++i)
				for (int imesh = i * meshnumi + 1; imesh <= (i + 1) * meshnumi; ++imesh)
					for (int j = 0; j < row; ++j)
						for (int jmesh = j * meshnumj + 1; jmesh <= (j + 1) * meshnumj; ++jmesh)
							Xs[hgroup][lgroup][imesh][jmesh] = Data2[grid[i][j]][hgroup][lgroup];

	// removal X section [cm-1] [group][i mesh][j mesh]
	vector<vector<vector<double>>> Xr(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	vector<vector<vector<double>>> Xs_sum(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	for (int i = 0; i < g; ++i)
		for (int j = 1; j <= Ni; ++j) {
			for (int k = 1; k <= Nj; ++k) {
				for (int l = 0; l < g; ++l)
					Xs_sum[i][j][k] += Xs[l][i][j][k];
				Xr[i][j][k] = Xa[i][j][k] + Xs_sum[i][j][k];
			}
		}


	// initial effective k
	double keff_old = 1.0f;
	double keff_new = 0.0f;

	// flux vector define [group][i mesh][j mesh]
	vector<vector<vector<double>>> phi_old(g, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f)));
	vector<vector<vector<double>>> phi_new(g, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f)));
	for (int i = 0; i < g; ++i) {
		for (int j = 1; j <= Ni; ++j)
			for (int k = 1; k <= Nj; ++k) {
				phi_old[i][j][k] = 1.0f;
				phi_new[i][j][k] = 1.0f;
			}
	}

	// beta vector define [group][i mesh][j mesh]
	vector<vector<vector<double>>> betai(g, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f))); // i direction beta vector
	vector<vector<vector<double>>> betaj(g, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f))); // j direction beta vector
	for (int i = 0; i < g; ++i) {
		for (int k = 1; k <= Nj; ++k) {
			betai[i][0][k] = r[2] / 2; // i direction West boundary beta
			betai[i][Ni + 1][k] = r[3] / 2; // i direction East boundary beta
		}
		for (int j = 1; j <= Ni; ++j) {
			betaj[i][j][0] = r[0] / 2; // j direction North boundary beta
			betaj[i][j][Nj + 1] = r[1] / 2; // j direction South boundary beta
			for (int k = 1; k <= Nj; ++k) {
				betai[i][j][k] = D[i][j][k] / meshsizei; // i direction beta
				betaj[i][j][k] = D[i][j][k] / meshsizej; // j direction beta
			}
		}
	}

	// D tilder vector define [group][i mesh][j mesh]
	vector<vector<vector<double>>> D_tilderN(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	vector<vector<vector<double>>> D_tilderS(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	vector<vector<vector<double>>> D_tilderW(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	vector<vector<vector<double>>> D_tilderE(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	for (int i = 0; i < g; ++i)
		for (int j = 1; j <= Ni; ++j)
			for (int k = 1; k <= Nj; ++k) {
				D_tilderN[i][j][k] = 2 * betaj[i][j][k] * betaj[i][j][k - 1] / (betaj[i][j][k] + betaj[i][j][k - 1]);
				D_tilderS[i][j][k] = 2 * betaj[i][j][k] * betaj[i][j][k + 1] / (betaj[i][j][k] + betaj[i][j][k + 1]);
				D_tilderW[i][j][k] = 2 * betai[i][j][k] * betai[i][j - 1][k] / (betai[i][j][k] + betai[i][j - 1][k]);
				D_tilderE[i][j][k] = 2 * betai[i][j][k] * betai[i][j + 1][k] / (betai[i][j][k] + betai[i][j + 1][k]);
			}

	// Matrix M define diagonal, bidiagonal [group][i mesh][j mesh]
	vector<vector<vector<double>>> diagonal(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	vector<vector<vector<double>>> bidiagonalN(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	vector<vector<vector<double>>> bidiagonalS(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	vector<vector<vector<double>>> bidiagonalW(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	vector<vector<vector<double>>> bidiagonalE(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
	for (int i = 0; i < g; ++i)
		for (int j = 1; j <= Ni; ++j)
			for (int k = 1; k <= Nj; ++k) {
				diagonal[i][j][k] = D_tilderE[i][j][k] / meshsizei + D_tilderW[i][j][k] / meshsizei + D_tilderN[i][j][k] / meshsizej + D_tilderS[i][j][k] / meshsizej + Xr[i][j][k];
				bidiagonalN[i][j][k] = -D_tilderN[i][j][k] / meshsizej;
				bidiagonalS[i][j][k] = -D_tilderS[i][j][k] / meshsizej;
				bidiagonalW[i][j][k] = -D_tilderW[i][j][k] / meshsizei;
				bidiagonalE[i][j][k] = -D_tilderE[i][j][k] / meshsizei;
			}


	// initial error
	double errk = 1.0f; // keff error
	double errp = 1.0f; // flux error
	double errp1 = 1.0f;
	double errp2 = 1.0f;

	// outer iteration counts
	int outiter = 1;

	// outer iteration
	while (errk > crik || errp > crip) {
		// Source vector define [high E group][low E group][i mesh][j mesh]
		vector<vector<vector<double>>> source(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
		vector<vector<vector<double>>> fission(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
		vector<vector<vector<double>>> scattering(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
		for (int i = 0; i < g; ++i) {
			for (int j = 1; j <= Ni; ++j)
				for (int k = 1; k <= Nj; ++k) {
					for (int l = 0; l < g; ++l) {
						fission[i][j][k] += 1 / keff_old * Xf[l][j][k] * phi_old[l][j][k] * xg[i];
						scattering[i][j][k] += Xs[i][l][j][k] * phi_old[l][j][k];
					}
					source[i][j][k] = fission[i][j][k] + scattering[i][j][k];
				}
		}

		// flux iteration (Gauss Seidel Method)
		for (int counts = 1; counts <= inneriter; ++counts)
			for (int i = 0; i < g; ++i)
				for (int j = 1; j <= Ni; ++j)
					for (int k = 1; k <= Nj; ++k)
						phi_new[i][j][k] = (source[i][j][k] - bidiagonalE[i][j][k] * phi_new[i][j + 1][k] - bidiagonalW[i][j][k] * phi_new[i][j - 1][k] - bidiagonalN[i][j][k] * phi_new[i][j][k - 1] - bidiagonalS[i][j][k] * phi_new[i][j][k + 1]) / diagonal[i][j][k];

		// keff calculation
		double numerator = 0.0f;
		double denominator = 0.0f;
		for (int j = 1; j <= Ni; ++j)
			for (int k = 1; k <= Nj; ++k)
				for (int i = 0; i < g; ++i) {
					numerator += Xf[i][j][k] * phi_new[i][j][k];
					denominator += Xf[i][j][k] * phi_old[i][j][k];
				}
		keff_new = keff_old * (numerator / denominator);

		// keff error calculation
		errk = abs((keff_new - keff_old) / keff_old);


		// flux error calculation
		double errp1, errp2;
		for (int i = 0; i < g; ++i)
			for (int j = 1; j < Ni; ++j)
				for (int k = 1; k <= Nj; ++k) {
					errp1 = abs((phi_new[i][j][k] - phi_old[i][j][k]) / phi_old[i][j][k]);
					errp2 = abs((phi_new[i][j + 1][k] - phi_old[i][j + 1][k]) / phi_old[i][j + 1][k]);
					errp = max(errp1, errp2);
				}

		// reset flux
		for (int i = 0; i < g; ++i)
			for (int j = 1; j <= Ni; ++j)
				for (int k = 1; k <= Nj; ++k)
					phi_old[i][j][k] = phi_new[i][j][k];

		// reset keff
		keff_old = keff_new;

		// outer iteration
		++outiter;
	}

	// check end time
	double end_time = clock();

	// duration calculation
	double duration = (end_time - start_time) * 1000 / CLOCKS_PER_SEC;

	// show result
	cout << "Mesh zone type [row X column] : " << Nj << " X " << Ni << "\n"
		<< "keff : " << keff_new << ", error : " << min(errk, errp) << "\n"
		<< "duration : " << duration << "ms" << ", iteration counts : " << outiter << "\n\n";

		// write
		ofstream fileoutput;
		fileoutput.open("output.txt");
		// write flux distribution for each group
		for (int i = 0; i < g; ++i) {
			fileoutput << "group " << i + 1 << "\n";
			for (int j = 1; j <= Ni; ++j) {
				for (int k = 1; k <= Nj; ++k) {
					fileoutput << phi_new[i][j][k] << "\t";
				}
				fileoutput << "\n";
			}
			fileoutput << "\n\n";
		}

}
