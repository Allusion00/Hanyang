#include "define.h"

// 필요한 크기의 vector를 생성후 0.0으로 채운 뒤, input.file에서 정보 읽은 후 저장하는 방식

vector<vector<vector<double>>> Data1; // D, Xs(Fission), Xs(Absorption)
vector<vector<vector<double>>> Data2; // Xs(Only Down Scattering
string line;
int main() {
	ifstream fileinput;
	fileinput.open("input.txt");
	getline(fileinput, line);
	int region = 0; // # of Region
	fileinput >> region;
	cout << "# of Region is : " << region << "\n"; // 출력
	int g = 0; // # of Energy group
	fileinput >> g;
	cout << "# of group is : " << g << "\n"; // 출력
	int n = 3; // # of types of Data1

	// Read Diffision coefficient, fission X-section, absorption X-section
	Data1 = vector<vector<vector<double>>>(region, vector<vector<double>>(g, vector<double>(n, 0.0f)));
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
	
	// Read scattering X-section (no up scattering)
	Data2 = vector<vector<vector<double>>>(region, vector<vector<double>>(g, vector<double>(g, 0.0f)));
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
	
	int NC = 10;
	int NR = 5;

	// Diffusion coefficient [cm]
	vector<vector<double>> D(g, vector<double>(NC + NR + 1, 0.0f));
	// fission X section [cm-1]
	vector<vector<double>> Xf(g, vector<double>(NC + NR + 1, 0.0f));
	// absorption X section [cm-1]
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
		cout << "group " << i+1 << " 의 D : ";
		for (int j = 0; j <= NC + NR; ++j) {
			cout << D[0][j] << " ";
		}
		cout << "\ngroup " << i + 1 << " 의 Xf : ";
		for (int j = 0; j <= NC + NR; ++j) {
			cout << Xf[0][j] << " ";
		}
		cout << "\ngroup " << i + 1 << " 의 Xa : ";
		for (int j = 0; j <= NC + NR; ++j) {
			cout << Xa[0][j] << " ";
		}
		cout << "\n";
	}
	
	// scattering X section [cm-1]
	vector<vector<vector<vector<double>>>> Xs(region, vector<vector<vector<double>>>(g, vector<vector<double>>(g, vector<double>(NC + NR + 1, 0.0f))));
	for (int i = 0; i < region; ++i) {
		for (int j = 0; j < g; ++j) {
			for (int k = 0; k < g; ++k) {
				for (int f = 0; f < NC + NR + 1; ++f) {
					cout << Xs[i][j][k][f] << " "; // 출력
				}
				cout << "\n"; // 출력
			}
			cout << "\n"; // 출력
		}
		cout << "\n"; // 출력
	}

}
