#include "define.h"
class Data {
private:
	ifstream inFP;
	ofstream ouFP;
public:

	double Me;
	int region = 0;
	int Group = 0;
	double s = 0;
	int c = 0;
	int a = 0, b = 0;
	int i = 0, j = 0;
	int p = 0, q = 0;
	string line;
	char Trash[50];
	void read(double& Size, vector<vector<int>>& Gm, vector<vector<int>>& Mesh, vector<int>& Boun, vector<vector<vector<double>>>& DFA, vector<vector<vector<double>>>& Xs, vector<double>& Xi, vector<double>& Conv) {
		inFP.open("input.file");
		getline(inFP, line);
		inFP >> region;
		inFP >> Group;
		inFP >> s;

		Size = s;
		getline(inFP, line);
		getline(inFP, line);
		inFP >> a;
		inFP >> b;
		for (i = 0; i < a; i++) {
			Gm.push_back(vector<int>{});
		}

		for (i = 0; i < a; i++) {
			for (j = 0; j < b; j++) {
				inFP >> Me;
				Gm[i].push_back(Me);
			}
		} //geometry

		for (i = 0; i < a; i++) {
			for (j = 0; j < b; j++) {
				cout << Gm[i][j] << " ";
			}
			cout << endl;
		} //geometry



		getline(inFP, line);
		getline(inFP, line);

		int* H_num_Data = nullptr;
		H_num_Data = new int[2];
		H_num_Data[0] = b;
		H_num_Data[1] = a;


		for (i = 0; i < 2; i++) {
			Mesh.push_back(vector<int>{});
		}
		for (i = 0; i < 2; i++) {
			for (j = 0; j < H_num_Data[i]; j++) {
				inFP >> Me;
				cout << Me << " ";
				Mesh[i].push_back(Me);
			}
		}
		getline(inFP, line);
		getline(inFP, line);

		for (i = 0; i < 4; i++) {
			inFP >> Me;

			Boun.push_back(Me);
		}
		cout << endl;
		for (i = 0; i < 4; i++) {
			cout << Boun[i] << endl;
		}
		getline(inFP, line);
		getline(inFP, line);

		for (i = 0; i < region; i++) {
			DFA.push_back(vector<vector<double>>(Group - 1)); //???
		}
		for (i = 0; i < region; i++) {
			DFA[i].push_back(vector<double>{});   //DFA data
		}

		for (i = 0; i < region; i++) {
			inFP >> Trash;
			for (j = 0; j < Group; j++) {
				for (p = 0; p < 3; p++) { //
					inFP >> Me;
					DFA[i][j].push_back(Me);
				}
			}
		}

		for (i = 0; i < DFA.size(); i++) {
			for (j = 0; j < DFA[0].size(); j++) {
				for (p = 0; p < DFA[0][0].size(); p++) {
					cout << DFA[i][j][p] << " ";
				}
				cout << endl;
			}
		}


		getline(inFP, line);
		getline(inFP, line);
		for (i = 0; i < region; i++) {
			Xs.push_back(vector<vector<double>>(Group - 1)); //???
		}
		for (i = 0; i < region; i++) {
			Xs[i].push_back(vector<double>{});   //Scattering data
		}
		for (i = 0; i < region; i++) {
			inFP >> Trash;
			for (j = 0; j < Group; j++) {
				for (a = 0; a < Group; a++) { //
					inFP >> Me;
					Xs[i][j].push_back(Me);
				}
			}
		}
		for (i = 0; i < Xs.size(); i++) {
			for (j = 0; j < Xs[0].size(); j++) {
				for (p = 0; p < Xs[0][0].size(); p++) {
					cout << Xs[i][j][p] << "  ";
				}
				cout << endl;
			}
		}

		getline(inFP, line);
		getline(inFP, line);

		for (i = 0; i < Group; i++) {
			inFP >> Me;
			cout << Me << endl;
			Xi.push_back(Me);
		}


		getline(inFP, line);
		getline(inFP, line);
		for (i = 0; i < 3; i++) {
			inFP >> Me;
			cout << Me << endl;
			Conv.push_back(Me);
		}

		inFP.close();
	}
	void Write(double*** A, int G, int Nt_V, int Nt_H, double keff, double e_keff, double e_source, vector<vector<double>>error) {
		ouFP.open("outout.file");
		ouFP << "K : " << keff << ", Relative error(keff, source) : " << e_keff << ", " << e_source << endl << endl;



		ouFP << "|| Keff Error Distribution || " << endl;
		for (i = 0; i < error[0].size(); i++) {
			ouFP << error[0][i] << " ";
		}

		ouFP << endl << endl << endl;

		ouFP << "|| Source Error Distribution || " << endl;
		for (i = 0; i < error[1].size(); i++) {
			ouFP << error[1][i] << " ";
		}


		for (Group = 0; Group < G; Group++) {
			ouFP << endl << endl;
			ouFP << "||" << Group + 1 << " Group Flux Distribution || " << endl;
			for (i = 1; i < Nt_V - 1; i++) {
				for (j = 1; j < Nt_H - 1; j++) {
					ouFP << A[Group][i][j] << " ";
				}
				ouFP << endl;
			}
		}

		ouFP.close();
	}
};
