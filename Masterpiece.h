class For2D {
public:
	int i = 0, j = 0, p = 0, q = 0;
	int k = 0, l = 0, a = 0, b = 0;
	double convcrit_keff = 0.0;
	double convcrit_source = 0.0;
	double error_keff = 0.0;
	double error_source_max = 0.0;
	int Iter = 0;
	int* hv_Data = nullptr;
	int* hh_Data = nullptr;
	vector<vector<int>>Chopping_Data;
	double** h_Data = nullptr;
	double*** DV = nullptr, *** DH = nullptr;
	double*** DV1 = nullptr, *** DH1 = nullptr;
	double*** DV2 = nullptr, *** DH2 = nullptr;  // 계산을 최적화 하기 위해서 만듬
	double*** D_sum = nullptr;
	double* Xi = nullptr;
	double*** x_old = nullptr, *** x_new = nullptr;
	double*** DFA = nullptr;
	double*** Xsa = nullptr;
	double** Beta = nullptr;
	int** Geometry = nullptr;
	vector<vector<int>> Geometry_vec;
	int size_h = 0;
	int size_v = 0;
	int num = 0;
	int G_num = 0;
	int Group = 0;
	int Region = 0;
	For2D(double size, vector<vector<int>> Geo, vector<vector<int>> Mesh, vector<int> Boundary, vector<vector<vector<double>>> DFA_, vector<vector<vector<double>>> Xsa_, vector<double>Xi_data, vector<double> Convergence) {
		hh_Data = new int[Mesh[0].size()];
		hv_Data = new int[Mesh[1].size()];

		for (i = 0; i < Mesh[0].size(); i++) {
			hh_Data[i] = Mesh[0][i];
		}
		for (i = 0; i < Mesh[0].size(); i++) {
			size_h += hh_Data[i];
		}
		for (i = 0; i < Mesh[1].size(); i++) {
			hv_Data[i] = Mesh[1][i];
		}
		for (i = 0; i < Mesh[1].size(); i++) {
			size_v += hv_Data[i];
		}
		size_h += 2;
		size_v += 2; // Boundary를 위함


		for (i = 0; i < 2; i++) {
			Chopping_Data.push_back(vector<int>{});
		}

		for (i = 0; i < Mesh[0].size(); i++) {
			for (j = 0; j < hh_Data[i]; j++) {
				Chopping_Data[0].push_back(hh_Data[i]);
			}
		}
		for (i = 0; i < Mesh[1].size(); i++) {
			for (j = 0; j < hv_Data[i]; j++) {
				Chopping_Data[1].push_back(hv_Data[i]);
			}
		}


		//for (j = 0; j < Chopping_Data[0].size(); j++) {
		//	cout << Chopping_Data[0][j] << "  ";  //수평
		//}
		//cout << endl << endl << endl << endl;
		//for (j = 0; j < Chopping_Data[1].size(); j++) {
		//	cout << Chopping_Data[1][j] << "  ";  //수직
		//}

		h_Data = new double* [2];
		h_Data[0] = new double[Chopping_Data[0].size()];  // 수평방향
		h_Data[1] = new double[Chopping_Data[1].size()]; //수직방향



		for (i = 0; i < Chopping_Data[0].size(); i++) {
			h_Data[0][i] = size / Chopping_Data[0][i]; //수평방향
		}
		for (i = 0; i < Chopping_Data[1].size(); i++) {
			h_Data[1][i] = size / Chopping_Data[1][i]; //수직방향
		}



		//cout << endl << endl << endl;
		//for (i = 0; i < Chopping_Data[0].size(); i++) {
		//	cout << h_Data[0][i] << " ";
		//}
		//cout << endl << endl << endl;
		//for (i = 0; i < Chopping_Data[1].size(); i++) {
		//	cout << h_Data[1][i] << " ";
		//}
		//
		//cout << endl << endl << endl << endl;



		G_num = DFA_[0].size();

		Geometry = new int* [size_v];
		for (i = 0; i < size_v; i++) {
			Geometry[i] = new int[size_h];
		}

		for (j = 0; j < size_h; j++) {
			Geometry[0][j] = Boundary[0];
		}
		for (j = 0; j < size_h; j++) {
			Geometry[size_v - 1][j] = Boundary[1];
		}
		for (i = 0; i < size_v; i++) {
			Geometry[i][0] = Boundary[2];
		}
		for (i = 0; i < size_v; i++) {
			Geometry[i][size_h - 1] = Boundary[3];   //boundary
		}


		for (i = 0; i < size_v - 2; i++) {
			Geometry_vec.push_back(vector<int>{});
		}



		k = 0;
		for (p = 0; p < Mesh[1].size(); p++) { // 코스 메쉬 개수 수직방향
			for (i = 0; i < hv_Data[p]; i++) {
				for (q = 0; q < Mesh[0].size(); q++) { //코스메쉬 개수 수평방향
					for (j = 0; j < hh_Data[q]; j++) {
						Geometry_vec[k + i].push_back(Geo[p][q]);
					}
				}
			}
			k += hv_Data[p];
		}


		for (i = 0; i < size_v - 2; i++) {
			for (j = 0; j < size_h - 2; j++) {
				Geometry[i + 1][j + 1] = Geometry_vec[i][j];
			}

		}

		cout << endl << endl << endl << endl;
		for (i = 1; i < size_v - 1; i++) {
			for (j = 1; j < size_h - 1; j++) {
				cout << Geometry[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;




		DFA = new double** [20];   //boundary 를 넣기 위해서 10, 11, 12
		for (i = 0; i < 20; i++) {
			DFA[i] = new double* [G_num];
		}
		for (i = 0; i < 20; i++) {
			for (j = 0; j < G_num; j++) {
				DFA[i][j] = new double[DFA_[0][0].size()];
			}

		}
		for (i = 0; i < DFA_.size(); i++) {
			for (j = 0; j < G_num; j++) {
				for (p = 0; p < DFA_[0][0].size(); p++) {
					DFA[i][j][p] = DFA_[i][j][p];
				}
			}
		} // XS Data

		Xsa = new double** [Xsa_.size()];
		for (i = 0; i < Xsa_.size(); i++) {
			Xsa[i] = new double* [G_num];
		}
		for (i = 0; i < Xsa_.size(); i++) {
			for (j = 0; j < G_num; j++) {
				Xsa[i][j] = new double[G_num];
			}

		}
		for (i = 0; i < Xsa_.size(); i++) {
			for (j = 0; j < G_num; j++) {
				for (p = 0; p < G_num; p++) {
					Xsa[i][j][p] = Xsa_[i][j][p];
				}
			}
		} // Sacttering


		for (j = 0; j < G_num; j++) {
			DFA[10][j][0] = 1.0e-10;
			DFA[11][j][0] = 1.0e+10;
			DFA[12][j][0] = 0.25;

		}//  Boundary


		DV = new double** [G_num]; //group 개수
		for (i = 0; i < G_num; i++) {
			DV[i] = new double* [size_v - 1]; //전체 갯수 +2 -1
		}
		for (i = 0; i < G_num; i++) {
			for (j = 0; j < size_v - 1; j++) {
				DV[i][j] = new double[size_h - 2];
			}
		}

		DV1 = new double** [G_num]; //group 개수
		for (i = 0; i < G_num; i++) {
			DV1[i] = new double* [size_v - 2];
		}
		for (i = 0; i < G_num; i++) {
			for (j = 0; j < size_v - 1; j++) {
				DV1[i][j] = new double[size_h - 2];
			}
		}

		DV2 = new double** [G_num]; //group 개수
		for (i = 0; i < G_num; i++) {
			DV2[i] = new double* [size_v - 2];
		}
		for (i = 0; i < G_num; i++) {
			for (j = 0; j < size_v - 1; j++) {
				DV2[i][j] = new double[size_h - 2];
			}
		}


		DH = new double** [G_num];
		for (i = 0; i < G_num; i++) {
			DH[i] = new double* [size_v - 2];  //전체 갯수 +2 -1 
		}
		for (i = 0; i < G_num; i++) {
			for (j = 0; j < size_v - 2; j++) {
				DH[i][j] = new double[size_h - 1];
			}
		}

		DH1 = new double** [G_num];
		for (i = 0; i < G_num; i++) {
			DH1[i] = new double* [size_v - 2];
		}
		for (i = 0; i < G_num; i++) {
			for (j = 0; j < size_v - 2; j++) {
				DH1[i][j] = new double[size_h - 2];
			}
		}

		DH2 = new double** [G_num];
		for (i = 0; i < G_num; i++) {
			DH2[i] = new double* [size_v - 2];
		}
		for (i = 0; i < G_num; i++) {
			for (j = 0; j < size_v - 2; j++) {
				DH2[i][j] = new double[size_h - 2];
			}
		}


		for (Group = 0; Group < G_num; Group++) {    //2
			for (i = 0; i < size_v - 1; i++) {        //51
				for (j = 0; j < size_h - 2; j++) {     //50
					if (Geometry[i][j + 1] == 10 || Geometry[i][j + 1] == 11 || Geometry[i][j + 1] == 12) {
						DV[Group][i][j] = 2 * (DFA[Geometry[i][j + 1]][Group][0] * DFA[Geometry[i + 1][j + 1]][Group][0] / h_Data[1][i]) / (DFA[Geometry[i][j + 1]][Group][0] + DFA[Geometry[i + 1][j + 1]][Group][0] / h_Data[1][i]);

					}
					else if (Geometry[i + 1][j + 1] == 10 || Geometry[i + 1][j + 1] == 11 || Geometry[i + 1][j + 1] == 12) {
						DV[Group][i][j] = 2 * (DFA[Geometry[i][j + 1]][Group][0] / h_Data[1][i - 1] * DFA[Geometry[i + 1][j + 1]][Group][0]) / (DFA[Geometry[i][j + 1]][Group][0] / h_Data[1][i - 1] + DFA[Geometry[i + 1][j + 1]][Group][0]);
					}
					else {
						DV[Group][i][j] = 2 * (DFA[Geometry[i][j + 1]][Group][0] / h_Data[1][i - 1] * DFA[Geometry[i + 1][j + 1]][Group][0] / h_Data[1][i]) / (DFA[Geometry[i][j + 1]][Group][0] / h_Data[1][i - 1] + DFA[Geometry[i + 1][j + 1]][Group][0] / h_Data[1][i]);
					}
				}
			}
		}


		for (Group = 0; Group < G_num; Group++) {       //2
			for (i = 0; i < size_v - 2; i++) {                //50 
				for (j = 0; j < size_h - 1; j++) {     //51
					if (Geometry[i + 1][j] == 10 || Geometry[i + 1][j] == 11 || Geometry[i + 1][j] == 12) {
						DH[Group][i][j] = 2 * (DFA[Geometry[i + 1][j]][Group][0] * DFA[Geometry[i + 1][j + 1]][Group][0] / h_Data[0][j]) / (DFA[Geometry[i + 1][j]][Group][0] + DFA[Geometry[i + 1][j + 1]][Group][0] / h_Data[0][j]);
					}
					else if (Geometry[i + 1][j + 1] == 10 || Geometry[i + 1][j + 1] == 11 || Geometry[i + 1][j + 1] == 12) {
						DH[Group][i][j] = 2 * (DFA[Geometry[i + 1][j]][Group][0] / h_Data[0][j - 1] * DFA[Geometry[i + 1][j + 1]][Group][0]) / (DFA[Geometry[i + 1][j]][Group][0] / h_Data[0][j - 1] + DFA[Geometry[i + 1][j + 1]][Group][0]);
					}
					else {
						DH[Group][i][j] = 2 * (DFA[Geometry[i + 1][j]][Group][0] / h_Data[0][j - 1] * DFA[Geometry[i + 1][j + 1]][Group][0] / h_Data[0][j]) / (DFA[Geometry[i + 1][j]][Group][0] / h_Data[0][j - 1] + DFA[Geometry[i + 1][j + 1]][Group][0] / h_Data[0][j]);
					}
				}
			}
		}


		D_sum = new double** [G_num];
		for (i = 0; i < G_num; i++) { //boundary 제외한 수
			D_sum[i] = new double* [size_v - 2];
		}
		for (i = 0; i < G_num; i++) {
			for (j = 0; j < size_v - 2; j++) {
				D_sum[i][j] = new double[size_h - 2];
			}
		}

		for (Group = 0; Group < G_num; Group++) {
			for (i = 0; i < size_v - 2; i++) {
				for (j = 0; j < size_h - 2; j++) {

					DV1[Group][i][j] = DV[Group][i][j] / h_Data[1][i];
					DV2[Group][i][j] = DV[Group][i + 1][j] / h_Data[1][i];
					DH1[Group][i][j] = DH[Group][i][j] / h_Data[0][j];
					DH2[Group][i][j] = DH[Group][i][j + 1] / h_Data[0][j];
				}
			}
		}

		double Xr = 0.0;
		for (Group = 0; Group < G_num; Group++) {
			for (i = 0; i < size_v - 2; i++) {
				for (j = 0; j < size_h - 2; j++) {
					for (p = 0; p < G_num; p++) { // 산란단면적 + 흡수단면적 = 그룹수랑 같음
						if (p != Group) {
							Xr += Xsa[Geometry[i + 1][j + 1]][Group][p];
						}
					}
					Xr += DFA[Geometry[i + 1][j + 1]][Group][2];
					D_sum[Group][i][j] = DV[Group][i][j] / h_Data[1][i] + DV[Group][i + 1][j] / h_Data[1][i] + DH[Group][i][j] / h_Data[0][j] + DH[Group][i][j + 1] / h_Data[0][j] + Xr;
					Xr = 0.0; //초기화
				}
			}
		}


		//cout << fixed;
		//cout.precision(3);

		//for (Group = 0;  Group < G_num; Group++) {
		//	for (i = 0; i < size_v - 2; i++) {
		//		for (j = 0; j < size_h - 1; j++) {
		//			cout << DH1[Group][i][j] << " ";
		//		}
		//		cout << endl;
		//	}
		//	cout << endl<<endl<<endl;
		//}





		Xi = new double[Group];
		for (i = 0; i < Group; i++) {
			Xi[i] = Xi_data[i];
		}

		Iter = Convergence[0];
		convcrit_keff = Convergence[1];
		convcrit_source = Convergence[2]; //Convergence Criteria
	}

	void Gauss_Seidal() {

		x_new = new double** [G_num];
		x_old = new double** [G_num];
		for (Group = 0; Group < G_num; Group++) {
			x_new[Group] = new double* [size_v];
			x_old[Group] = new double* [size_v];
		}
		for (Group = 0; Group < G_num; Group++) {
			for (i = 0; i < size_v; i++) {
				x_new[Group][i] = new double[size_h];
				x_old[Group][i] = new double[size_h];
			}
		}

		for (Group = 0; Group < G_num; Group++) {
			for (i = 0; i < size_v; i++) {
				for (j = 0; j < size_h; j++) {
					x_new[Group][i][j] = 0.0;
					x_old[Group][i][j] = 0.0; //바운더리 0.0
				}

			}
		}
		for (Group = 0; Group < G_num; Group++) {
			for (i = 1; i < size_v - 1; i++) {
				for (j = 1; j < size_h - 1; j++) {
					x_new[Group][i][j] = 1.0;
					x_old[Group][i][j] = 1.0;  //초기값 1.0
				}
			}
		}

		double* Fsum_new = nullptr;
		double* Fsum_old = nullptr;

		Fsum_new = new double[G_num];
		Fsum_old = new double[G_num];
		for (i = 0; i < G_num; i++) {
			Fsum_new[i] = 0.0;
			Fsum_old[i] = 0.0;
		}

		int k = 0;
		for (k = 0; k < G_num; k++) {
			for (i = 0; i < size_v - 2; i++) {
				for (j = 0; j < size_h - 2; j++) {
					Fsum_old[k] += x_old[k][i + 1][j + 1] * DFA[Geometry[i + 1][j + 1]][k][1];

				}
			}
		}
		double All_Fsum_old = 0.0;
		double All_Fsum_new = 0.0;
		for (k = 0; k < G_num; k++) {
			All_Fsum_old += Fsum_old[k];
		}



		vector<vector<double>>error;
		for (i = 0; i < 2; i++) {
			error.push_back(vector<double>{});
		}

		double** error_s_old = nullptr;
		error_s_old = new double* [size_v - 2];
		for (i = 0; i < size_v - 2; i++) {
			error_s_old[i] = new double[size_h - 2];
		}
		double** error_s_new = nullptr;
		error_s_new = new double* [size_v - 2];
		for (i = 0; i < size_v - 2; i++) {
			error_s_new[i] = new double[size_h - 2];
		}


		for (i = 0; i < size_v - 2; i++) {
			for (j = 0; j < size_h - 2; j++) {
				error_s_old[i][j] = 0.0; //초기화
				for (k = 0; k < G_num; k++) {
					error_s_old[i][j] += DFA[Geometry[i + 1][j + 1]][k][1] * x_old[k][i + 1][j + 1];
				}
			}
		}


		error_keff = 1.0;
		error_source_max = 1.0;
		double error_source = 0.0;
		double Fsource = 0.0;
		double Ssource = 0.0;
		double k_new = 1.0;
		double k_old = 1.0;
		int Iteration = 0;
		int counts = 0;
		int L = 1;
		double** source_sum = nullptr;
		source_sum = new double* [size_v - 2];
		for (i = 0; i < size_v - 2; i++) {
			source_sum[i] = new double[size_h - 2];
		}

		clock_t start, end;
		double result;
		start = clock();
		while (error_keff > convcrit_keff || error_source_max > convcrit_source) {
			Iteration++;
			for (Group = 0; Group < G_num; Group++) {

				for (i = 0; i < size_v - 2; i++) {
					for (j = 0; j < size_h - 2; j++) {
						Fsource = 0.0; //초기화
						for (k = 0; k < G_num; k++) {
							Fsource += DFA[Geometry[i + 1][j + 1]][k][1] * x_old[k][i + 1][j + 1] * Xi[Group];
						}
						Fsource = Fsource / k_old;

						Ssource = 0.0; //초기화
						for (k = 0; k < G_num; k++) {
							if (k != Group) {
								Ssource += Xsa[Geometry[i + 1][j + 1]][k][Group] * x_old[k][i + 1][j + 1];
							}
						}
						source_sum[i][j] = Fsource + Ssource;
						//inner Iteration 전에 소스계산
					}
				}
				for (counts = 0; counts < Iter; counts++) {
					for (i = 0; i < size_v - 2; i++) {
						for (j = 0; j < size_h - 2; j++) {
							x_new[Group][i + 1][j + 1] = (source_sum[i][j]
								+ DV1[Group][i][j] * x_new[Group][i][j + 1]
								+ DV2[Group][i][j] * x_new[Group][i + 2][j + 1]
								+ DH1[Group][i][j] * x_new[Group][i + 1][j]
								+ DH2[Group][i][j] * x_new[Group][i + 1][j + 2]) / D_sum[Group][i][j];
						}
					}
				}

				for (i = 0; i < size_v; i++) {
					for (j = 0; j < size_h; j++) {
						x_old[Group][i][j] = x_new[Group][i][j];
					}
				}

			} // Outer Iteration


			for (i = 0; i < size_v - 2; i++) {
				for (j = 0; j < size_h - 2; j++) {
					error_s_new[i][j] = 0.0; //초기화
					for (k = 0; k < G_num; k++) {
						error_s_new[i][j] += DFA[Geometry[i + 1][j + 1]][k][1] * x_new[k][i + 1][j + 1];
					}
				}
			}  // For source error


			error_source_max = 0.0;
			for (i = 0; i < size_v - 2; i++) {
				for (j = 0; j < size_h - 2; j++) {
					error_source = abs((error_s_new[i][j] - error_s_old[i][j]) / error_s_old[i][j]);
					if (error_source > error_source_max) {
						error_source_max = error_source;
					}
				}
			}

			for (i = 0; i < size_v - 2; i++) {
				for (j = 0; j < size_h - 2; j++) {
					error_s_old[i][j] = error_s_new[i][j];

				}
			}  // For source error


			for (i = 0; i < G_num; i++) {
				Fsum_new[i] = 0.0;  // 값 초기화
			}

			for (k = 0; k < G_num; k++) {
				for (i = 0; i < size_v - 2; i++) {
					for (j = 0; j < size_h - 2; j++) {
						Fsum_new[k] += x_new[k][i + 1][j + 1] * DFA[Geometry[i + 1][j + 1]][k][1] * h_Data[1][i] * h_Data[0][j];
					}
				}
			}

			All_Fsum_new = 0.0; //초기화
			for (k = 0; k < G_num; k++) {
				All_Fsum_new += Fsum_new[k];
			}

			k_new = k_new * (All_Fsum_new / All_Fsum_old);
			error_keff = abs((k_new - k_old) / k_old); //Error 계산
			k_old = k_new;
			All_Fsum_old = All_Fsum_new;


			error[0].push_back(error_keff);
			error[1].push_back(error_source_max); //For Error Distribution





			if (L == 1000) {
				cout << "Iteration : " << Iteration << "  K :  " << k_new << endl;
				cout << "Relative error(keff, Soruce) : " << error_keff << " , " << error_source_max << endl << endl;
				L = 1;
			}
			else if (L != 1000) {
				L++;
			}



		}




		/*
		for (Group = 0; Group < G_num; Group++) {
			for (i = 0; i < size_v; i++) {
				for (j = 0; j < size_h; j++) {
					cout << x_old[Group][i][j] << "  ";
				}
				cout << endl;
			}
		}*/
		end = clock();
		result = end - start;
		cout << endl << endl;
		cout << "Iteration : " << Iteration << "  K :  " << k_new << endl;
		cout << "Relative error(keff, Soruce) : " << error_keff << " , " << error_source_max << endl << endl;

		cout << fixed;
		cout.precision(5);
		cout << " 수행시간 : " << result << "ms";


		/*for (Group = 0; Group < G_num; Group++) {
			for (i = 0; i < size_v; i++) {
				for (j = 0; j < size_h; j++) {
					cout << x_new[Group][i][j] << " ";
				}
				cout << end;
			}
		}*/


		double Normal = 200 / (3.2 * 1.0e-11 * All_Fsum_new);



		for (Group = 0; Group < G_num; Group++) {
			for (i = 0; i < size_v; i++) {
				for (j = 0; j < size_h; j++) {
					x_new[Group][i][j] = Normal * x_new[Group][i][j];
				}
			}
		}

		for (i = 0; i < G_num; i++) {
			Fsum_new[i] = 0.0;  // 값 초기화
		}
		for (k = 0; k < G_num; k++) {
			for (i = 0; i < size_v - 2; i++) {
				for (j = 0; j < size_h - 2; j++) {
					Fsum_new[k] += x_new[k][i + 1][j + 1];

				}
			}
		}

		cout << endl << endl << endl << endl;
		cout << fixed;
		cout.precision(10);
		for (i = 0; i < G_num; i++) {
			cout << Fsum_new[i] << " ";
		}


		Data Flux;
		Flux.Write(x_new, G_num, size_v, size_h, k_new, error_keff, error_source_max, error);

	}

};
