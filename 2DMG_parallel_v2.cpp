//#include <iostream>		// I/O���꿡 �ʿ��� �⺻���� ���
//#include <fstream>		// ���� ����� �� �� �ְ� ����� ���̺귯��
//#include <string>
//#include <vector>
//#include <omp.h>
//using namespace std;	// �̸��� �浹�� �����ְ�, ���α׷��� �������� �����ִ� �ڵ�
//// Visual 2008���ʹ� iostream.h�� �⺻���� �������� �ʱ� ������ #include <iostream>�� �ϰų� cout�� ������ �� �ڵ带 ����Ѵ�.
//
//// �ʿ��� ũ���� vector�� ������ 0.0���� ä�� ��, input.file���� ���� ���� �� �����ϴ� ���
//// check start time
//double start_time = clock();
//
//int main() {
//    ifstream fileinput; // define ifstream as fileinput
//
//    // Select group by cin
//    int groupinput;
//    cout << "Select group between 1, 2, 7, 0 (C5G7), 22 (Kang) : ";
//    cin >> groupinput;
//    if (groupinput == 1)
//        fileinput.open("2D1G input.txt");
//    else if (groupinput == 2)
//        fileinput.open("2D2G input.txt");
//    else if (groupinput == 7)
//        fileinput.open("2D7G input.txt");
//    else if (groupinput == 0)
//        fileinput.open("C5G7 input.txt");
//    else if (groupinput == 22)
//        fileinput.open("Kang input.txt");
//    else {
//        cout << "Please select group between 1, 2, 7";
//        return 0;
//    }
//
//    string line; // define line as string character for getline
//
//    // Read CARD1 (# of Region, # of Energy group)
//    cout << "\n[CARD1]\n";
//    getline(fileinput, line);
//    int region, group; // # of Region, # of Energy group
//    fileinput >> region >> group;
//    cout << "# of Region : " << region << "\n" << "# of group : " << group << "\n\n"; // output for checking CARD1
//
//    // Read CARD2 (Grid size [row column], Geometry size [row column], Geometry)
//    cout << "[CARD2]\n";
//    getline(fileinput, line);
//    getline(fileinput, line);
//    double gridsizei, gridsizej;
//    fileinput >> gridsizei >> gridsizej; // Grid size [cm]
//    int rowsize, columnsize;
//    fileinput >> rowsize >> columnsize;
//    vector<vector<int>> grid(rowsize, vector<int>(columnsize, 0));
//    for (int row = 0; row < rowsize; ++row) {
//        getline(fileinput, line);
//        for (int column = 0; column < columnsize; ++column) {
//            fileinput >> grid[row][column];
//            cout << grid[row][column] << " ";
//        }
//        cout << "\n";
//    }
//
//    // Read CARD3 (# of fine meshes for each grid [horizontal vertical])
//    cout << "\n[CARD3]\n";
//    getline(fileinput, line);
//    getline(fileinput, line);
//    vector<int> horizontal(columnsize, 0);
//    vector<int> vertical(rowsize, 0);
//    int Ni = 0;
//    int Nj = 0;
//    for (int column = 0; column < columnsize; ++column) {
//        fileinput >> horizontal[column];
//    }
//    getline(fileinput, line);
//    for (int row = 0; row < rowsize; ++row) {
//        fileinput >> vertical[row];
//    }
//    vector<int> horizontalsum(columnsize + 1, 0);
//    vector<int> verticalsum(rowsize + 1, 0);
//    for (int column = 1; column <= columnsize; ++column)
//        for (int i = 0; i <= column - 1; ++i)
//            horizontalsum[column] += horizontal[i];
//    for (int row = 1; row <= rowsize; ++row)
//        for (int i = 0; i <= row - 1; ++i)
//            verticalsum[row] += vertical[i];
//    Ni = horizontalsum[columnsize];
//    Nj = verticalsum[rowsize];
//    vector<double> hi(Ni + 1, 0.0f);
//    vector<double> hj(Nj + 1, 0.0f);
//    for (int column = 0; column < columnsize; ++column)
//        for (int i = horizontalsum[column] + 1; i <= horizontalsum[column + 1]; ++i)
//            for (int row = 0; row < rowsize; ++row)
//                for (int j = verticalsum[row] + 1; j <= verticalsum[row + 1]; ++j) {
//                    hi[i] = gridsizei / horizontal[column];
//                    hj[j] = gridsizej / vertical[row];
//                }
//    cout << "Mesh zone type [row X column] : " << Nj << " X " << Ni << "\n";
//
//    // Read CARD4 (Boundary Condition[Left Right] : 10 = reflective, 11 = flux zero, 12 = vacuum)
//    cout << "\n[CARD4]\n";
//    getline(fileinput, line);
//    getline(fileinput, line);
//    vector<int> BC(4, 0);
//    for (int i = 0; i < 4; ++i)
//        fileinput >> BC[i];
//    vector<double> b(4, 0.0f);
//    for (int i = 0; i < 4; ++i) {
//        if (BC[i] == 10)
//            b[i] = 10E-10f; // reflective
//        else if (BC[i] == 11)
//            b[i] = 10E+10f; // flux zero
//        else if (BC[i] == 12)
//            b[i] = 0.5; // vacuum
//        else
//            return 0;
//    }
//    vector<string> cardinal{ "North", "South", "West", "East" };
//    // output for checking CARD2
//    for (int i = 0; i < 4; ++i)
//        cout << cardinal[i] << " Boundary Condition : " << BC[i] << " (" << b[i] << ")\n\n";
//
//    // Read CARD5 (Convergence criteria : # of inner iterations, relative error of keff, max relative error of flux)
//    cout << "[CARD5]\n";
//    getline(fileinput, line);
//    getline(fileinput, line);
//    int inneriter;
//    double crik, crip;
//    fileinput >> inneriter >> crik >> crip;
//    // output for checking CARD3
//    cout << "# of inner iterations : " << inneriter << "\nkeff error : " << crik << "\nsource group error : " << crip << "\n\n";
//
//    // Read CARD6 (D, Xs(Fission), Xs(Absorption)) [region][group][data type]
//    cout << "[CARD6]\n";
//    int data = 3; // # of types of CARD6
//    vector<vector<vector<double>>> CARD6(region, vector<vector<double>>(group, vector<double>(data, 0.0f)));
//    getline(fileinput, line);
//    getline(fileinput, line);
//    for (int r = 0; r < region; ++r) {
//        for (int g = 0; g < group; ++g) {
//            getline(fileinput, line);
//            for (int i = 0; i < data; ++i) {
//                fileinput >> CARD6[r][g][i];
//                cout << CARD6[r][g][i] << " "; // output for checking CARD6
//            }
//            cout << "\n"; // output for checking CARD6
//        }
//        getline(fileinput, line);
//        cout << "\n"; // output for checking CARD6
//    }
//
//
//    // Read CARD7 (Xs(Only Down Scattering)) [region][High E group][Low E group]
//    cout << "[CARD7]\n";
//    vector<vector<vector<double>>> CARD7(region, vector<vector<double>>(group, vector<double>(group, 0.0f)));
//    getline(fileinput, line);
//    for (int r = 0; r < region; ++r) {
//        for (int lg = 0; lg < group; ++lg) {
//            getline(fileinput, line);
//            for (int hg = 0; hg < group; ++hg) {
//                fileinput >> CARD7[r][hg][lg];
//                cout << CARD7[r][lg][hg] << " "; // output for checking CARD7
//            }
//            cout << "\n"; // output for checking CARD7
//        }
//        getline(fileinput, line);
//        cout << "\n"; // output for checking CARD7
//    }
//
//    // Read CARD8 (xg : fraction of fission neutrons)
//    cout << "[CARD8]\n";
//    // xg [group]
//    vector<double> xg(group, 0.0f);
//    getline(fileinput, line);
//    for (int g = 0; g < group; ++g) {
//        fileinput >> xg[g];
//        cout << "group " << g + 1 << " �� fission neutron fraction : " << xg[g] << "\n"; // output for checking CARD8
//    }
//    cout << "\n\n"; // output for checking CARD8
//
//    // Diffusion coefficient [cm] [group][i mesh][j mesh]
//    vector<vector<vector<double>>> D(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    // fission X section [cm-1] [group][i mesh][j mesh]
//    vector<vector<vector<double>>> Xf(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    // absorption X section [cm-1] [group][i mesh][j mesh]
//    vector<vector<vector<double>>> Xa(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    for (int column = 0; column < columnsize; ++column)
//        for (int i = horizontalsum[column] + 1; i <= horizontalsum[column + 1]; ++i)
//            for (int row = 0; row < rowsize; ++row)
//                for (int j = verticalsum[row] + 1; j <= verticalsum[row + 1]; ++j)
//                    for (int g = 0; g < group; ++g) {
//                        D[i][j][g] = CARD6[grid[row][column]][g][0];
//                        Xf[i][j][g] = CARD6[grid[row][column]][g][1];
//                        Xa[i][j][g] = CARD6[grid[row][column]][g][2];
//                    }
//
//    // scattering X section [cm-1] [i mesh][j mesh][high E group][low E group]
//    vector<vector<vector<vector<double>>>> Xs(Ni + 1, vector<vector<vector<double>>>(Nj + 1, vector<vector<double>>(group, vector<double>(group, 0.0f))));
//    for (int column = 0; column < columnsize; ++column)
//        for (int i = horizontalsum[column] + 1; i <= horizontalsum[column + 1]; ++i)
//            for (int row = 0; row < rowsize; ++row)
//                for (int j = verticalsum[row] + 1; j <= verticalsum[row + 1]; ++j)
//                    for (int hg = 0; hg < group; ++hg)
//                        for (int lg = 0; lg < group; ++lg)
//                            Xs[i][j][hg][lg] = CARD7[grid[row][column]][hg][lg];
//
//    // removal X section [cm-1] [i mesh][j mesh][group]
//    vector<vector<vector<double>>> Xr(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    for (int i = 1; i <= Ni; ++i)
//        for (int j = 1; j <= Nj; ++j)
//            for (int g1 = 0; g1 < group; ++g1) {
//                for (int g2 = 0; g2 < group; ++g2) {
//                    if (g1 == g2)
//                        continue;
//                    Xr[i][j][g1] += Xs[i][j][g2][g1];
//                }
//                Xr[i][j][g1] += Xa[i][j][g1];
//            }
//
//
//    // initial effective k
//    double keff_old = 1.0f;
//    double keff_new = 0.0f;
//
//    // flux vector define [i mesh][j mesh][group]
//    vector<vector<vector<double>>> phi_old(Ni + 2, vector<vector<double>>(Nj + 2, vector<double>(group, 0.0f)));
//    vector<vector<vector<double>>> phi_new(Ni + 2, vector<vector<double>>(Nj + 2, vector<double>(group, 0.0f)));
//    for (int i = 1; i <= Ni; ++i)
//        for (int j = 1; j <= Nj; ++j)
//            for (int g = 0; g < group; ++g) {
//                phi_old[i][j][g] = 1.0f;
//                phi_new[i][j][g] = 1.0f;
//            }
//
//    // beta vector define [i mesh][j mesh][group]
//    vector<vector<vector<double>>> betai(Ni + 2, vector<vector<double>>(Nj + 2, vector<double>(group, 0.0f))); // i direction beta vector
//    vector<vector<vector<double>>> betaj(Ni + 2, vector<vector<double>>(Nj + 2, vector<double>(group, 0.0f))); // j direction beta vector
//    for (int j = 1; j <= Nj; ++j)
//        for (int g = 0; g < group; ++g) {
//            betai[0][j][g] = b[2] / 2; // i direction West boundary beta
//            betai[Ni + 1][j][g] = b[3] / 2; // i direction East boundary beta
//        }
//    for (int i = 1; i <= Ni; ++i)
//        for (int g = 0; g < group; ++g) {
//            betaj[i][0][g] = b[0] / 2; // j direction North boundary beta
//            betaj[i][Nj + 1][g] = b[1] / 2; // j direction South boundary beta   
//        }
//    for (int i = 1; i <= Ni; ++i)
//        for (int j = 1; j <= Nj; ++j)
//            for (int g = 0; g < group; ++g) {
//                betai[i][j][g] = D[i][j][g] / hi[i]; // i direction beta
//                betaj[i][j][g] = D[i][j][g] / hj[j]; // j direction beta
//            }
//
//    // D tilder vector define [i mesh][j mesh][group]
//    vector<vector<vector<double>>> D_tilderi(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    vector<vector<vector<double>>> D_tilderj(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    for (int i = 0; i <= Ni; ++i)
//        for (int j = 0; j <= Nj; ++j)
//            for (int g = 0; g < group; ++g) {
//                D_tilderi[i][j][g] = 2 * betai[i][j][g] * betai[i + 1][j][g] / (betai[i][j][g] + betai[i + 1][j][g]);
//                D_tilderj[i][j][g] = 2 * betaj[i][j][g] * betaj[i][j + 1][g] / (betaj[i][j][g] + betaj[i][j + 1][g]);
//            }
//
//    // Matrix M define diagonal, bidiagonal [i mesh][j mesh][group]
//    vector<vector<vector<double>>> diagonal(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    vector<vector<vector<double>>> bidiagonali(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    vector<vector<vector<double>>> bidiagonalj(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//    for (int i = 1; i <= Ni; ++i)
//        for (int j = 1; j <= Nj; ++j)
//            for (int g = 0; g < group; ++g) {
//                diagonal[i][j][g] = D_tilderi[i][j][g] / hi[i] + D_tilderi[i - 1][j][g] / hi[i] + D_tilderj[i][j - 1][g] / hj[j] + D_tilderj[i][j][g] / hj[j] + Xr[i][j][g];
//                bidiagonali[i][j][g] = -D_tilderi[i][j][g] / hi[i];
//                bidiagonalj[i][j][g] = -D_tilderj[i][j][g] / hj[j];
//            }
//
//
//    // initial error
//    double errk = 1.0f; // keff error
//    double errp = 1.0f; // flux error
//    double errp1 = 1.0f;
//    double errp2 = 1.0f;
//
//    // outer iteration counts
//    int outiter = 1;
//
//    // outer iteration
//    while (errk > crik || errp > crip) {
//        // Source vector define [high E group][low E group][i mesh][j mesh]
//        vector<vector<vector<double>>> source(Ni + 1, vector<vector<double>>(Nj + 1, vector<double>(group, 0.0f)));
//        #pragma omp parallel for
//        for (int i = 1; i <= Ni; ++i)
//            for (int j = 1; j <= Nj; ++j)
//                for (int g1 = 0; g1 < group; ++g1) {
//                    for (int g2 = 0; g2 < group; ++g2) {
//                        source[i][j][g1] += 1 / keff_old * Xf[i][j][g2] * phi_old[i][j][g2] * xg[g1];
//                        if (g1 == g2)
//                            continue;
//                        source[i][j][g1] += Xs[i][j][g1][g2] * phi_old[i][j][g2];
//                    }
//                }
//        
//
//        // flux iteration (Gauss Seidel Method)
//        for (int counts = 1; counts <= inneriter; ++counts)
//            #pragma omp parallel for
//            for (int g = 0; g < group; ++g)
//                for (int i = 1; i <= Ni; ++i)
//                    for (int j = 1; j <= Nj; ++j)
//                        phi_new[i][j][g] = (source[i][j][g] - bidiagonali[i][j][g] * phi_new[i + 1][j][g] - bidiagonali[i - 1][j][g] * phi_new[i - 1][j][g] - bidiagonalj[i][j - 1][g] * phi_new[i][j - 1][g] - bidiagonalj[i][j][g] * phi_new[i][j + 1][g]) / diagonal[i][j][g];
//
//        // keff calculation
//        double numerator = 0.0f;
//        double denominator = 0.0f;
//        for (int i = 1; i <= Ni; ++i)
//            for (int j = 1; j <= Nj; ++j)
//                for (int g = 0; g < group; ++g) {
//                    numerator += Xf[i][j][g] * phi_new[i][j][g];
//                    denominator += Xf[i][j][g] * phi_old[i][j][g];
//                }
//        keff_new = keff_old * (numerator / denominator);
//
//        // keff error calculation
//        errk = abs((keff_new - keff_old) / keff_old);
//
//
//        // flux error calculation
//        double errp1, errp2;
//        for (int i = 1; i < Ni; ++i)
//            for (int j = 1; j <= Nj; ++j)
//                for (int g = 0; g < group; ++g) {
//                    errp1 = abs((phi_new[i][j][g] - phi_old[i][j][g]) / phi_old[i][j][g]);
//                    errp2 = abs((phi_new[i + 1][j][g] - phi_old[i + 1][j][g]) / phi_old[i + 1][j][g]);
//                    errp = max(errp1, errp2);
//                }
//
//        // reset flux
//        for (int i = 1; i <= Ni; ++i)
//            for (int j = 1; j <= Nj; ++j)
//                for (int g = 0; g < group; ++g)
//                    phi_old[i][j][g] = phi_new[i][j][g];
//
//        // reset keff
//        keff_old = keff_new;
//
//        // outer iteration
//        ++outiter;
//    }
//
//    // check end time
//    double end_time = clock();
//
//    // duration calculation
//    double duration = (end_time - start_time) * 1000 / CLOCKS_PER_SEC;
//
//    // show result
//    cout << "Mesh zone type [row X column] : " << Nj << " X " << Ni << "\n"
//        << "keff : " << keff_new << "\n"
//        << "keff error : " << errk << ", flux error : " << errp << "\n"
//        << "duration : " << duration << "ms" << ", iteration counts : " << outiter << "\n\n";
//
//    // write
//    ofstream fileoutput;
//    fileoutput.open("output.txt");
//    // write flux distribution for each group
//    for (int g = 0; g < group; ++g) {
//        fileoutput << "group " << g + 1 << "\n";
//        for (int i = 1; i <= Ni; ++i) {
//            for (int j = 1; j <= Nj; ++j) {
//                fileoutput << phi_new[i][j][g] << "\t";
//            }
//            fileoutput << "\n";
//        }
//        fileoutput << "\n\n";
//    }
//
//}