//#include <iostream>		// I/O연산에 필요한 기본적인 기능
//#include <fstream>		// 파일 입출력 할 수 있게 만드는 라이브러리
//#include <string>
//#include <vector>
//using namespace std;	// 이름의 충돌을 막아주고, 프로그램의 가독성을 높여주는 코드
//// Visual 2008부터는 iostream.h를 기본으로 지원하지 않기 때문에 #include <iostream>을 하거나 cout를 쓰려면 이 코드를 써야한다.
//
//// 필요한 크기의 vector를 생성후 0.0으로 채운 뒤, input.file에서 정보 읽은 후 저장하는 방식
//// check start time
//double start_time = clock();
//
//int main() {
//    ifstream fileinput; // define ifstream as fileinput
//    fileinput.open("Assembly test input.txt");
//    string line; // define line as string character for getline
//
//    // Read CARD1 (# of Region, # of Energy group)
//    cout << "\n[CARD1]\n";
//    getline(fileinput, line);
//    int region, group, category; // # of Region, # of Energy group, # of Assembly category
//    fileinput >> region >> group >> category;
//    cout << "# of Region : " << region << "\n" << "# of group : " << group << "\n\n"; // output for checking CARD1
//
//    // Read CARD2 (Grid size [row column], Geometry size [row column], Assembly placement)
//    cout << "[CARD2]\n";
//    getline(fileinput, line);
//    getline(fileinput, line);
//    double gridsizei, gridsizej;
//    fileinput >> gridsizei >> gridsizej; // Grid size [cm]
//    int georowsize, geocolumnsize;
//    fileinput >> georowsize >> geocolumnsize;
//    vector<vector<int>> Geometry(georowsize, vector<int>(geocolumnsize, 0));
//    for (int row = 0; row < georowsize; ++row) {
//        getline(fileinput, line);
//        for (int column = 0; column < geocolumnsize; ++column) {
//            fileinput >> Geometry[row][column];
//            cout << Geometry[row][column] << " ";
//        }
//        cout << "\n";
//    }
//
//    // Read CARD3
//    cout << "\n[CARD3]\n";
//    getline(fileinput, line);
//    getline(fileinput, line);
//    int assemrowsize, assemcolumnsize;
//    fileinput >> assemrowsize >> assemcolumnsize;
//    vector<vector<vector<int>>> Assembly(category, vector<vector<int>>(assemrowsize, vector<int>(assemcolumnsize, 0)));
//    getline(fileinput, line);
//    getline(fileinput, line);
//    for (int c = 0; c < category; ++c) {
//        for (int row = 0; row < assemrowsize; ++row) {
//            for (int column = 0; column < assemcolumnsize; ++column) {
//                fileinput >> Assembly[c][row][column];
//                cout << Assembly[c][row][column] << " ";
//            }
//            cout << "\n";
//            getline(fileinput, line);
//        }
//        cout << "\n";
//        getline(fileinput, line);
//    }
//
//    int columnsize = assemcolumnsize * geocolumnsize;
//    int rowsize = assemrowsize * georowsize;
//    vector<vector<int>> grid(rowsize, vector<int>(columnsize, 0));
//    for (int geocolumn = 0; geocolumn < geocolumnsize; ++geocolumn) {
//        for (int column = geocolumn * assemcolumnsize, assemcolumn = 0; column < (geocolumn + 1) * assemcolumnsize, assemcolumn < assemcolumnsize; ++column, ++assemcolumn) {
//            for (int georow = 0; georow < georowsize; ++georow) {
//                for (int row = georow * assemrowsize, assemrow = 0; row < (georow + 1) * assemrowsize, assemrow < assemrowsize; ++row, ++assemrow) {
//                    grid[row][column] = Assembly[Geometry[georow][geocolumn]][assemrow][assemcolumn];
//                    cout << grid[row][column] << " ";
//                }
//            }
//            cout << "\n";
//        }
//    }
//
//    // Read CARD4 (# of fine meshes for each grid [horizontal vertical])
//    cout << "\n[CARD4]\n";
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
//
//    for (int i = 0; i <= 51; ++i)
//        cout << hi[i] << " " << hj[i] << "\n";
//    cout << "Mesh zone type [row X column] : " << Nj << " X " << Ni << "\n";
//
//    // Read CARD5 (Boundary Condition[Left Right] : 0 = reflective, 1 = flux zero, 2 = vacuum)
//    cout << "\n[CARD5]\n";
//    getline(fileinput, line);
//    getline(fileinput, line);
//    vector<int> BC(4, 0);
//    for (int i = 0; i < 4; ++i)
//        fileinput >> BC[i];
//    vector<double> b(4, 0.0f);
//    for (int i = 0; i < 4; ++i) {
//        if (BC[i] == 0)
//            b[i] = 10E-10f; // reflective
//        else if (BC[i] == 1)
//            b[i] = 10E+10f; // flux zero
//        else if (BC[i] == 2)
//            b[i] = 0.5; // vacuum
//        else
//            return 0;
//    }
//    vector<string> cardinal{ "North", "South", "West", "East" };
//    // output for checking CARD5
//    for (int i = 0; i < 4; ++i)
//        cout << cardinal[i] << " Boundary Condition : " << BC[i] << " (" << b[i] << ")\n\n";
//
//    // Read CARD6 (Convergence criteria : # of inner iterations, relative error of keff, max relative error of flux)
//    cout << "[CARD6]\n";
//    getline(fileinput, line);
//    getline(fileinput, line);
//    int inneriter;
//    double crik, crip;
//    fileinput >> inneriter >> crik >> crip;
//    // output for checking CARD6
//    cout << "# of inner iterations : " << inneriter << "\nkeff error : " << crik << "\nsource group error : " << crip << "\n\n";
//
//    // Read CARD7 (D, Xs(Fission), Xs(Absorption)) [region][group][data type]
//    cout << "[CARD7]\n";
//    int data = 3; // # of types of CARD7
//    vector<vector<vector<double>>> CARD7(region, vector<vector<double>>(group, vector<double>(data, 0.0f)));
//    getline(fileinput, line);
//    getline(fileinput, line);
//    for (int r = 0; r < region; ++r) {
//        for (int g = 0; g < group; ++g) {
//            getline(fileinput, line);
//            for (int i = 0; i < data; ++i) {
//                fileinput >> CARD7[r][g][i];
//                cout << CARD7[r][g][i] << " "; // output for checking CARD7
//            }
//            cout << "\n"; // output for checking CARD7
//        }
//        getline(fileinput, line);
//        cout << "\n"; // output for checking CARD7
//    }
//
//
//    // Read CARD8 (Xs(Only Down Scattering)) [region][High E group][Low E group]
//    cout << "[CARD8]\n";
//    vector<vector<vector<double>>> CARD8(region, vector<vector<double>>(group, vector<double>(group, 0.0f)));
//    getline(fileinput, line);
//    for (int r = 0; r < region; ++r) {
//        for (int lg = 0; lg < group; ++lg) {
//            getline(fileinput, line);
//            for (int hg = 0; hg < group; ++hg) {
//                fileinput >> CARD8[r][hg][lg];
//                cout << CARD8[r][lg][hg] << " "; // output for checking CARD8
//            }
//            cout << "\n"; // output for checking CARD8
//        }
//        getline(fileinput, line);
//        cout << "\n"; // output for checking CARD8
//    }
//
//    // Read CARD9 (xg : fraction of fission neutrons)
//    cout << "[CARD9]\n";
//    // xg [group]
//    vector<double> xg(group, 0.0f);
//    getline(fileinput, line);
//    for (int g = 0; g < group; ++g) {
//        fileinput >> xg[g];
//        cout << "group " << g + 1 << " 의 fission neutron fraction : " << xg[g] << "\n"; // output for checking CARD9
//    }
//    cout << "\n\n"; // output for checking CARD9
//
//    // Diffusion coefficient [cm] [group][i mesh][j mesh]
//    vector<vector<vector<double>>> D(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    // fission X section [cm-1] [group][i mesh][j mesh]
//    vector<vector<vector<double>>> Xf(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    // absorption X section [cm-1] [group][i mesh][j mesh]
//    vector<vector<vector<double>>> Xa(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    for (int g = 0; g < group; ++g)
//        for (int column = 0; column < columnsize; ++column)
//            for (int i = horizontalsum[column] + 1; i <= horizontalsum[column + 1]; ++i)
//                for (int row = 0; row < rowsize; ++row)
//                    for (int j = verticalsum[row] + 1; j <= verticalsum[row + 1]; ++j) {
//                        D[g][i][j] = CARD7[grid[row][column]][g][0];
//                        Xf[g][i][j] = CARD7[grid[row][column]][g][1];
//                        Xa[g][i][j] = CARD7[grid[row][column]][g][2];
//                    }
//
//    // scattering X section [cm-1] [high E group][low E group][i mesh][j mesh]
//    vector<vector<vector<vector<double>>>> Xs(group, vector<vector<vector<double>>>(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f))));
//    for (int hg = 0; hg < group; ++hg)
//        for (int lg = 0; lg < group; ++lg)
//            for (int column = 0; column < columnsize; ++column)
//                for (int i = horizontalsum[column] + 1; i <= horizontalsum[column + 1]; ++i)
//                    for (int row = 0; row < rowsize; ++row)
//                        for (int j = verticalsum[row] + 1; j <= verticalsum[row + 1]; ++j)
//                            Xs[hg][lg][i][j] = CARD8[grid[row][column]][hg][lg];
//
//    // removal X section [cm-1] [group][i mesh][j mesh]
//    vector<vector<vector<double>>> Xr(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    for (int g1 = 0; g1 < group; ++g1)
//        for (int i = 1; i <= Ni; ++i)
//            for (int j = 1; j <= Nj; ++j) {
//                for (int g2 = 0; g2 < group; ++g2) {
//                    if (g1 == g2)
//                        continue;
//                    Xr[g1][i][j] += Xs[g2][g1][i][j];
//                }
//                Xr[g1][i][j] += Xa[g1][i][j];
//            }
//
//
//    // initial effective k
//    double keff_old = 1.0f;
//    double keff_new = 0.0f;
//
//    // flux vector define [group][i mesh][j mesh]
//    vector<vector<vector<double>>> phi_old(group, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f)));
//    vector<vector<vector<double>>> phi_new(group, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f)));
//    for (int g = 0; g < group; ++g)
//        for (int i = 1; i <= Ni; ++i)
//            for (int j = 1; j <= Nj; ++j) {
//                phi_old[g][i][j] = 1.0f;
//                phi_new[g][i][j] = 1.0f;
//            }
//
//    // beta vector define [group][i mesh][j mesh]
//    vector<vector<vector<double>>> betai(group, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f))); // i direction beta vector
//    vector<vector<vector<double>>> betaj(group, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f))); // j direction beta vector
//    for (int g = 0; g < group; ++g) {
//        for (int j = 1; j <= Nj; ++j) {
//            betai[g][0][j] = b[2] / 2; // i direction West boundary beta
//            betai[g][Ni + 1][j] = b[3] / 2; // i direction East boundary beta
//        }
//        for (int i = 1; i <= Ni; ++i) {
//            betaj[g][i][0] = b[0] / 2; // j direction North boundary beta
//            betaj[g][i][Nj + 1] = b[1] / 2; // j direction South boundary beta
//            for (int j = 1; j <= Nj; ++j) {
//                betai[g][i][j] = D[g][i][j] / hi[i]; // i direction beta
//                betaj[g][i][j] = D[g][i][j] / hj[j]; // j direction beta
//            }
//        }
//    }
//
//    // D tilder vector define [group][i mesh][j mesh]
//    vector<vector<vector<double>>> D_tilderi(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    vector<vector<vector<double>>> D_tilderj(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    for (int g = 0; g < group; ++g)
//        for (int i = 0; i <= Ni; ++i)
//            for (int j = 0; j <= Nj; ++j) {
//                D_tilderi[g][i][j] = 2 * betai[g][i][j] * betai[g][i + 1][j] / (betai[g][i][j] + betai[g][i + 1][j]);
//                D_tilderj[g][i][j] = 2 * betaj[g][i][j] * betaj[g][i][j + 1] / (betaj[g][i][j] + betaj[g][i][j + 1]);
//            }
//
//    // Matrix M define diagonal, bidiagonal [group][i mesh][j mesh]
//    vector<vector<vector<double>>> diagonal(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    vector<vector<vector<double>>> bidiagonali(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    vector<vector<vector<double>>> bidiagonalj(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//    for (int g = 0; g < group; ++g)
//        for (int i = 1; i <= Ni; ++i)
//            for (int j = 1; j <= Nj; ++j) {
//                diagonal[g][i][j] = D_tilderi[g][i][j] / hi[i] + D_tilderi[g][i - 1][j] / hi[i] + D_tilderj[g][i][j - 1] / hj[j] + D_tilderj[g][i][j] / hj[j] + Xr[g][i][j];
//                bidiagonali[g][i][j] = -D_tilderi[g][i][j] / hi[i];
//                bidiagonalj[g][i][j] = -D_tilderj[g][i][j] / hj[j];
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
//        vector<vector<vector<double>>> source(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
//        for (int g1 = 0; g1 < group; ++g1) {
//            for (int i = 1; i <= Ni; ++i)
//                for (int j = 1; j <= Nj; ++j) {
//                    for (int g2 = 0; g2 < group; ++g2) {
//                        source[g1][i][j] += 1 / keff_old * Xf[g2][i][j] * phi_old[g2][i][j] * xg[g1];
//                        if (g1 == g2)
//                            continue;
//                        source[g1][i][j] += Xs[g1][g2][i][j] * phi_old[g2][i][j];
//                    }
//                }
//        }
//
//        // flux iteration (Gauss Seidel Method)
//        for (int counts = 1; counts <= inneriter; ++counts)
//            for (int g = 0; g < group; ++g)
//                for (int i = 1; i <= Ni; ++i)
//                    for (int j = 1; j <= Nj; ++j)
//                        phi_new[g][i][j] = (source[g][i][j] - bidiagonali[g][i][j] * phi_new[g][i + 1][j] - bidiagonali[g][i - 1][j] * phi_new[g][i - 1][j] - bidiagonalj[g][i][j - 1] * phi_new[g][i][j - 1] - bidiagonalj[g][i][j] * phi_new[g][i][j + 1]) / diagonal[g][i][j];
//
//        // keff calculation
//        double numerator = 0.0f;
//        double denominator = 0.0f;
//        for (int i = 1; i <= Ni; ++i)
//            for (int j = 1; j <= Nj; ++j)
//                for (int g = 0; g < group; ++g) {
//                    numerator += Xf[g][i][j] * phi_new[g][i][j];
//                    denominator += Xf[g][i][j] * phi_old[g][i][j];
//                }
//        keff_new = keff_old * (numerator / denominator);
//
//        // keff error calculation
//        errk = abs((keff_new - keff_old) / keff_old);
//
//
//        // flux error calculation
//        double errp1, errp2;
//        for (int g = 0; g < group; ++g)
//            for (int i = 1; i < Ni; ++i)
//                for (int j = 1; j <= Nj; ++j) {
//                    errp1 = abs((phi_new[g][i][j] - phi_old[g][i][j]) / phi_old[g][i][j]);
//                    errp2 = abs((phi_new[g][i + 1][j] - phi_old[g][i + 1][j]) / phi_old[g][i + 1][j]);
//                    errp = max(errp1, errp2);
//                }
//
//        // reset flux
//        for (int g = 0; g < group; ++g)
//            for (int i = 1; i <= Ni; ++i)
//                for (int j = 1; j <= Nj; ++j)
//                    phi_old[g][i][j] = phi_new[g][i][j];
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
//                fileoutput << phi_new[g][i][j] << "\t";
//            }
//            fileoutput << "\n";
//        }
//        fileoutput << "\n\n";
//    }
//
//}