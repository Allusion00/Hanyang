#include <iostream>		// I/O연산에 필요한 기본적인 기능
#include <fstream>		// 파일 입출력 할 수 있게 만드는 라이브러리
#include <string>
#include <vector>
#include <algorithm>
using namespace std;	// 이름의 충돌을 막아주고, 프로그램의 가독성을 높여주는 코드
// Visual 2008부터는 iostream.h를 기본으로 지원하지 않기 때문에 #include <iostream>을 하거나 cout를 쓰려면 이 코드를 써야한다.

// 필요한 크기의 vector를 생성후 0.0으로 채운 뒤, input.file에서 정보 읽은 후 저장하는 방식
// check start time
double start_time = clock();

int main() {
    ifstream fileinput; // define ifstream as fileinput

    // Input file select option
    int groupinput;
    cout << "[Select input file]\n1. 2D1G\n2. 2D2G\n3. 2D7G\n4. C5G7\n5. Kang\n6. BIBLIS\n7. IAEA PWR\n";
    cin >> groupinput;
    if (groupinput == 1)
        fileinput.open("2D1G input.txt");
    else if (groupinput == 2)
        fileinput.open("2D2G input.txt");
    else if (groupinput == 3)
        fileinput.open("2D7G input.txt");
    else if (groupinput == 4)
        fileinput.open("C5G7 input.txt");
    else if (groupinput == 5)
        fileinput.open("Kang input.txt");
    else if (groupinput == 6)
        fileinput.open("BIBLIS input.txt");
    else if (groupinput == 7)
        fileinput.open("IAEA PWR input.txt");
    else {
        cout << "Please select between examples";
        return 0;
    }

    string line; // define line as string character for getline

    // Read CARD1 (# of Region, # of Energy group)
    cout << "\n[CARD1]\n";
    getline(fileinput, line);
    int region, group; // # of Region, # of Energy group
    fileinput >> region >> group;
    cout << "# of Region : " << region << "\n" << "# of group : " << group << "\n\n"; // output for checking CARD1

    // Read CARD2 (Grid size [row column], Geometry size [row column], Geometry)
    cout << "[CARD2]\n";
    getline(fileinput, line);
    getline(fileinput, line);
    double gridsizei, gridsizej;
    fileinput >> gridsizei >> gridsizej; // Grid size [cm]
    int rowsize, columnsize;
    fileinput >> rowsize >> columnsize;
    getline(fileinput, line);

    vector<string> grid_input(rowsize);
    for (int row = 0; row < rowsize; ++row) {
        getline(fileinput, grid_input[row]);
        grid_input[row].erase(remove(grid_input[row].begin(), grid_input[row].end(), ' '), grid_input[row].end());
        cout << grid_input[row] << "\n";
    }

    vector<vector<int>> grid(rowsize, vector<int>(columnsize));
    for (int row = 0; row < rowsize; ++row) {
        for (int column = 0; column < columnsize; ++column) {
            if (grid_input[row][column] == 'f')
                grid[row][column] = 11;
            else if (grid_input[row][column] == 'v')
                grid[row][column] = 12;
            else
                grid[row][column] = int(grid_input[row][column] - '0');
            cout << grid[row][column] << " ";
        }
        cout << "\n";
    }


    // Read CARD3 (# of fine meshes for each grid [horizontal vertical])
    cout << "\n[CARD3]\n";
    getline(fileinput, line);
    vector<int> horizontal(columnsize, 0);
    vector<int> vertical(rowsize, 0);
    int Ni = 0;
    int Nj = 0;
    for (int column = 0; column < columnsize; ++column) {
        fileinput >> horizontal[column];
    }
    getline(fileinput, line);
    for (int row = 0; row < rowsize; ++row) {
        fileinput >> vertical[row];
    }
    vector<int> horizontalsum(columnsize + 1, 0);
    vector<int> verticalsum(rowsize + 1, 0);
    for (int column = 1; column <= columnsize; ++column)
        for (int i = 0; i <= column - 1; ++i)
            horizontalsum[column] += horizontal[i];
    for (int row = 1; row <= rowsize; ++row)
        for (int i = 0; i <= row - 1; ++i)
            verticalsum[row] += vertical[i];
    Ni = horizontalsum[columnsize];
    Nj = verticalsum[rowsize];
    vector<double> hi(Ni + 1, 0.0f);
    vector<double> hj(Nj + 1, 0.0f);
    for (int column = 0; column < columnsize; ++column)
        for (int i = horizontalsum[column] + 1; i <= horizontalsum[column + 1]; ++i)
            for (int row = 0; row < rowsize; ++row)
                for (int j = verticalsum[row] + 1; j <= verticalsum[row + 1]; ++j) {
                    hi[i] = gridsizei / horizontal[column];
                    hj[j] = gridsizej / vertical[row];
                }
    cout << "Mesh zone type [row X column] : " << Nj << " X " << Ni << "\n";

    // Read CARD4 (Boundary Condition[Left Right] : 10 = reflective, 11 = flux zero, 12 = vacuum)
    cout << "\n[CARD4]\n";
    getline(fileinput, line);
    getline(fileinput, line);
    vector<int> BC(4, 0);
    for (int i = 0; i < 4; ++i)
        fileinput >> BC[i];
    vector<double> b(4, 0.0f);
    for (int i = 0; i < 4; ++i) {
        if (BC[i] == 10)
            b[i] = 10E-10f; // reflective
        else if (BC[i] == 11)
            b[i] = 10E+10f; // flux zero
        else if (BC[i] == 12)
            b[i] = 0.5; // vacuum
        else
            return 0;
    }
    vector<string> cardinal{ "North", "South", "West", "East" };
    // output for checking CARD2
    for (int i = 0; i < 4; ++i)
        cout << cardinal[i] << " Boundary Condition : " << BC[i] << " (" << b[i] << ")\n\n";

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
    int data = 3; // # of types of CARD6
    vector<vector<vector<double>>> CARD6(region, vector<vector<double>>(group, vector<double>(data, 0.0f)));
    getline(fileinput, line);
    getline(fileinput, line);
    for (int r = 0; r < region; ++r) {
        for (int g = 0; g < group; ++g) {
            getline(fileinput, line);
            for (int i = 0; i < data; ++i) {
                fileinput >> CARD6[r][g][i];
                cout << CARD6[r][g][i] << " "; // output for checking CARD6
            }
            cout << "\n"; // output for checking CARD6
        }
        getline(fileinput, line);
        cout << "\n"; // output for checking CARD6
    }


    // Read CARD7 (Xs(Only Down Scattering)) [region][High E group][Low E group]
    cout << "[CARD7]\n";
    vector<vector<vector<double>>> CARD7(region, vector<vector<double>>(group, vector<double>(group, 0.0f)));
    getline(fileinput, line);
    for (int r = 0; r < region; ++r) {
        for (int lg = 0; lg < group; ++lg) {
            getline(fileinput, line);
            for (int hg = 0; hg < group; ++hg) {
                fileinput >> CARD7[r][hg][lg];
                cout << CARD7[r][lg][hg] << " "; // output for checking CARD7
            }
            cout << "\n"; // output for checking CARD7
        }
        getline(fileinput, line);
        cout << "\n"; // output for checking CARD7
    }

    // Read CARD8 (xg : fraction of fission neutrons)
    cout << "[CARD8]\n";
    // xg [group]
    vector<double> xg(group, 0.0f);
    getline(fileinput, line);
    for (int g = 0; g < group; ++g) {
        fileinput >> xg[g];
        cout << "group " << g + 1 << " 의 fission neutron fraction : " << xg[g] << "\n"; // output for checking CARD8
    }
    cout << "\n\n"; // output for checking CARD8


    // Diffusion coefficient [cm] [group][i mesh][j mesh]
    vector<vector<vector<double>>> D(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    // fission X section [cm-1] [group][i mesh][j mesh]
    vector<vector<vector<double>>> Xf(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    // absorption X section [cm-1] [group][i mesh][j mesh]
    vector<vector<vector<double>>> Xa(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    // Mesh data [group][imesh][jmesh]
    vector<vector<vector<int>>> Data(group, vector<vector<int>>(Ni + 1, vector<int>(Nj + 1, 1)));
    for (int g = 0; g < group; ++g)
        for (int column = 0; column < columnsize; ++column)
            for (int i = horizontalsum[column] + 1; i <= horizontalsum[column + 1]; ++i)
                for (int row = 0; row < rowsize; ++row)
                    for (int j = verticalsum[row] + 1; j <= verticalsum[row + 1]; ++j) {
                        if (grid[row][column] == 11 || grid[row][column] == 12) {
                            Data[g][i][j] = 0;
                            continue;
                        }
                        D[g][i][j] = CARD6[grid[row][column]][g][0];
                        Xf[g][i][j] = CARD6[grid[row][column]][g][1];
                        Xa[g][i][j] = CARD6[grid[row][column]][g][2];
                    }

    // scattering X section [cm-1] [high E group][low E group][i mesh][j mesh]
    vector<vector<vector<vector<double>>>> Xs(group, vector<vector<vector<double>>>(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f))));
    for (int hg = 0; hg < group; ++hg)
        for (int lg = 0; lg < group; ++lg)
            for (int column = 0; column < columnsize; ++column)
                for (int i = horizontalsum[column] + 1; i <= horizontalsum[column + 1]; ++i)
                    for (int row = 0; row < rowsize; ++row)
                        for (int j = verticalsum[row] + 1; j <= verticalsum[row + 1]; ++j) {
                            if (grid[row][column] == 11 || grid[row][column] == 12)
                                continue;
                            Xs[hg][lg][i][j] = CARD7[grid[row][column]][hg][lg];
                        }

    // removal X section [cm-1] [group][i mesh][j mesh]
    vector<vector<vector<double>>> Xr(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    for (int g1 = 0; g1 < group; ++g1)
        for (int i = 1; i <= Ni; ++i)
            for (int j = 1; j <= Nj; ++j) {
                for (int g2 = 0; g2 < group; ++g2) {
                    if (g1 == g2)
                        continue;
                    Xr[g1][i][j] += Xs[g2][g1][i][j];
                }
                Xr[g1][i][j] += Xa[g1][i][j];
            }


    // initial effective k
    double keff_old = 1.0f;
    double keff_new = 0.0f;

    // flux vector define [group][i mesh][j mesh]
    vector<vector<vector<double>>> phi_old(group, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f)));
    vector<vector<vector<double>>> phi_new(group, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f)));
    for (int g = 0; g < group; ++g)
        for (int i = 1; i <= Ni; ++i)
            for (int j = 1; j <= Nj; ++j) {
                phi_old[g][i][j] = 1.0f;
                phi_new[g][i][j] = 1.0f;
            }

    // beta vector define [group][i mesh][j mesh]
    vector<vector<vector<double>>> betai(group, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f))); // i direction beta vector
    vector<vector<vector<double>>> betaj(group, vector<vector<double>>(Ni + 2, vector<double>(Nj + 2, 0.0f))); // j direction beta vector
    for (int g = 0; g < group; ++g) {
        for (int j = 1; j <= Nj; ++j) {
            betai[g][0][j] = b[2] / 2; // i direction West boundary beta
            betai[g][Ni + 1][j] = b[3] / 2; // i direction East boundary beta
        }
        for (int i = 1; i <= Ni; ++i) {
            betaj[g][i][0] = b[0] / 2; // j direction North boundary beta
            betaj[g][i][Nj + 1] = b[1] / 2; // j direction South boundary beta
            for (int j = 1; j <= Nj; ++j) {
                if (Data[g][i][j] == 0)
                    continue;
                betai[g][i][j] = D[g][i][j] / hi[i]; // i direction beta
                betaj[g][i][j] = D[g][i][j] / hj[j]; // j direction beta
            }
        }
    }

    // D tilder vector define [group][i mesh][j mesh]
    vector<vector<vector<double>>> D_tilderi(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    vector<vector<vector<double>>> D_tilderj(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    for (int g = 0; g < group; ++g)
        for (int i = 0; i <= Ni; ++i)
            for (int j = 0; j <= Nj; ++j) {
                if (Data[g][i][j] == 0)
                    continue;
                D_tilderi[g][i][j] = 2 * D[g][i][j] * betai[g][i + 1][j] / (betai[g][i][j] + betai[g][i + 1][j]);
                D_tilderj[g][i][j] = 2 * betaj[g][i][j] * betaj[g][i][j + 1] / (betaj[g][i][j] + betaj[g][i][j + 1]);
            }

    // Matrix M define diagonal, bidiagonal [group][i mesh][j mesh]
    vector<vector<vector<double>>> diagonal(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    vector<vector<vector<double>>> bidiagonali(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    vector<vector<vector<double>>> bidiagonalj(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    for (int g = 0; g < group; ++g)
        for (int i = 1; i <= Ni; ++i)
            for (int j = 1; j <= Nj; ++j) {
                if (Data[g][i][j] == 0)
                    continue;
                diagonal[g][i][j] = D_tilderi[g][i][j] / hi[i] + D_tilderi[g][i - 1][j] / hi[i] + D_tilderj[g][i][j - 1] / hj[j] + D_tilderj[g][i][j] / hj[j] + Xr[g][i][j];
                bidiagonali[g][i][j] = -D_tilderi[g][i][j] / hi[i];
                bidiagonalj[g][i][j] = -D_tilderj[g][i][j] / hj[j];
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
        vector<vector<vector<double>>> source(group, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
        for (int g1 = 0; g1 < group; ++g1) {
            for (int i = 1; i <= Ni; ++i)
                for (int j = 1; j <= Nj; ++j) {
                    for (int g2 = 0; g2 < group; ++g2) {
                        if (Data[g2][i][j] == 0)
                            continue;
                        source[g1][i][j] += 1 / keff_old * Xf[g2][i][j] * phi_old[g2][i][j] * xg[g1];
                        if (g1 == g2)
                            continue;
                        source[g1][i][j] += Xs[g1][g2][i][j] * phi_old[g2][i][j];
                    }
                }
        }

        // flux iteration (Gauss Seidel Method)
        for (int counts = 1; counts <= inneriter; ++counts)
            for (int g = 0; g < group; ++g)
                for (int i = 1; i <= Ni; ++i)
                    for (int j = 1; j <= Nj; ++j) {
                        if (Data[g][i][j] == 0)
                            continue;
                        phi_new[g][i][j] = (source[g][i][j] - bidiagonali[g][i][j] * phi_new[g][i + 1][j] - bidiagonali[g][i - 1][j] * phi_new[g][i - 1][j] - bidiagonalj[g][i][j - 1] * phi_new[g][i][j - 1] - bidiagonalj[g][i][j] * phi_new[g][i][j + 1]) / diagonal[g][i][j];
                    }


        // keff calculation
        double numerator = 0.0f;
        double denominator = 0.0f;
        for (int i = 1; i <= Ni; ++i)
            for (int j = 1; j <= Nj; ++j)
                for (int g = 0; g < group; ++g) {
                    if (Data[g][i][j] == 0)
                        continue;
                    numerator += Xf[g][i][j] * phi_new[g][i][j];
                    denominator += Xf[g][i][j] * phi_old[g][i][j];
                }
        keff_new = keff_old * (numerator / denominator);

        // keff error calculation
        errk = abs((keff_new - keff_old) / keff_old);


        // flux error calculation
        double errp1, errp2;
        for (int g = 0; g < group; ++g)
            for (int i = 1; i < Ni; ++i)
                for (int j = 1; j <= Nj; ++j) {
                    if (Data[g][i][j] == 0)
                        continue;
                    errp1 = abs((phi_new[g][i][j] - phi_old[g][i][j]) / phi_old[g][i][j]);
                    errp2 = abs((phi_new[g][i + 1][j] - phi_old[g][i + 1][j]) / phi_old[g][i + 1][j]);
                    errp = max(errp1, errp2);
                }

        // reset flux
        for (int g = 0; g < group; ++g)
            for (int i = 1; i <= Ni; ++i)
                for (int j = 1; j <= Nj; ++j)
                    phi_old[g][i][j] = phi_new[g][i][j];

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
        << "keff : " << keff_new << "\n"
        << "keff error : " << errk << ", flux error : " << errp << "\n"
        << "duration : " << duration << "ms" << ", iteration counts : " << outiter << "\n\n";

    // write
    ofstream fileoutput;
    fileoutput.open("output.txt");
    // write flux distribution for each group
    for (int g = 0; g < group; ++g) {
        fileoutput << "group " << g + 1 << "\n";
        for (int i = 1; i <= Ni; ++i) {
            for (int j = 1; j <= Nj; ++j) {
                fileoutput << phi_new[g][i][j] << "\t";
            }
            fileoutput << "\n";
        }
        fileoutput << "\n\n";
    }

}