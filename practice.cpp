#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main() {
    // Open the file
    ifstream inputFile("input.txt");

    // Read the file line by line
    string line;

    // read row, column data for region vector
    getline(inputFile, line);
    int row, column;
    inputFile >> row >> column;

    // define grid vector [row][column] = [material]
    vector<vector<double>> grid(row, vector<double>(column, 0));
    for (int i = 0; i < row; ++i) {
        getline(inputFile, line);
        for (int j = 0; j < column; ++j) {
            inputFile >> grid[i][j];
            cout << grid[i][j] << " ";
        }
        cout << "\n";
    }


    // define mesh vector
    double meshsize = 10;
    int meshnum = 10 / meshsize;
    int Ni = meshnum * column + 1;
    int Nj = meshnum * row + 1;
    

    // Read CARD6 (D, Xs(Fission), Xs(Absorption))
    cout << "[CARD6]\n";
    int material = 2; // UO2, MOX4.3
    int n = 3; // # of types of Data1
    int g = 7; // 7 group input
    vector<vector<vector<double>>> Data1(material, vector<vector<double>>(g, vector<double>(n, 0.0f)));
    getline(inputFile, line);
    getline(inputFile, line);
    for (int i = 0; i < material; ++i) {
        for (int j = 0; j < g; ++j) {
            getline(inputFile, line);
            for (int k = 0; k < n; ++k) {
                inputFile >> Data1[i][j][k];
                cout << Data1[i][j][k] << " "; // output for checking CARD6
            }
            cout << "\n"; // output for checking CARD6
        }
        getline(inputFile, line);
        cout << "\n"; // output for checking CARD6
    }


    
    // Diffusion coefficient [cm] [group][i mesh][j mesh]
    vector<vector<vector<double>>> D(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    vector<vector<vector<double>>> Xf(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    vector<vector<vector<double>>> Xa(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    for (int group = 0; group < g; ++group)
        for (int i = 0; i < column; ++i)
            for (int imesh = i * meshnum + 1; imesh <= (i + 1) * meshnum; ++imesh)
                for (int j = 0; j < row; ++j)
                    for (int jmesh = j * meshnum + 1; jmesh <= (j + 1) * meshnum; ++jmesh) {
                        D[group][imesh][jmesh] = Data1[grid[i][j]][group][0];
                        Xf[group][imesh][jmesh] = Data1[grid[i][j]][group][1];
                        Xa[group][imesh][jmesh] = Data1[grid[i][j]][group][2];
                    }

    // Read CARD7 (Xs(Only Down Scattering))
    cout << "[CARD7]\n";
    vector<vector<vector<double>>> Data2(material, vector<vector<double>>(g, vector<double>(g, 0.0f)));
    getline(inputFile, line);
    for (int i = 0; i < material; ++i) {
        for (int j = 0; j < g; ++j) {
            getline(inputFile, line);
            for (int k = 0; k < g; ++k) {
                inputFile >> Data2[i][k][j];
                cout << Data2[i][j][k] << " "; // output for checking CARD7
            }
            cout << "\n"; // output for checking CARD7
        }
        getline(inputFile, line);
        cout << "\n"; // output for checking CARD7
    }

    // Scattering X section [cm-1] [High E group][Low E group][i mesh][j mesh]
    vector<vector<vector<vector<double>>>> Xs(g, vector<vector<vector<double>>>(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f))));
    for (int hgroup = 0; hgroup < g; ++hgroup)
        for (int lgroup = 0; lgroup < g; ++lgroup)
            for (int i = 0; i < column; ++i)
                for (int imesh = i * meshnum + 1; imesh <= (i + 1) * meshnum; ++imesh)
                    for (int j = 0; j < row; ++j)
                        for (int jmesh = j * meshnum + 1; jmesh <= (j + 1) * meshnum; ++jmesh)
                            Xs[hgroup][lgroup][imesh][jmesh] = Data2[grid[i][j]][hgroup][lgroup];


    // Removal X section [cm-1] [group][i mesh][j mesh]
    vector<vector<vector<double>>> Xr(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    vector<vector<vector<double>>> Xs_sum(g, vector<vector<double>>(Ni + 1, vector<double>(Nj + 1, 0.0f)));
    for (int i = 0; i < g; ++i)
        for (int j = 1; j <= i; ++j)
            for (int k = 1; k <= Nj; ++k) {
                for (int l = 0; l < g; ++l)
                    Xs_sum[i][j][k] += Xs[l][i][j][k];
                Xr[i][j][k] = Xa[i][j][k] + Xs_sum[i][j][k];
            }

}
