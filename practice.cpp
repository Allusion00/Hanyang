//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <string>
//
//using namespace std;
//
//int main() {
//    // Open the file
//    ifstream inputFile("input.txt");
//
//    // Read the file line by line
//    string line;
//
//    // read row, column data for region vector
//    getline(inputFile, line);
//    int row, column;
//    inputFile >> row >> column;
//
//    // define region vector [row][column] = [material]
//    vector<vector<double>> region(row, vector<double>(column, 0));
//    for (int i = 0; i < row; ++i) {
//        getline(inputFile, line);
//        for (int j = 0; j < column; ++j) {
//            inputFile >> region[i][j];
//            cout << region[i][j] << " ";
//        }
//        cout << "\n";
//    }
//
//    // define mesh vector
//    double meshsize = 0.5;
//    int meshnum = 10 / meshsize;
//    int imesh = meshnum * column + 1;
//    int jmesh = meshnum * row + 1;
//    vector<vector<double>> mesh(imesh, vector<double>(jmesh, 0.0f));
//    
//}
