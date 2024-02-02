//#include <iostream>
//#include <vector>
//
//using namespace std;
//
//vector<vector<vector<double>>> geo(const vector<vector<int>>& mesh, int x, int y) {
//    // Initialize geo with specific size (region, x, y) and initial values of 0.0
//    vector<vector<vector<double>>> result(mesh.size(), vector<vector<double>>(x, vector<double>(y, 0.0)));
//
//    // Convert each 2D vector from mesh to 3D and store in result
//    for (size_t i = 0; i < mesh.size(); ++i) {
//        for (size_t j = 0; j < mesh[i].size(); ++j) {
//            result[i][j][2] = mesh[i][j]; // Assuming z is the third dimension
//        }
//    }
//
//    return result;
//}
//
//int main() {
//    int region = 7; // Define the number of 2D vectors (depth)
//    int x = 5;      // Define the number of rows
//    int y = 3;      // Define the number of columns
//
//    // Sample mesh
//    vector<vector<int>> mesh = {
//        {1, 2, 3, 4, 5},
//        {6, 7, 8, 9, 10},
//        {11, 12, 13, 14, 15},
//        {16, 17, 18, 19, 20},
//        {21, 22, 23, 24, 25},
//        {26, 27, 28, 29, 30},
//        {31, 32, 33, 34, 35}
//    };
//
//    // Call the geo function with mesh and desired dimensions
//    vector<vector<vector<double>>> result = geo(mesh, x, y);
//
//    // Display the contents of the result vector
//    for (const auto& v3d : result) {
//        for (const auto& v2d : v3d) {
//            for (const auto& v : v2d) {
//                cout << v << " ";
//            }
//            cout << endl;
//        }
//        cout << endl;
//    }
//
//    return 0;
//}
