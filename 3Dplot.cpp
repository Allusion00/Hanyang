//#include <iostream>
//#include <fstream>
//#include <vector>
//
//using namespace std;
//
//// Define the dimensions of the neutron flux grid
//const int num_groups = 10;
//const int num_i_meshes = 20;
//const int num_j_meshes = 30;
//
//// Function to output neutron flux data to a file
//void outputFluxData(const vector<vector<vector<double>>>& phi) {
//    ofstream outputFile("neutron_flux_data.txt");
//    if (!outputFile.is_open()) {
//        cerr << "Error opening output file." << endl;
//        return;
//    }
//
//    // Output flux data in a format that can be easily plotted
//    for (int i = 0; i < num_i_meshes; ++i) {
//        for (int j = 0; j < num_j_meshes; ++j) {
//            // Output flux values for each group at mesh (i, j)
//            for (int g = 0; g < num_groups; ++g) {
//                outputFile << i << " " << j << " " << g << " " << phi[g][i][j] << endl;
//            }
//        }
//    }
//
//    outputFile.close();
//}
//
//int main() {
//    // Create a 3D vector to hold neutron flux values
//    vector<vector<vector<double>>> phi(num_groups, vector<vector<double>>(num_i_meshes, vector<double>(num_j_meshes)));
//
//    // Populate the neutron flux data (for demonstration purposes)
//    // You would typically read or calculate this data from elsewhere in your code
//    for (int g = 0; g < num_groups; ++g) {
//        for (int i = 0; i < num_i_meshes; ++i) {
//            for (int j = 0; j < num_j_meshes; ++j) {
//                // Example: Assign some arbitrary flux values
//                phi[g][i][j] = 0.1 * g + 0.01 * i + 0.001 * j;
//            }
//        }
//    }
//
//    // Output the neutron flux data to a file
//    outputFluxData(phi);
//
//    return 0;
//}
