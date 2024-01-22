#include "define.h"
#include "Data.h"
#include "Masterpiece.h"

int main() {
	double Size = 0.0;
	vector<vector<int>> Geometry_Data;
	vector<vector<int>> Mesh_Data;
	vector<int> Boundary_Data;
	vector<vector<vector<double>>> DFA_Data;
	vector<vector<vector<double>>> Xs_Data;
	vector<double> Xi_Data;
	vector<double> Convcrit_Data;
	Data A;



	A.read(Size, Geometry_Data, Mesh_Data, Boundary_Data, DFA_Data, Xs_Data, Xi_Data, Convcrit_Data); // Input file 

	For2D B(Size, Geometry_Data, Mesh_Data, Boundary_Data, DFA_Data, Xs_Data, Xi_Data, Convcrit_Data);

	B.Gauss_Seidal();


	return 0;
}
