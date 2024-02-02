#include <iostream>
#include <vector>

using namespace std;

int main() {
	int row = 4;
	int column = 5;
	vector<vector<int>> region(row, vector<int>(column, 0));

	int meshnum = 20;
	int imesh = meshnum * column + 1;
	int jmesh = meshnum * row + 1;
	vector<vector<double>> mesh(imesh, vector<double>(jmesh, 0.0f));

	
}