#include <iostream>

#include "readdata.h"

using namespace std;

int main () {

	string filename = "data/ssir_S10_T100_R3_B0.0012_r0.2_TRUE.dat";

	int xdim, ydim;

	int * data = getData(filename, &xdim, &ydim);

	cout << "xdim:	" << xdim << endl;
	cout << "ydim:	" << ydim << endl;

	int nLines = 10;

	for (int y = 0; y < nLines; y++) {
		for (int x = 0; x < nLines; x++) {
			cout << data[y*xdim + x] << "\t";
		}
		cout << endl;
	}

}