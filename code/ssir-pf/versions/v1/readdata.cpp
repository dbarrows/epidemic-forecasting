#include <string.h>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

int getNumLines(string filename) {

	int numLines = 0;
	int ch;

    FILE *infile = fopen(filename.c_str(), "r");

    while ( EOF != (ch = getc(infile)) )
        if (ch == '\n')
            ++numLines;
    
    return numLines;

}

int getDim(string line) {

	int dim = 0;
	int dummyVal;

	stringstream lineStream(line);

	while (lineStream >> dummyVal)
		++dim;

	return dim;

}

int * getData(string filename, int * xdim, int * ydim) {

	int xlen, ylen;
	int dctr = 0;

	string line;
	ifstream fileStream(filename.c_str());

	getline(fileStream, line);

	xlen = getDim(line);
	ylen = getNumLines(filename);

	int memsize = xlen*ylen*sizeof(int);

	int * data = (int *) malloc (memsize);

	int val;

	// process first line
	stringstream lineStream(line);
    while(lineStream >> val) {
        data[dctr] = val;
        dctr++;
    }

    // process remaining lines
    while( getline(fileStream, line) ) {

        stringstream lineStream(line);
        while(lineStream >> val) {
            data[dctr] = val;
            dctr++;
        }

    }

    // pass back dimensions
    if (xdim != NULL && ydim != NULL) {
	    *xdim = xlen;
		*ydim = ylen;
	}

    return data;

}