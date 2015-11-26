#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include "rand.h"
#include "timer.h"

#define TS 		100			// time to simulate over
#define R0 		3.0			// infectiousness
#define r 		2e-1		// recovery rate
#define N 		500 		// population size
#define B 		R0*r/N		// transmission factor
#define pnoise 	1.0			// process noise case multiplier
#define merr 	10.0  		// measurement error, 1% of population
#define ALP 	1e-2		// cooling parameter for neighbour interactions
#define SpW		10			// number of Euler steps to take per week
#define eta		0.5			// Beta drift cooling
#define berr 	1e-3 		// Beta shocks
#define Phi 	0.5			// Difficulty of transmission between cells, 0 = Impossible, 1 = no barrier

struct SIR {
	float S;
	float I;
	float R;
};

struct Param {
	float Beta;
	float alpha;
	float Beta_bar;
};

void exp_euler_SSIR(float h, SIR * population, Param * params, float * nInfec_last, int cID, int neibIDs[][8], int * nNeibs);
double dist (int x1, int y1, int x2, int y2);

using namespace std;

int main (int argc, char ** argv) {

	// Parse input arg(s) ***************

	int dim;
	if (argc == 2)
		dim = atoi(argv[1]);
	else {
		cout << "Missing dimension argument" << endl;
		return 1;
	}

	// **********************************


	// Parameters and such **********************************************************************************

	int init_infected = 1;							// Number of initial infections in seed cells

	int T = max( (int) round(10*dim/R0), TS );		// Time to simulate over
	int nCells = dim*dim;							// Number of cells (a DIMxDIM square)

	SIR population[nCells];							// SIR compartment counts for each cell
	float nInfec_last[nCells];						// stores last infection counts between propagation steps
	int adjMat[nCells][8];							// Each row holds the indices of that cell's neighbours

	srand (time(NULL));								// Seed PRNG with system time

	int nInfecCells = max(dim/20, 1);				// Initial number of cells to infect
	int out_start_loc[nInfecCells];					// Locations to start infections

	float neiInf[nCells][8];						// used to store current neighbour infected counts

	// ******************************************************************************************************


	// Print simulation stats*************************

	cout << "Simulation stats" << endl;
	cout << "----------------" << endl;
	cout << "Dims      " << dim << "x" << dim << endl;
	cout << "Time      " << T << endl;
	cout << "R0        " << R0 << endl;
	cout << "r (rec)   " << r << endl;

	// ***********************************************


	// Generate Cell parameters ***************

	Param params[nCells];
	for (int cell = 0; cell < nCells; cell++) {
		params[cell].Beta 	= B + (float) B/10.0*randn();
		params[cell].alpha 	= ALP+ (float) ALP/10.0*randn();
		params[cell].Beta_bar = params[cell].Beta;
	}

	// ****************************************


	// populate adjacentcy matrix **************************************************************************

	int inds_check[8] 	= {-(dim+1), -(dim), -(dim-1), -1, 1, (dim-1), (dim), (dim+1)};
	int nNeibVec[nCells];

	for (int cell = 0; cell < nCells; cell++) {

		int neibCtr = 0;

		int i_cell = cell / dim;
		int j_cell = cell % dim;

		for (int j = 0; j < 8; j++) {

			int nei_ind = cell + inds_check[j];
			int i_ref = nei_ind / dim;
			int j_ref = nei_ind % dim;

			if (0 <= nei_ind && nei_ind < nCells && abs(i_cell - i_ref) <= 1 && abs(j_cell - j_ref) <= 1 ) {
				adjMat[cell][neibCtr] = nei_ind;
				neibCtr++;
			}

			if (neibCtr < 8) {
				for (int i = neibCtr; i < 8; i++)
					adjMat[cell][i] = -1;
			}

		}

		nNeibVec[cell] = neibCtr;

	}

	// *****************************************************************************************************


	// Draw [nInfecCells] random numbers in the range [0,nCells],
	// and check for duplicates
	int f_ctr = 0;
	while (f_ctr < nInfecCells) {

		int canCell = rand() % nCells;
		int used_flag = 0;
		for (int j = 0; j < f_ctr; j++) {
			if (canCell == out_start_loc[j])
				used_flag = 1;
		}
		if (!used_flag){
			out_start_loc[f_ctr] = canCell;
			f_ctr++;
		}

	}

	// write each starting state as no infection
	for (int i = 0; i < nCells; i++) {

		population[i].S = N;
		population[i].I = 0;
		population[i].R = 0;

	}

	// overwrite no infection state with infection state
	for (int i = 0; i < nInfecCells; i++) {

		population[out_start_loc[i]].S = N - init_infected;
		population[out_start_loc[i]].I = init_infected;
		population[out_start_loc[i]].R = 0;

	}


	// Data generation and output ***************************************************************

	// Output filename strings for data, true model state, parameter values 
	stringstream outfile_d, outfile_t, out_p;

	outfile_d << "ssir";
	outfile_t << "ssir";

	out_p << "_S" << dim;
	out_p << "_T" << T;
	out_p << "_R" << R0;
	out_p << "_B" << B;
	out_p << "_r" << r;

	outfile_d << out_p.str() << "_DATA.dat";
	outfile_t << out_p.str() << "_TRUE.dat";

	ofstream file_d, file_t;
	file_d.open( outfile_d.str().c_str() );
	file_t.open( outfile_t.str().c_str() );

	for (int t = 0; t < T; t++) {

		for (int st = 0; st < SpW; st++) {

			for (int cell = 0; cell < nCells; cell++) {
				nInfec_last[cell] = population[cell].I;
			}

			// propagate cells forward
			for (int cell = 0; cell < nCells; cell++)
				exp_euler_SSIR(1.0/SpW, population, params, nInfec_last, cell, adjMat, nNeibVec);

		}


		// perturb states (observation noise) and save as data for the week
		for (int cell = 0; cell < nCells; cell++) {

			int truepoint = round( population[cell].I );
			int datapoint = round( truepoint + merr*randn() );
			if(datapoint < 0)
				datapoint = 0;

			file_t << truepoint;
			file_d << datapoint;

			if ( (cell-(dim-1)) % dim == 0) {
				file_t << endl;
				file_d << endl;
			} else {
				file_t << " ";
				file_d << " ";
			}

		}

	}

	file_t.close();
	file_d.close();

	// ******************************************************************************************


	return 0;

}


void exp_euler_SSIR(float h, SIR * population, Param * params, float * nInfec_last, int cID, int neibIDs[][8], int * nNeibs) {

	float Beta 		= params[cID].Beta;
	float alpha 	= params[cID].alpha;
	float Beta_bar 	= params[cID].Beta_bar;

	float S = population[cID].S;
	float I = population[cID].I;
	float R = population[cID].R;

	// save initial states
	float S0 = S;
	float I0 = I;

	// number of neighbours
	int n = nNeibs[cID];

	// calculate the sum of Beta*I products from each neighbour 
	float neiInf_sum = 0;
	for (int i = 0; i < n; i++) {
		int nID = neibIDs[cID][i];
		neiInf_sum += params[nID].Beta * nInfec_last[nID];
	}

	// update beta, make sure it doesn't go negative
	Beta = Beta_bar + eta*(Beta - Beta_bar) + berr*randn();
	if (Beta < 0)
		Beta = 0;

	// get derivatives
	float dS = - ( 1 - Phi*( (float) n/(n+1) ) )*Beta*S*I - ( (float) Phi/(n+1) )*S*(neiInf_sum);
	float dI = - dS - r*I;
	float dR = r*I;

	// step forward by h time units
	S += h*dS;
	I += h*dI;
	R += h*dR;

	population[cID].S = S;
	population[cID].I = I;
	population[cID].R = R;

	params[cID].Beta = Beta;

}

double dist (int x1, int y1, int x2, int y2) {

    return sqrt ( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

}