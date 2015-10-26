#include <iostream>
#include <cmath>
#include <fstream>
#include "rand.h"
#include "timer.h"

#define TS 		100			// time to simulate over
#define R0 		3.0			// infectiousness
#define r 		2e-1		// recovery rate
#define N 		500 		// population size
#define B 		R0*r/N		// transmission factor
#define pnoise 	1.0			// process noise case multiplier
#define merr 	N/100  		// measurement error, 1% of population
#define ALP 	1e-2		// cooling parameter for neighbour interactions
#define SpW		10			// number of Euler steps to take per week
#define eta		0.5			// Beta drift cooling
#define berr 	1e-3 		// Beta shocks

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
int timeval_subtract (double *result, struct timeval *x, struct timeval *y);
int check_float(float x,float y);
int poissrand(int lambda);
float randu();
float randn();
double dist (int x1, int y1, int x2, int y2);

using namespace std;

int main (int argc, char ** argv) {

	// Parse input arg(s) and check for errors ***************

	/*
	cout << "Got " << argc-1 << " argument(s)" << endl;
	for (int i = 1; i < argc; i++)
		cout << " [" << i << "] " << argv[i] << endl;
	*/

	int dim;
	if (argc == 2)
		dim = atoi(argv[1]);
	else {
		cout << "Missing dimension argument" << endl;
		return 1;
	}

	// *******************************************************


	// Parameters and such ***********************************

	int init_infected = 1;							// Number of initial infections in seed cells

	int T = max( (int) round(10*dim/R0), TS );			// Time to simulate over
	int nCells = dim*dim;							// Number of cells (a DIMxDIM square)

	SIR population[nCells];							// SIR compartment counts for each cell
	float nInfec_last[nCells];						// stores last infection counts between propagation steps
	int adjMat[nCells][8];							// Each row holds the indices of that cell's neighbours

	srand (time(NULL));								// Seed PRNG with system time

	int nInfecCells = max(dim/20, 1);				// Initial number of cells to infect
	int out_start_loc[nInfecCells];					// Locations to start infections

	float neiInf[nCells][8];						// used to store current neighbour infected counts

	// *******************************************************

	// Print simulation stats

	cout << "Simulation stats" << endl;
	cout << "----------------" << endl;
	cout << "Dims      " << dim << "x" << dim << endl;
	cout << "Time      " << T << endl;
	cout << "R0        " << R0 << endl;
	cout << "r (rec)   " << r << endl;


	// Generate Cell parameters ******************************

	Param params[nCells];
	for (int cell = 0; cell < nCells; cell++) {
		params[cell].Beta 	= B;
		params[cell].alpha 	= ALP;
		params[cell].Beta_bar = B;
	}

		/*
		// Circle of SLOW transmissibility
		int radius = round( (float) dim / 4.0 );
		for (int cell = 0; cell < nCells; cell++) {

			int i_cell = cell / dim;
			int j_cell = cell % dim;

			float distance = dist(i_cell,j_cell,dim/2,dim/2);

			if ( (radius-1) <  distance && distance < (radius+1)) {
				paramMat[cell].Beta 	= (float) B 	/ 100.0;
				paramMat[cell].alpha	= (float) ALP 	/ 100.0;
			}
		}
		*/

		// Line of SLOW transmissibility
		/*
		for (int cell = 0; cell < nCells; cell++) {

			int i_cell = cell / dim;
			int j_cell = cell % dim;

			float distance = dist(i_cell,j_cell,dim/2,dim/2);

			if ( j_cell == dim/2 ) {
				paramMat[cell].Beta 	= (float) B 	/ 2.0;
				paramMat[cell].alpha	= (float) ALP 	/ 2.0;
			}
		}
	*/

	// *******************************************************


	// populate adjacentcy matrix ****************************

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

	/*
	for (int cell = 0; cell < nCells; cell++) {
		cout << nNeibVec[cell] << " :\t";
		for (int i = 0; i < 8; i++)
			cout << adjMat[cell][i] << " ";
		cout << endl;
	}

	return 0;
	*/

	// *******************************************************


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


	// Data generation and output **************************

	ofstream file;
	file.open("ssir_data.dat");

	for (int t = 0; t < T; t++) {

		for (int st = 0; st < SpW; st++) {

			for (int cell = 0; cell < nCells; cell++) {
				nInfec_last[cell] = population[cell].I;
			}

			// propagate cells forward
			for (int cell = 0; cell < nCells; cell++)
				//exp_euler_SSIR( 1.0/SpW, 0.0, 1.0, &population[cell], &paramMat[cell], &neiInf[cell][0], nNeibVec[cell]);
				exp_euler_SSIR(1.0/SpW, population, params, nInfec_last, cell, adjMat, nNeibVec);

		}


		// perturb states (observation noise) and save as data for the week
		for (int cell = 0; cell < nCells; cell++) {

			int datapoint =  round(population[cell].I) + merr*randn();
			if(datapoint < 0)
				datapoint = 0;

			file << datapoint;
			if ( (cell-(dim-1)) % dim == 0)
				file << endl;
			else
				file << " ";

		}

	}

	file.close();

	// *******************************************************


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

	// update beta
	Beta = Beta_bar + eta * ( Beta - Beta_bar ) + berr * randn();

	// get derivatives
	float dS = - (1-n*alpha)*Beta*S*I - alpha*S*(neiInf_sum);
	float dI = - dS - r*I;
	float dR = r*I;

	// step forward by h time units
	S += h*dS;
	I += h*dI;
	R += h*dR;

	/*
	// Generate a number in [-pnoise,pnoise].
	// This will be noise representing more / less susceptibles getting sick
	float iShift = pnoise * poissrand( h*dI );
	if ( (I0+iShift) <= N && 0 <= (S0-iShift) ) {
		S = S0 - iShift;
		I = I0 + iShift;
	}
	*/

	population[cID].S = S;
	population[cID].I = I;
	population[cID].R = R;

	params[cID].Beta = Beta;

}

double dist (int x1, int y1, int x2, int y2) {

    return sqrt ( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

}