#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <random>

#include "rand.h"
#include "timer.h"
#include "readdata.h"

//#define TS 		100			// time to simulate over
//#define R0 		3.0			// infectiousness
//#define r 		2e-1		// recovery rate
#define N 		500 		// population size
//#define B 		R0*r/N		// transmission factor
//#define pnoise 	1.0			// process noise case multiplier
//#define merr 	N/100  		// measurement error, 1% of population
//#define ALP 	1e-2		// cooling parameter for neighbour interactions
#define SpW		10			// number of Euler steps to take per week
//#define eta		0.5			// Beta drift cooling
#define berr 	1e-3 		// Beta shocks
//#define Phi 	0.5			// Difficulty of transmission between cells, 0 = Impossible, 1 = no barrier

#define PI 		3.141592654

struct SIR {
	float S;
	float I;
	float R;
};

struct Param {
	float R0;
	float r;
	float sigma;
	float phi;
	SIR * states0;
};

void exp_euler_SSIR(float h, SIR * population, Param * params, float * nInfec_last, int cID, int neibIDs[][8], int * nNeibs);
double dist (int x1, int y1, int x2, int y2);

float prior(Param params, int * data, float nCells) {

	float pR0, pr, psigma, pphi, pI;
	float diff;
	float sum;

	diff = log(params.R0) 	 - 1;	pR0 	= log( 1.0/( params.R0	  *sqrt(2.0*PI) ) * expf( - diff*diff / 2 ) );
	diff = log(params.r)  	 - 1;	pr 		= log( 1.0/( params.r     *sqrt(2.0*PI) ) * expf( - diff*diff / 2 ) );
	diff = log(params.sigma) - 1;	psigma 	= log( 1.0/( params.sigma *sqrt(2.0*PI) ) * expf( - diff*diff / 2 ) );
	//diff = log(params.phi) 	 - 1;	pphi 	= log( 1.0/( params.phi   *sqrt(2.0*PI) ) * expf( - diff*diff / 2 ) );

	pI = 0.0;
	for (int cell = 0; cell < nCells; cell++) {
		diff = params.states0[cell].I - data[cell];
		pI += log( 1.0/( params.sigma *sqrt(2.0*PI) ) * expf( - diff*diff / (2*params.sigma*params.sigma) ) );
	}

	sum = pR0 + pr + psigma + pI;

	return sum;
	
}

float likelihood(float* data_hat, int * data, Param params, int npoints) {

	float sum = 0;
	float diff;

	float sigma = params.sigma;

	for (int p = 0; p < npoints; p++) {
		diff = data_hat[p] - data[p];
		sum += log( 1.0/(sigma*sqrt(2.0*PI)) * expf( - diff*diff / (2.0*sigma*sigma) ) );
	}

	return sum;

}

float posterior(float * data_hat, int * data, Param params, int nCells, int npoints) {

	return ( prior(params, data, nCells) + likelihood(data_hat, data, params, npoints) );

}

void proposal(Param * newparams, Param oldparams, int nCells) {

	float R0can, rcan, sigcan, phican, Scan, Ican, Rcan;

	R0can = oldparams.R0 + 0.5*randn();
	while(R0can < 0)
		R0can = oldparams.R0 + 0.5*randn();
	newparams -> R0 = R0can;

	rcan = oldparams.r + 0.05*randn();
	while(rcan < 0)
		rcan = oldparams.r + 0.05*randn();
	newparams -> r = rcan;

	sigcan = oldparams.sigma + 0.5*randn();
	while(sigcan < 0)
		sigcan = oldparams.sigma + 0.5*randn();
	newparams -> sigma 	= sigcan;

	phican = randu();
	newparams -> phi = phican;

	for (int cell = 0; cell < nCells; cell++) {

		Ican = oldparams.states0[cell].I + 0.1*randn();
		while(Ican < 0)
			Ican = oldparams.states0[cell].I + 0.1*randn();

		newparams->states0[cell].S = N - Ican;
		newparams->states0[cell].I = Ican;
		newparams->states0[cell].R = 0.0;

	}

}

using namespace std;

int main (int argc, char ** argv) {

	int dim, ydim;

	double restime;
	struct timeval tdr0, tdr1, tdrMaster;

	// Parse arguments **********************************************

	if (argc < 3) {
		std::cout << "Not enough arguments" << std::endl;
		return 0;
	}

	std::string arg1(argv[1]);
	std::string arg2(argv[2]);

	std::cout << "Arguments:" << std::endl;
	std::cout << "    [1] " << arg1 << std::endl;
	std::cout << "    [2] " << arg2 << std::endl;

	// **************************************************************


	// Read data ****************************************************

	std::cout << "Getting data" << std::endl;

	int * trueCounts 	= getData(arg2, &dim, &ydim);
	int * data 			= getData(arg1, NULL, NULL);

	// **************************************************************

	int T = ydim / dim;
	int nCells = dim*dim;

	// Parameters and such **********************************************************************************

	SIR 	population[nCells];									// SIR compartment counts for each cell
	float 	nInfec_last[nCells];								// stores last infection counts between propagation steps
	int 	adjMat[nCells][8];									// Each row holds the indices of that cell's neighbours
	float 	neiInf[nCells][8];									// used to store current neighbour infected counts

	float 	truestates[T*nCells];								// generated system "true" states
	float  	datastates[T*nCells]; 								// perturbed generated "true" states

	Param 	params, params_old, params_new;						// all system parameters
	params.states0 	= (SIR*) malloc (nCells*sizeof(SIR)); 		// allocate space for initial SIR states
	params_old.states0 = (SIR*) malloc (nCells*sizeof(SIR)); 	// allocate space for initial SIR states
	params_new.states0 = (SIR*) malloc (nCells*sizeof(SIR)); 	// allocate space for initial SIR states

	srand(time(NULL));											// Seed PRNG with system time

	// ******************************************************************************************************


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


	// Data generation and output ***************************************************************

	int numsamples 	= 0;
	int BURNMAX 	= 1000;
	int SAMPLEMAX 	= 1000;

	Param 	chain[SAMPLEMAX];
	Param 	params_last[nCells];

	// seed values

	// PARAM SET 1

	params.R0 		= 3.0;
	params.r 		= 0.2;
	params.sigma 	= 5.0;
	params.phi 		= 0.5;

	int refcell = 45; 	// cell with index 45 has the first case
	for (int cell = 0; cell < nCells; cell++) {
		params.states0[cell].S = 500;
		params.states0[cell].I = 0;
		params.states0[cell].R = 0;
	}
	params.states0[refcell].S = 499;
	params.states0[refcell].I = 1;

	int nacc 	= 0;
	int nrej 	= 0;
	int MAXACC 	= 10;

	float posterior_old, posterior_new;

	// INITIALIZATION

	// initialize system with parameter initial state estimates
	for (int cell = 0; cell < nCells; cell++) {
		population[cell].S = params.states0[cell].S;
		population[cell].I = params.states0[cell].I;
		population[cell].R = params.states0[cell].R;
	}

	// Evolve

	for (int t = 0; t < T; t++) {

		for (int st = 0; st < SpW; st++) {

			for (int cell = 0; cell < nCells; cell++) {
				nInfec_last[cell] = population[cell].I;
			}

			// propagate cells forward
			for (int cell = 0; cell < nCells; cell++)
				exp_euler_SSIR(1.0/SpW, population, &params, nInfec_last, cell, adjMat, nNeibVec);

		}

		// perturb states (observation noise) and save as data for the week
		for (int cell = 0; cell < nCells; cell++) {

			int truepoint = round( population[cell].I );
			int datapoint = round( truepoint + params.sigma*randn() );
			if(datapoint < 0)
				datapoint = 0;

			truestates[t*nCells+cell] = truepoint;
			datastates[t*nCells+cell] = datapoint;

		}

	}

	params_old.R0 		= params.R0;
	params_old.r 		= params.r;
	params_old.sigma 	= params.sigma;
	params_old.phi 		= params.phi;
	params_old.R0 		= params.R0;

	for (int cell = 0; cell < nCells; cell++) {

		params_old.states0[cell].S = params.states0[cell].S;
		params_old.states0[cell].I = params.states0[cell].I;
		params_old.states0[cell].R = params.states0[cell].R;

	}

	posterior_old = posterior(truestates, data, params_old, nCells, nCells*T);

	//cout << "R0 \t r \t sigma \t posterior" << endl;
	//cout << "-- \t - \t ----- \t ---------" << endl;

	int ctr = 0;

	while (nacc < MAXACC) {

		proposal(&params_new, params_old, nCells);

		// initialize system with parameter initial state estimates
		for (int cell = 0; cell < nCells; cell++) {
			population[cell].S = params_new.states0[cell].S;
			population[cell].I = params_new.states0[cell].I;
			population[cell].R = params_new.states0[cell].R;
		}

		// Evolve

		for (int t = 0; t < T; t++) {

			for (int st = 0; st < SpW; st++) {

				for (int cell = 0; cell < nCells; cell++) {
					nInfec_last[cell] = population[cell].I;
				}

				// propagate cells forward
				for (int cell = 0; cell < nCells; cell++)
					exp_euler_SSIR(1.0/SpW, population, &params_new, nInfec_last, cell, adjMat, nNeibVec);

			}

			// perturb states (observation noise) and save as data for the week
			for (int cell = 0; cell < nCells; cell++) {

				int truepoint = round( population[cell].I );
				int datapoint = round( truepoint + params_new.sigma*randn() );
				if(datapoint < 0)
					datapoint = 0;

				truestates[t*nCells+cell] = truepoint;
				datastates[t*nCells+cell] = datapoint;

			}

		}

		float posterior_new = posterior(truestates, data, params_new, nCells, nCells*T);

		float P = expf(posterior_new - posterior_old);

		if ( randu() < P ) {

			// step is accepted

			params_old.R0 		= params_new.R0;
			params_old.r 		= params_new.r;
			params_old.sigma 	= params_new.sigma;
			params_old.phi 		= params_new.phi;
			params_old.R0 		= params_new.R0;

			for (int cell = 0; cell < nCells; cell++) {

				params_old.states0[cell].S = params_new.states0[cell].S;
				params_old.states0[cell].I = params_new.states0[cell].I;
				params_old.states0[cell].R = params_new.states0[cell].R;

			}

			posterior_old = posterior_new;

			nacc++;

		} else {

			// step is rejected

			nrej++;

		}

		ctr++;

		if ( (ctr%100) == 0 ) {
			cout << ctr << " Acc " << nacc << " Nrej " << nrej << endl;
		}

	}

	cout << "Accepted: " << nacc << endl;
	cout << "Rejected: " << nrej << endl;


	//for (int t = 0; t < T; t++)
	//	cout << t << ":\t" << trueCounts[t*nCells+refcell] << "\t" << truestates[t][refcell] << endl;

	//float lik_test = likelihood(float* data_hat, float * data, Param params, int npoints)
	
		


	


	// LATER

	// Initial parameter values

		// Proposal distribution sampling

		// Initialize population states to parameter values

		// EVOLVE SYSTEM

			// perturb states (observation noise) and save as data for the week

		// COMPUTE LIKELIHOOD

	// ******************************************************************************************

	return 0;

}


void exp_euler_SSIR(float h, SIR * population, Param * params, float * nInfec_last, int cID, int neibIDs[][8], int * nNeibs) {

	float R0 		= params->R0;
	float r 		= params->r;
	float Beta 		= R0*r / N;
	float phi 		= params->phi;

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
		neiInf_sum += Beta * nInfec_last[nID];
	}

	// get derivatives
	float dS = - ( 1 - phi*( (float) n/(n+1) ) )*Beta*S*I - ( (float) phi/(n+1) )*S*(neiInf_sum);
	float dI = - dS - r*I;
	float dR = r*I;

	// step forward by h time units
	S += h*dS;
	I += h*dI;
	R += h*dR;

	population[cID].S = S;
	population[cID].I = I;
	population[cID].R = R;

}

double dist (int x1, int y1, int x2, int y2) {

    return sqrt ( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );

}