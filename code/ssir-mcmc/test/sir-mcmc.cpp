#include <iostream>
#include <cmath>
#include <string>
#include <random>

#include "timer.h"
#include "rand.h"

#define N 	500
#define PI 	3.141592654

struct State {
	float S;
	float I;
	float R;
};

struct Param {
	float R0;
	float r;
	float sigma;
	float S;
	float I;
	float R;
};

void sir_integrate(State * states, Param params, int nsteps, int T) {

	float R0 	= params.R0;
	float r 	= params.r;
	float Beta 	= R0*r / N;

	float S = params.S;
	float I = params.I;
	float R = params.R;

	float h = 1.0/nsteps;

	states[0].S = S;
	states[0].I = I;
	states[0].R = R;

	for (int t = 1; t < T; t++) {

		for (int s = 0; s < nsteps; s++) {

			S += h*( - S*I*Beta );
			I += h*( S*I*Beta - r*I );
			R += h*( r*I );

		}

		states[t].S = S;
		states[t].I = I;
		states[t].R = R;

	}

}

float prior(Param params) {

	float pR0, pr, psigma, pI;
	float diff;
	float sum;

	diff = log(params.R0) 	 - 1;	pR0 	= log( 1.0/( params.R0	  *sqrt(2.0*PI) ) * expf( - diff*diff / 2 ) ); 	//R0
	diff = log(params.r)  	 - 1;	pr 		= log( 1.0/( params.r     *sqrt(2.0*PI) ) * expf( - diff*diff / 2 ) ); 	//r
	diff = log(params.sigma) - 1;	psigma 	= log( 1.0/( params.sigma *sqrt(2.0*PI) ) * expf( - diff*diff / 2 ) ); 	//sigma
	diff = params.I 	     - 5; 	pI 		= log( 1.0/( params.sigma *sqrt(2.0*PI) ) *
																expf( - diff*diff / (2*params.sigma*params.sigma) ) );	//I
	
	sum = pR0 + pr + psigma + pI;

	return sum;
	
}

float likelihood(State * states, float * data, int npoints, Param params) {

	float sum = 0;
	float diff;

	float sigma = params.sigma;

	for (int p = 0; p < npoints; p++) {
		diff = states[p].I - data[p];
		sum += log( 1.0/(sigma*sqrt(2.0*PI)) * expf( - diff*diff / (2.0*sigma*sigma) ) );
	}

	return sum;

}

float posterior(State * states, float * data, int npoints, Param params) {

	return ( prior(params) + likelihood(states, data, npoints, params) );

}

void proposal(Param * newparams, Param oldparams) {

	float R0can, rcan, sigcan, Scan, Ican, Rcan;

	float litrat = 0.045;
	float bigrat = 0.45;

	R0can = oldparams.R0 + bigrat*randn();
	while(R0can < 0)
		R0can = oldparams.R0 + bigrat*randn();

	rcan = oldparams.r + litrat*randn();
	while(rcan < 0)
		rcan = oldparams.r + litrat*randn();

	sigcan = oldparams.sigma + bigrat*randn();
	while(sigcan < 0)
		sigcan = oldparams.sigma + bigrat*randn();

	Ican = oldparams.I + bigrat*randn();
	while(Ican < 0)
		Ican = oldparams.I + bigrat*randn();

	Scan = N - Ican;
	Rcan = 0;

	newparams -> R0 	= R0can;
	newparams -> r  	= rcan;
	newparams -> sigma 	= sigcan;
	newparams -> S 		= Scan;
	newparams -> I 		= Ican;
	newparams -> R 		= Rcan;

}

void printParams(Param params) {

	std::cout << "R0\tr\t\tsigma\tS\tI\tR" << std::endl;
	std::cout << "--\t-\t\t-----\t-\t-\t-" << std::endl;
	std::cout 	<< params.R0 << "\t"
				<< params.r << "\t\t"
				<< params.sigma << "\t"
				<< params.S << "\t"
				<< params.I << "\t"
				<< params.R << std::endl;

}

int main (int argc, char ** argv) {

	srand(time(NULL));								// Seed PRNG with system time

	int dataflag = 0;
	int diagflag = 0;

	int BURNMAX 	= 1000000;
	int SAMPLEMAX 	= 200000;

	if (argc == 2) {

		std::string mode = argv[1];

		if (mode == "data") {
			dataflag = 1;
		} else if (mode == "diagnostic") {
			diagflag = 1;
		} else {
			std::cout << "Unrecognized mode argument, needs to be [data|diagnostic]" << std:: endl; 
			return 0;
		}

	} else {

		std::cout << "Argument requirements: either no argument or [data|diagnostic]" << std:: endl; 
		return 0;

	}

	int T 		= 100;
	int nsteps 	= 10;

	Param params_real = {3.0, 0.1, 5, 495, 5, 0};
	Param params_new, params_old;

	State states_real[T];
	State states[T];

	Param chain[SAMPLEMAX];

	// TRUE states
	sir_integrate(states_real, params_real, nsteps, T);

	// DATA states (fake)
	float data[T];
	float datapoint;
	for (int t = 0; t < T; t++) {
		datapoint = states_real[t].I + 5*randn();
		while (datapoint < 0)
			datapoint = states_real[t].I + 5*randn();
		data[t] = datapoint;
	}

	// MCMC loop

	int numsamples = 0;
	int numrej = 0;
	int numacc = 0;

	// starting parameters
	params_old.R0 		= 1.0;	//params_real.R0;
	params_old.r 		= 0.3;	//params_real.r;
	params_old.sigma 	= 3.0; 	//params_real.sigma;
	params_old.S 		= 497.0;	//params_real.S;
	params_old.I 		= 3.0;	//params_real.I;
	params_old.R 		= 0.0;	//params_real.R;

	// starting states and parameters
	sir_integrate(states, params_old, nsteps, T);
	float posterior_old = posterior(states, data, T, params_old);


	while( numsamples < (BURNMAX + SAMPLEMAX) ) {

		proposal(&params_new, params_old);
		sir_integrate(states, params_new, nsteps, T);

		float posterior_new = posterior(states, data, T, params_new);

		float P = expf( posterior_new - posterior_old );

		if ( randu() < P ) {

			// Step is accepted

			params_old.R0 		= params_new.R0;
			params_old.r 		= params_new.r;
			params_old.sigma 	= params_new.sigma;
			params_old.S 		= params_new.S;
			params_old.I 		= params_new.I;
			params_old.R 		= params_new.R;

			posterior_old = posterior_new;

			numacc++;

		} else {

			// Step is rejected
			numrej++;

		}

		if (BURNMAX <= numsamples) {
			chain[numsamples - BURNMAX].R0 		= params_old.R0;
			chain[numsamples - BURNMAX].r 		= params_old.r;
			chain[numsamples - BURNMAX].sigma 	= params_old.sigma;
			chain[numsamples - BURNMAX].S 		= params_old.S;
			chain[numsamples - BURNMAX].I 		= params_old.I;
			chain[numsamples - BURNMAX].R 		= params_old.R;
		}

		numsamples++;

	}


	if (diagflag) {

		float meanR0 	= 0;
		float meanr 	= 0;
		float meansig 	= 0;
		float meanS 	= 0;
		float meanI 	= 0;
		float meanR 	= 0;

		for (int i = 0; i < SAMPLEMAX; i++) {

			meanR0 	+= chain[i].R0  	/ SAMPLEMAX;
			meanr 	+= chain[i].r   	/ SAMPLEMAX;
			meansig += chain[i].sigma 	/ SAMPLEMAX;
			meanS 	+= chain[i].S   	/ SAMPLEMAX;
			meanI 	+= chain[i].I   	/ SAMPLEMAX;
			meanR 	+= chain[i].R   	/ SAMPLEMAX;

		}

		std::cout << "Steps accepted " << numacc << std::endl;
		std::cout << "Steps rejected " << numrej << std::endl;

		std::cout << "R0   " << meanR0 	<< std::endl;
		std::cout << "r    " << meanr 	<< std::endl;
		std::cout << "sig  " << meansig << std::endl;
		std::cout << "S    " << meanS 	<< std::endl;
		std::cout << "I    " << meanI 	<< std::endl;
		std::cout << "R    " << meanR 	<< std::endl;

	}

	if (dataflag) {

		for (int i = 0; i < SAMPLEMAX; i++) 
			std::cout 	<< chain[i].R0 		<< " "
						<< chain[i].r 		<< " "
						<< chain[i].sigma 	<< " "
						<< chain[i].S 		<< " "
						<< chain[i].I 		<< " "
						<< chain[i].R 		<< std::endl;

	}



}

