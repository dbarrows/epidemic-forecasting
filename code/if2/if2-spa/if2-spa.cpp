/*	Author: Dexter Barrows
	Github: dbarrows.github.io

	*/

#include <stdio.h> 
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>

//#include "rand.h"
//#include "timer.h"

#define Treal 		100			// time to simulate over
#define R0true 		3.0			// infectiousness
#define rtrue 		0.1			// recovery rate
#define Nreal 		500.0 		// population size
#define etatrue 	0.5			// real drift attraction strength
#define berrtrue	0.5			// real beta drift noise
#define phitrue 	0.5 		// real connectivity strength
#define merr 		10.0  		// expected measurement error
#define I0 			5.0			// Initial infected individuals

#define PSC  		0.5 		// perturbation scale factor for more sensitive parameters

#include <Rcpp.h>
using namespace Rcpp;

struct Particle {
	double R0;
	double r;
	double sigma;
	double eta;
	double berr;
	double phi;
	double * S;
	double * I;
	double * R;
	double * B;
	double * Iinit;
};


int timeval_subtract (double *result, struct timeval *x, struct timeval *y);
int check_double(double x,double y);
void initializeParticles(Particle ** particles, int NP, int nloc, int N);
void exp_euler_SSIR(double h, double t0, double tn, int N, Particle * particle,
                    NumericVector neinum, NumericMatrix neibmat, int nloc) ;
void copyParticle(Particle * dst, Particle * src, int nloc);
void perturbParticles(Particle * particles, int N, int NP, int nloc, int passnum, double coolrate);
double randu();
double randn();

// [[Rcpp::export]]
Rcpp::List if2_spa(NumericMatrix data, int T, int N, int NP, int nPasses, double coolrate, NumericVector neinum, NumericMatrix neibmat, int nloc) {

	NumericMatrix paramdata(NP, 6); 	// for R0, r, sigma, eta, berr, phi
	NumericMatrix initInfec(nloc, NP); 	// for Iinit
	NumericMatrix infecmeans(nloc, T); 	// mean infection counts for each location
	NumericMatrix finalstate(nloc, 4); 	// SIRB means for each location

	srand(time(NULL));		// Seed PRNG with system time

	double w[NP]; 			// particle weights

	// initialize particles
	printf("Initializing particle states\n");
	Particle * particles = NULL; 		// particle estimates for current step
	Particle * particles_old = NULL; 	// intermediate particle states for resampling
	initializeParticles(&particles, NP, nloc, N);
	initializeParticles(&particles_old, NP, nloc, N);

	/*
	// copy particle test
	copyParticle(&particles[0], &particles_old[0], nloc);

	// perturb particle test
	perturbParticles(particles, N, NP, nloc, 1, coolrate);

	// evolution test
	// reset particle system evolution states
	for (int n = 0; n < NP; n++) {
		for (int loc = 0; loc < nloc; loc++) {
			particles[n].S[loc] = N - particles[n].Iinit[loc];
			particles[n].I[loc] = particles[n].Iinit[loc];
			particles[n].R[loc] = 0.0;
			particles[n].B[loc] = (double) particles[n].R0 * particles[n].r / N;
		}
	}
	printf("Before S:%f | I:%f | R:%f\n", particles[0].S[0], particles[0].I[0], particles[0].R[0]);
	exp_euler_SSIR(1.0/7.0, 0.0, 1.0, N, &particles[0], neinum, neibmat, nloc);
	printf("After S:%f | I:%f | R:%f\n", particles[0].S[0], particles[0].I[0], particles[0].R[0]);
	*/

	// START PASSES THROUGH DATA

	printf("Starting filter\n");
	printf("---------------\n");
	printf("Pass\n");


	for (int pass = 0; pass < nPasses; pass++) {

		printf("...%d / %d\n", pass, nPasses);

		// reset particle system evolution states
		for (int n = 0; n < NP; n++) {
			for (int loc = 0; loc < nloc; loc++) {
				particles[n].S[loc] = N - particles[n].Iinit[loc];
				particles[n].I[loc] = particles[n].Iinit[loc];
				particles[n].R[loc] = 0.0;
				particles[n].B[loc] = (double) particles[n].R0 * particles[n].r / N;
			}
		}

		if (pass == (nPasses-1)) {
			double means[nloc];
			for (int loc = 0; loc < nloc; loc++) {
				means[loc] = 0.0;
				for (int n = 0; n < NP; n++) {
					means[loc] += particles[n].I[loc] / NP;	
				}
				infecmeans(loc, 0) = means[loc];
			}
		}

		for (int t = 1; t < T; t++) {

			// generate individual predictions and weight
			for (int n = 0; n < NP; n++) {

				exp_euler_SSIR(1.0/7.0, 0.0, 1.0, N, &particles[n], neinum, neibmat, nloc);

				double merr_par = particles[n].sigma;

				w[n] = 1.0;
				for (int loc = 0; loc < nloc; loc++) {
					double y_diff 	= data(loc, t) - particles[n].I[loc];
					w[n] *= 1.0/(merr_par*sqrt(2.0*M_PI)) * exp( - y_diff*y_diff / (2.0*merr_par*merr_par) );
				}

			}

			// cumulative sum
			for (int n = 1; n < NP; n++) {
				w[n] += w[n-1];
			}

			// save particle states to resample from
			for (int n = 0; n < NP; n++){
				copyParticle(&particles_old[n], &particles[n], nloc);
			}

			// resampling
			for (int n = 0; n < NP; n++) {

				double w_r = randu() * w[NP-1];
				int i = 0;
				while (w_r > w[i]) {
					i++;
				}

				// i is now the index to copy state from
				copyParticle(&particles[n], &particles_old[i], nloc);

			}
			
			// between-iteration perturbations, not after last time step
			if (t < (T-1))
			    perturbParticles(particles, N, NP, nloc, pass, coolrate);

			if (pass == (nPasses-1)) {
				double means[nloc];
				for (int loc = 0; loc < nloc; loc++) {
					means[loc] = 0.0;
					for (int n = 0; n < NP; n++) {
						means[loc] += particles[n].I[loc] / NP;	
					}
					infecmeans(loc, t) = means[loc];
				}
			}

		}

		// between-pass perturbations, not after last pass
		if (pass < (nPasses + 1))
		    perturbParticles(particles, N, NP, nloc, pass, coolrate);

	}

	// pack parameter data (minus initial conditions)
	for (int n = 0; n < NP; n++) {
		paramdata(n, 0) = particles[n].R0;
		paramdata(n, 1) = particles[n].r;
		paramdata(n, 2) = particles[n].sigma;
		paramdata(n, 3) = particles[n].eta;
		paramdata(n, 4) = particles[n].berr;
		paramdata(n, 5) = particles[n].phi;
	}

	// Pack initial condition data
	for (int n = 0; n < NP; n++) {
		for (int loc = 0; loc < nloc; loc++) {
			initInfec(loc, n) = particles[n].Iinit[loc];
		}
	}

	// Pack final state means data
	double Smeans[nloc], Imeans[nloc], Rmeans[nloc], Bmeans[nloc];
	for (int loc = 0; loc < nloc; loc++) {
		Smeans[loc] = 0.0;
		Imeans[loc] = 0.0;
		Rmeans[loc] = 0.0;
		Bmeans[loc] = 0.0;
		for (int n = 0; n < NP; n++) {
			Smeans[loc] += particles[n].S[loc] / NP;
			Imeans[loc] += particles[n].I[loc] / NP;
			Rmeans[loc] += particles[n].R[loc] / NP;
			Bmeans[loc] += particles[n].B[loc] / NP;
		}
		finalstate(loc, 0) = Smeans[loc];
		finalstate(loc, 1) = Imeans[loc];
		finalstate(loc, 2) = Rmeans[loc];
		finalstate(loc, 3) = Bmeans[loc];
	}


	return Rcpp::List::create(	Rcpp::Named("paramdata") = paramdata,
	                          	Rcpp::Named("initInfec") = initInfec,
                             	Rcpp::Named("infecmeans") = infecmeans,
                             	Rcpp::Named("finalstate") = finalstate);



}


/*	Use the Explicit Euler integration scheme to integrate SIR model forward in time
	double h 	- time step size
	double t0 	- start time
	double tn 	- stop time
	double * y 	- current system state; a three-component vector representing [S I R], susceptible-infected-recovered

	*/
void exp_euler_SSIR(double h, double t0, double tn, int N, Particle * particle,
                    NumericVector neinum, NumericMatrix neibmat, int nloc) {

	int num_steps = floor( (tn-t0) / h );

	double * S = particle->S;
	double * I = particle->I;
	double * R = particle->R;
	double * B = particle->B;

	// create last state vectors
	double S_last[nloc];
	double I_last[nloc];
	double R_last[nloc];
	double B_last[nloc];

	double R0 	= particle->R0;
	double r 	= particle->r;
	double B0 	= R0 * r / N;
	double eta 	= particle->eta;
	double berr = particle->berr;
	double phi  = particle->phi;

	//printf("sphi \t\t| ophi \t\t| BSI \t\t| rI \t\t| dS \t\t| dI \t\t| dR \t\t| S \t\t| I \t\t| R |\n");

	for(int t = 0; t < num_steps; t++) {

		for (int loc = 0; loc < nloc; loc++) {
			S_last[loc] = S[loc];
			I_last[loc] = I[loc];
			R_last[loc] = R[loc];
			B_last[loc] = B[loc];
		}

		for (int loc = 0; loc < nloc; loc++) {

			B[loc] = exp( log(B_last[loc]) + eta*(log(B0) - log(B_last[loc])) + berr*randn() );

			int n = neinum[loc];
        	double sphi = 1.0 - phi*( (double) n/(n+1.0) );
        	double ophi = phi/(n+1.0);

        	double nBIsum = 0.0;
        	for (int j = 0; j < n; j++)
        		nBIsum += B_last[(int) neibmat(loc, j) - 1] * I_last[(int) neibmat(loc, j) - 1];

        	double BSI = S_last[loc]*( sphi*B_last[loc]*I_last[loc] + ophi*nBIsum );
        	double rI  = r*I_last[loc];

			// get derivatives
			double dS = - BSI;
			double dI = BSI - rI;
			double dR = rI;

			// step forward by h
			S[loc] += h*dS;
			I[loc] += h*dI;
			R[loc] += h*dR;

			//if (loc == 1)
			//	printf("%f\t|%f\t|%f\t|%f\t|%f\t|%f\t|%f\t|%f\t|%f\t|%f\t|\n", sphi, ophi, BSI, rI, dS, dI, dR, S[1], I[1], R[1]);

		}

	}

	/*particle->S = S;
	particle->I = I;
	particle->R = R;
	particle->B = B;*/

}

/* 	Initializes particles
	*/
void initializeParticles(Particle ** particles, int NP, int nloc, int N) {

	// allocate space for doubles
	*particles = (Particle*) malloc (NP*sizeof(Particle));

	// allocate space for arays inside particles
	for (int n = 0; n < NP; n++) {
		(*particles)[n].S = (double*) malloc(nloc*sizeof(double));
		(*particles)[n].I = (double*) malloc(nloc*sizeof(double));
		(*particles)[n].R = (double*) malloc(nloc*sizeof(double));
		(*particles)[n].B = (double*) malloc(nloc*sizeof(double));
		(*particles)[n].Iinit = (double*) malloc(nloc*sizeof(double));
	}

	// initialize all all parameters
	for (int n = 0; n < NP; n++) {

		double R0can, rcan, sigmacan, Iinitcan, etacan, berrcan, phican;

		do {
			R0can = R0true + R0true*randn();
		} while (R0can < 0);
		(*particles)[n].R0 = R0can;

		do {
			rcan = rtrue + rtrue*randn();
		} while (rcan < 0);
		(*particles)[n].r = rcan;

		for (int loc = 0; loc < nloc; loc++)
			(*particles)[n].B[loc] = (double) R0can * rcan / N;

		do {
			sigmacan = merr + merr*randn();
		} while (sigmacan < 0);
		(*particles)[n].sigma = sigmacan;

		do {
			etacan = etatrue + PSC*etatrue*randn();
		} while (etacan < 0 || etacan > 1);
		(*particles)[n].eta = etacan;

		do {
			berrcan = berrtrue + PSC*berrtrue*randn();
		} while (berrcan < 0);
		(*particles)[n].berr = berrcan;

		do {
			phican = phitrue + PSC*phitrue*randn();
		} while (phican <= 0 || phican >= 1);
		(*particles)[n].phi = phican;

		for (int loc = 0; loc < nloc; loc++) {
			do {
				Iinitcan = I0 + I0*randn();
			} while (Iinitcan < 0 || N < Iinitcan);
			(*particles)[n].Iinit[loc] = Iinitcan;
		}

	}

}

/*	Particle pertubation function to be run between iterations and passes

	*/
void perturbParticles(Particle * particles, int N, int NP, int nloc, int passnum, double coolrate) {

	//double coolcoef = exp( - (double) passnum / coolrate );
	double coolcoef = pow(coolrate, passnum);
	
    double spreadR0 	= coolcoef * R0true / 10.0;
    double spreadr 		= coolcoef * rtrue / 10.0;
    double spreadsigma 	= coolcoef * merr / 10.0;
    double spreadIinit 	= coolcoef * I0 / 10.0;
    double spreadeta 	= coolcoef * etatrue / 10.0;
    double spreadberr 	= coolcoef * berrtrue / 10.0;
    double spreadphi 	= coolcoef * phitrue / 10.0;

    double R0can, rcan, sigmacan, Iinitcan, etacan, berrcan, phican;

    for (int n = 0; n < NP; n++) {

    	do {
    		R0can = particles[n].R0 + spreadR0*randn();
    	} while (R0can < 0);
    	particles[n].R0 = R0can;

    	do {
    		rcan = particles[n].r + spreadr*randn();
    	} while (rcan < 0);
    	particles[n].r = rcan;

    	do {
    		sigmacan = particles[n].sigma + spreadsigma*randn();
    	} while (sigmacan < 0);
    	particles[n].sigma = sigmacan;

    	do {
    		etacan = particles[n].eta + PSC*spreadeta*randn();
    	} while (etacan < 0 || etacan > 1);
    	particles[n].eta = etacan;

    	do {
    		berrcan = particles[n].berr + PSC*spreadberr*randn();
    	} while (berrcan < 0);
    	particles[n].berr = berrcan;

    	do {
    		phican = particles[n].phi + PSC*spreadphi*randn();
    	} while (phican <= 0 || phican >= 1);
    	particles[n].phi = phican;

    	for (int loc = 0; loc < nloc; loc++) {
	    	do {
	    		Iinitcan = particles[n].Iinit[loc] + spreadIinit*randn();
	    	} while (Iinitcan < 0 || Iinitcan > 500);
	    	particles[n].Iinit[loc] = Iinitcan;
	    }
    }

}

/*	Convinience function for particle resampling process
	*/
void copyParticle(Particle * dst, Particle * src, int nloc) {

	dst->R0 	= src->R0;
	dst->r 		= src->r;
	dst->sigma 	= src->sigma;
	dst->eta 	= src->eta;
	dst->berr 	= src->berr;
	dst->phi 	= src->phi;

	for (int n = 0; n < nloc; n++) {
		dst->S[n]		= src->S[n];
		dst->I[n]		= src->I[n];
		dst->R[n] 		= src->R[n];
		dst->B[n]      	= src->B[n];
		dst->Iinit[n]  	= src->Iinit[n];
	}

}



double randu() {

	return (double) rand() / (double) RAND_MAX;

}

/*
void getStateMeans(State * state, Particle* particles, int NP) {

	double Smean = 0, Imean = 0, Rmean = 0;

	for (int n = 0; n < NP; n++) {
		Smean += particles[n].S;
		Imean += particles[n].I;
		Rmean += particles[n].R;
	}

	state->S = (double) Smean / NP;
	state->I = (double) Imean / NP;
	state->R = (double) Rmean / NP;

}
*/

/*	Return a normally distributed random number with mean 0 and standard deviation 1
	Uses the polar form of the Box-Muller transformation
	From http://www.design.caltech.edu/erik/Misc/Gaussian.html
	*/
double randn() {

	double x1, x2, w, y1;

	do {
		x1 = 2.0 * randu() - 1.0;
		x2 = 2.0 * randu() - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;

	return y1;

}