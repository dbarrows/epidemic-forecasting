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

/*	True parameter values included for scaling purposes
	In reality these would have to be guessed to an order of magnitude or so
	*/
#define Treal 		100			// time to simulate over
#define R0true 		3.0			// infectiousness
#define rtrue 		0.1			// recovery rate
#define Nreal 		500.0 		// population size
#define etatrue 	0.5			// real drift attraction strength
#define berrtrue	0.5			// real beta drift noise
#define merr 		10.0  		// expected measurement error
#define I0 			5.0			// Initial infected individuals

#define PSC  		0.5 		// scale factor for more sensitive parameters

#include <Rcpp.h>
//using namespace Rcpp;

struct State {
	double S;
	double I;
	double R;
};

struct Particle {
	double R0;
	double r;
	double sigma;
	double eta;
	double berr;
	double B;
	double S;
	double I;
	double R;
	double Sinit;
	double Iinit;
	double Rinit;
};

struct ParticleInfo {
	double R0mean;		double R0sd;
	double rmean;		double rsd;
	double sigmamean;	double sigmasd;
	double etamean;		double etasd;
	double berrmean;	double berrsd;
	double Sinitmean; 	double Sinitsd;
	double Iinitmean;	double Iinitsd;
	double Rinitmean;	double Rinitsd;
};


int timeval_subtract (double *result, struct timeval *x, struct timeval *y);
int check_double(double x,double y);
void exp_euler_SIR(double h, double t0, double tn, int N, Particle * particle);
void copyParticle(Particle * dst, Particle * src);
void perturbParticles(Particle * particles, int N, int NP, int passnum, double coolrate);
void particleDiagnostics(ParticleInfo * partInfo, Particle * particles, int NP);
void getStateMeans(State * state, Particle* particles, int NP);
double randu();
double randn();

// [[Rcpp::export]]
Rcpp::List if2_s(Rcpp::NumericVector data, int T, int N, int NP, double coolrate, Rcpp::NumericVector params) {

	// params will be in the order R0, r, sigma, eta, berr, Sinit, Iinit, Rinit

	Rcpp::NumericMatrix statemeans(T, 3);
	Rcpp::NumericMatrix statedata(NP, 4);

	srand(time(NULL));		// Seed PRNG with system time

	double w[NP]; 			// particle weights

	Particle particles[NP]; 	// particle estimates for current step
	Particle particles_old[NP]; // intermediate particle states for resampling

	printf("Initializing particle states\n");

	// Initialize particle parameter states (seeding)
	// Note that they will all be the same for this code as we are only filtering over the states
	for (int n = 0; n < NP; n++) {

		particles[n].R0 	= params[0];
		particles[n].r 		= params[1];
		particles[n].sigma 	= params[2];
		particles[n].eta 	= params[3];
		particles[n].berr 	= params[4];
		particles[n].Sinit 	= params[5];
		particles[n].Iinit 	= params[6];
		particles[n].Rinit 	= params[7];

	}

	// START PASS THROUGH DATA

	printf("Starting filter\n");
	printf("---------------\n");

	// seed initial states
	for (int n = 0; n < NP; n++) {

		double Iinitcan;

		double i_infec = particles[n].Iinit;

		do {
			Iinitcan = i_infec + i_infec*randn();
		} while (Iinitcan < 0 || N < Iinitcan);

		particles[n].S = N - Iinitcan;
		particles[n].I = Iinitcan;
		particles[n].R = 0.0;

		particles[n].B = (double) particles[n].R0 * particles[n].r / N;

	}

	State sMeans;
	getStateMeans(&sMeans, particles, NP);
	statemeans(0,0) = sMeans.S;
	statemeans(0,1) = sMeans.I;
	statemeans(0,2) = sMeans.R;

	// start run through data
	for (int t = 1; t < T; t++) {

		// generate individual predictions and weight
		for (int n = 0; n < NP; n++) {

			exp_euler_SIR(1.0/7.0, 0.0, 1.0, N, &particles[n]);

			double merr_par = particles[n].sigma;
			double y_diff 	= data[t] - particles[n].I;
			
			w[n] = 1.0/(merr_par*sqrt(2.0*M_PI)) * exp( - y_diff*y_diff / (2.0*merr_par*merr_par) );

		}

		// cumulative sum
		for (int n = 1; n < NP; n++) {
			w[n] += w[n-1];
		}

		// save particle states to resample from
		for (int n = 0; n < NP; n++){
			copyParticle(&particles_old[n], &particles[n]);
		}

		// resampling
		for (int n = 0; n < NP; n++) {

			double w_r = randu() * w[NP-1];
			int i = 0;
			while (w_r > w[i]) {
				i++;
			}

			// i is now the index to copy state from
			copyParticle(&particles[n], &particles_old[i]);

		}
		
		State sMeans;
		getStateMeans(&sMeans, particles, NP);
		statemeans(t,0) = sMeans.S;
		statemeans(t,1) = sMeans.I;
		statemeans(t,2) = sMeans.R;

	}

	// Get particle results to pass back to R
	for (int n = 0; n < NP; n++) {
		statedata(n, 0) = particles[n].S;
		statedata(n, 1) = particles[n].I;
		statedata(n, 2) = particles[n].R;
		statedata(n, 3) = particles[n].B;
	}



	return Rcpp::List::create(	Rcpp::Named("statemeans") = statemeans,
                             	Rcpp::Named("statedata") = statedata 	);

}


/*	Use the Explicit Euler integration scheme to integrate SIR model forward in time
	double h 	- time step size
	double t0 	- start time
	double tn 	- stop time
	double * y 	- current system state; a three-component vector representing [S I R], susceptible-infected-recovered

	*/
void exp_euler_SIR(double h, double t0, double tn, int N, Particle * particle) {

	int num_steps = floor( (tn-t0) / h );

	double S = particle->S;
	double I = particle->I;
	double R = particle->R;

	double R0 	= particle->R0;
	double r 	= particle->r;
	double B0 	= R0 * r / N;
	double eta 	= particle->eta;
	double berr  = particle->berr;

	double B = particle->B;

	for(int i = 0; i < num_steps; i++) {

		B = exp( log(B) + eta*(log(B0) - log(B)) + berr*randn() );

		// get derivatives
		double dS = - B*S*I;
		double dI = B*S*I - r*I;
		double dR = r*I;

		// step forward by h
		S += h*dS;
		I += h*dI;
		R += h*dR;

	}

	particle->S = S;
	particle->I = I;
	particle->R = R;
	particle->B = B;

}


/*	Particle pertubation function to be run between iterations and passes

	*/
void perturbParticles(Particle * particles, int N, int NP, int passnum, double coolrate) {

	//double coolcoef = exp( - (double) passnum / coolrate );
	double coolcoef = pow(coolrate, passnum);
	

    double spreadR0 	= coolcoef * R0true / 10.0;
    double spreadr 		= coolcoef * rtrue / 10.0;
    double spreadsigma 	= coolcoef * merr / 10.0;
    double spreadIinit 	= coolcoef * I0 / 10.0;
    double spreadeta 	= coolcoef * etatrue / 10.0;
    double spreadberr 	= coolcoef * berrtrue / 10.0;


    double R0can, rcan, sigmacan, Iinitcan, etacan, berrcan;

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
    		Iinitcan = particles[n].Iinit + spreadIinit*randn();
    	} while (Iinitcan < 0 || Iinitcan > 500);
    	particles[n].Iinit = Iinitcan;
    	particles[n].Sinit = N - Iinitcan;

    }

}


/*	Convinience function for particle resampling process

	*/
void copyParticle(Particle * dst, Particle * src) {

	dst->B      = src->B;
	dst->S 		= src->S;
	dst->I 		= src->I;
	dst->R 		= src->R;

}

void particleDiagnostics(ParticleInfo * partInfo, Particle * particles, int NP) {

	double 	R0mean 		= 0.0,
    		rmean 		= 0.0,
    		sigmamean 	= 0.0,
    		etamean 	= 0.0,
    		berrmean 	= 0.0,
    		Sinitmean 	= 0.0,
    		Iinitmean 	= 0.0,
    		Rinitmean 	= 0.0;

    // means

    for (int n = 0; n < NP; n++) {

    	R0mean 		+= particles[n].R0;
    	rmean 		+= particles[n].r;
    	etamean 	+= particles[n].eta,
    	berrmean 	+= particles[n].berr,
    	sigmamean 	+= particles[n].sigma;
    	Sinitmean 	+= particles[n].Sinit;
    	Iinitmean 	+= particles[n].Iinit;
    	Rinitmean 	+= particles[n].Rinit;

    }

    R0mean 		/= NP;
    rmean 		/= NP;
    sigmamean 	/= NP;
    etamean 	/= NP;
    berrmean 	/= NP;
    Sinitmean 	/= NP;
    Iinitmean 	/= NP;
    Rinitmean 	/= NP;

    // standard deviations

    double 	R0sd 	= 0.0,
    		rsd 	= 0.0,
    		sigmasd = 0.0,
    		etasd	= 0.0,
    		berrsd	= 0.0,
    		Sinitsd = 0.0,
    		Iinitsd = 0.0,
    		Rinitsd = 0.0;

    for (int n = 0; n < NP; n++) {

    	R0sd 	+= ( particles[n].R0 - R0mean ) * ( particles[n].R0 - R0mean );
    	rsd 	+= ( particles[n].r - rmean ) * ( particles[n].r - rmean );
    	sigmasd += ( particles[n].sigma - sigmamean ) * ( particles[n].sigma - sigmamean );
    	etasd 	+= ( particles[n].eta - etamean ) * ( particles[n].eta - etamean );
    	berrsd 	+= ( particles[n].berr - berrmean ) * ( particles[n].berr - berrmean );
    	Sinitsd += ( particles[n].Sinit - Sinitmean ) * ( particles[n].Sinit - Sinitmean );
    	Iinitsd += ( particles[n].Iinit - Iinitmean ) * ( particles[n].Iinit - Iinitmean );
    	Rinitsd += ( particles[n].Rinit - Rinitmean ) * ( particles[n].Rinit - Rinitmean );

    }

    R0sd 		/= NP;
    rsd 		/= NP;
    sigmasd 	/= NP;
    etasd 		/= NP;
    berrsd 		/= NP;
    Sinitsd 	/= NP;
    Iinitsd 	/= NP;
    Rinitsd 	/= NP;

    partInfo->R0mean 	= R0mean;
    partInfo->R0sd 		= R0sd;
    partInfo->rmean 	= rmean;
    partInfo->rsd 		= rsd;
    partInfo->sigmamean = sigmamean;
    partInfo->sigmasd 	= sigmasd;
    partInfo->etamean 	= etamean;
    partInfo->etasd		= etasd;
    partInfo->berrmean	= berrmean;
    partInfo->berrsd	= berrsd;
    partInfo->Sinitmean = Sinitmean;
    partInfo->Sinitsd 	= Sinitsd;
    partInfo->Iinitmean = Iinitmean;
    partInfo->Iinitsd 	= Iinitsd;
    partInfo->Rinitmean = Rinitmean;
    partInfo->Rinitsd 	= Rinitsd;

}

double randu() {

	return (double) rand() / (double) RAND_MAX;

}

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
