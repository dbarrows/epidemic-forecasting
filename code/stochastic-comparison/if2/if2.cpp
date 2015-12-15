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
#define berrtrue	0.3			// real beta drift noise
#define merr 		10.0  		// expected measurement error
#define I0 			5.0			// Initial infected individuals

#include <Rcpp.h>
using namespace Rcpp;

struct State {
	float S;
	float I;
	float R;
};

struct Particle {
	float R0;
	float r;
	float sigma;
	float eta;
	float berr;
	float B;
	float S;
	float I;
	float R;
	float Sinit;
	float Iinit;
	float Rinit;
};

struct ParticleInfo {
	float R0mean;		float R0sd;
	float rmean;		float rsd;
	float sigmamean;	float sigmasd;
	float etamean;		float etasd;
	float berrmean;		float berrsd;
	float Sinitmean; 	float Sinitsd;
	float Iinitmean;	float Iinitsd;
	float Rinitmean;	float Rinitsd;
};


int timeval_subtract (double *result, struct timeval *x, struct timeval *y);
int check_float(float x,float y);
void exp_euler_SIR(float h, float t0, float tn, int N, Particle * particle);
void copyParticle(Particle * dst, Particle * src);
void perturbParticles(Particle * particles, int N, int NP, int passnum, float coolrate);
void particleDiagnostics(ParticleInfo * partInfo, Particle * particles, int NP);
NumericMatrix if2(NumericVector * data, int T, int N);
float randu();
float randn();

// [[Rcpp::export]]
Rcpp::List if2(NumericVector data, int T, int N, int NP, int nPasses, double coolrate) {

	//int 	NP 			= 10000;
	//int 	nPasses 	= 30;
	//float 	coolrate 	= 0.9;

	int 	i_infec 	= I0;

	NumericMatrix paramdata(NP, 8);
	NumericMatrix means(nPasses, 8);

	srand(time(NULL));		// Seed PRNG with system time

	float w[NP]; 			// particle weights

	Particle particles[NP]; 	// particle estimates for current step
	Particle particles_old[NP]; // intermediate particle states for resampling

	printf("Initializing particle states\n");

	// initialize particle parameter states (seeding)
	for (int n = 0; n < NP; n++) {

		float R0can, rcan, sigmacan, Iinitcan, etacan, berrcan;

		do {
			R0can = R0true + R0true*randn();
		} while (R0can < 0);
		particles[n].R0 = R0can;

		do {
			rcan = rtrue + rtrue*randn();
		} while (rcan < 0);
		particles[n].r = rcan;

		particles[n].B = (float) R0can * rcan / N;

		do {
			sigmacan = merr + merr*randn();
		} while (sigmacan < 0);
		particles[n].sigma = sigmacan;

		do {
			etacan = etatrue + 0.5*etatrue*randn();
		} while (etacan < 0 || etacan > 1);
		particles[n].eta = etacan;

		do {
			berrcan = berrtrue + 0.5*berrtrue*randn();
		} while (berrcan < 0);
		particles[n].berr = berrcan;

		do {
			Iinitcan = i_infec + i_infec*randn();
		} while (Iinitcan < 0 || N < Iinitcan);
		particles[n].Sinit = N - Iinitcan;
		particles[n].Iinit = Iinitcan;
		particles[n].Rinit = 0.0;

	}

	// START PASSES THROUGH DATA

	printf("Starting filter\n");
	printf("---------------\n");
	printf("Pass\n");


	for (int pass = 0; pass < nPasses; pass++) {

		printf("...%d / %d\n", pass, nPasses);

		// reset particle system evolution states
		for (int n = 0; n < NP; n++) {

			particles[n].S = particles[n].Sinit;
			particles[n].I = particles[n].Iinit;
			particles[n].R = particles[n].Rinit;
			particles[n].B = (float) particles[n].R0 * particles[n].r / N;

		}

		for (int t = 1; t <= (T+1); t++) {

			// generate individual predictions and weight
			for (int n = 0; n < NP; n++) {

				exp_euler_SIR(1.0/7.0, 0.0, 1.0, N, &particles[n]);

				float merr_par 	= particles[n].sigma;
				float y_diff 	= data[t] - particles[n].I;
				
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

				float w_r = randu() * w[NP-1];
				int i = 0;
				while (w_r > w[i]) {
					i++;
				}

				// i is now the index to copy state from
				copyParticle(&particles[n], &particles_old[i]);

			}
			
			// between-iteration perturbations, not after last time step
			if (t < (T+1))
			    perturbParticles(particles, N, NP, pass, coolrate);

		}

		ParticleInfo pInfo;
		particleDiagnostics(&pInfo, particles, NP);

		means(pass, 0) = pInfo.R0mean;
		means(pass, 1) = pInfo.rmean;
		means(pass, 2) = pInfo.sigmamean;
		means(pass, 3) = pInfo.etamean;
		means(pass, 4) = pInfo.berrmean;
		means(pass, 5) = pInfo.Sinitmean;
		means(pass, 6) = pInfo.Iinitmean;
		means(pass, 7) = pInfo.Rinitmean;


		// between-pass perturbations, not after last pass
		if (pass < (nPasses + 1))
		    perturbParticles(particles, N, NP, pass, coolrate);


	}

	ParticleInfo pInfo;
	particleDiagnostics(&pInfo, particles, NP);

	printf("Parameter results (mean | sd)\n");
	printf("-----------------------------\n");
	printf("R0        %f %f\n", pInfo.R0mean, pInfo.R0sd);
	printf("r         %f %f\n", pInfo.rmean, pInfo.rsd);
	printf("sigma     %f %f\n", pInfo.sigmamean, pInfo.sigmasd);
	printf("eta       %f %f\n", pInfo.etamean, pInfo.etasd);
	printf("berr    %f %f\n", pInfo.berrmean, pInfo.berrsd);
	printf("S_init  %f %f\n", pInfo.Sinitmean, pInfo.Sinitsd);
	printf("I_init    %f %f\n", pInfo.Iinitmean, pInfo.Iinitsd);
	printf("R_init    %f %f\n", pInfo.Rinitmean, pInfo.Rinitsd);

	printf("\n");



	// Get particle results to pass back to R

	for (int n = 0; n < NP; n++) {

		paramdata(n, 0) = particles[n].R0;
		paramdata(n, 1) = particles[n].r;
		paramdata(n, 2) = particles[n].sigma;
		paramdata(n, 3) = particles[n].eta;
		paramdata(n, 4) = particles[n].berr;
		paramdata(n, 5) = particles[n].Sinit;
		paramdata(n, 6) = particles[n].Iinit;
		paramdata(n, 7) = particles[n].Rinit;

	}

	return Rcpp::List::create(	Rcpp::Named("paramdata") = paramdata, 
                             	Rcpp::Named("means") = means);

}


/*	Use the Explicit Euler integration scheme to integrate SIR model forward in time
	float h 	- time step size
	float t0 	- start time
	float tn 	- stop time
	float * y 	- current system state; a three-component vector representing [S I R], susceptible-infected-recovered

	*/
void exp_euler_SIR(float h, float t0, float tn, int N, Particle * particle) {

	int num_steps = floor( (tn-t0) / h );

	float S = particle->S;
	float I = particle->I;
	float R = particle->R;

	float R0 	= particle->R0;
	float r 	= particle->r;
	float B0 	= R0 * r / N;
	float eta 	= particle->eta;
	float berr  = particle->berr;

	float B = particle->B;

	for(int i = 0; i < num_steps; i++) {

		B = exp( log(B) + eta*(log(B0) - log(B)) + berr*randn() );

		// get derivatives
		float dS = - B*S*I;
		float dI = B*S*I - r*I;
		float dR = r*I;

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
void perturbParticles(Particle * particles, int N, int NP, int passnum, float coolrate) {

	//float coolcoef = exp( - (float) passnum / coolrate );
	float coolcoef = pow(coolrate, passnum);
	

    float spreadR0 		= coolcoef * R0true / 10.0;
    float spreadr 		= coolcoef * rtrue / 10.0;
    float spreadsigma 	= coolcoef * merr / 10.0;
    float spreadIinit 	= coolcoef * I0 / 10.0;
    float spreadeta 	= coolcoef * etatrue / 10.0;
    float spreadberr 	= coolcoef * berrtrue / 10.0;


    float R0can, rcan, sigmacan, Iinitcan, etacan, berrcan;

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
    		etacan = particles[n].eta + 0.5*spreadeta*randn();
    	} while (etacan < 0 || etacan > 1);
    	particles[n].eta = etacan;

    	do {
    		berrcan = particles[n].berr + 0.5*spreadberr*randn();
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

	dst->R0 	= src->R0;
	dst->r 		= src->r;
	dst->sigma 	= src->sigma;
	dst->eta 	= src->eta;
	dst->berr 	= src->berr;
	dst->B      = src->B;
	dst->S 		= src->S;
	dst->I 		= src->I;
	dst->R 		= src->R;
	dst->Sinit  = src->Sinit;
	dst->Iinit  = src->Iinit;
	dst->Rinit  = src->Rinit;

}

void particleDiagnostics(ParticleInfo * partInfo, Particle * particles, int NP) {

	float 	R0mean 		= 0.0,
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

    float 	R0sd 	= 0.0,
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

float randu() {

	return (float) rand() / (float) RAND_MAX;

}


/*	Return a normally distributed random number with mean 0 and standard deviation 1
	Uses the polar form of the Box-Muller transformation
	From http://www.design.caltech.edu/erik/Misc/Gaussian.html
	*/
float randn() {

	float x1, x2, w, y1;

	do {
		x1 = 2.0 * randu() - 1.0;
		x2 = 2.0 * randu() - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;

	return y1;

}
