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

#define Treal 	100			// time to simulate over
#define R0true 	3.0			// infectiousness
#define rtrue 	0.1			// recovery rate
#define Nreal 	500.0 		// population size
#define merr 	10.0  		// expected measurement error
#define I0 		5.0			// Initial infected individuals

#include <Rcpp.h>
using namespace Rcpp;


struct Particle {
	float R0;
	float r;
	float sigma;
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
	float Sinitmean; 	float Sinitsd;
	float Iinitmean;	float Iinitsd;
	float Rinitmean;	float Rinitsd;
};


int timeval_subtract (double *result, struct timeval *x, struct timeval *y);
int check_float(float x,float y);
void exp_euler_SIR(float h, float t0, float tn, int N, Particle * particle);
void copyParticle(Particle * dst, Particle * src);
void perturbParticles(Particle * particles, int N, int NP, int passnum, float coolrate);
bool isCollapsed(Particle * particles, int NP);
void particleDiagnostics(ParticleInfo * partInfo, Particle * particles, int NP);
NumericMatrix if2(NumericVector * data, int T, int N);
float randu();
float randn();

int main(int argc, char *argv[]) {

	int T, N;

	std::ifstream read("sysinfo.dat");
	read >> T;
	read >> N;

	printf("T: %d, N: %d\n", T, N);

	/*

	int 	T 			= Treal;
	int 	N 			= Nreal;
	float 	i_infec 	= I0;
	int 	Tlim 		= T;

	Particle particle_true;
	particle_true.R0 	= R0true;
	particle_true.r  	= rtrue;
	particle_true.sigma = merr;
	particle_true.S  	= N - i_infec;
	particle_true.I  	= i_infec;
	particle_true.R  	= 0;

	srand(time(NULL));								// Seed PRNG with system time

	float y_true[T];	// true number of infected peeps
	float y_noise[T];	// true number of infected peeps with observation noise

	printf("System parameters\n");
	printf("-----------------\n");
	printf("R0:    %f\n", R0true);
	printf("r:     %f\n", rtrue);
	printf("merr:  %f\n", merr);

	// generate our true trajectory and noisy observation data
	y_true[0] = particle_true.I;
	y_noise[0] = y_true[0] +  merr*randn();
	if (y_noise[0] < 0)
		y_noise[0] = 0;
	for (int i = 1; i < T; i++) {
		exp_euler_SIR( 1.0/100, 0.0, 1.0, N, &particle_true);
		y_true[i] = particle_true.I;
		y_noise[i] = y_true[i] + merr*randn();
		if (y_noise[i] < 0)
			y_noise[i] = 0;
	}

	*/

	//if2(y_noise, Tlim, N);


}

// [[Rcpp::export]]
NumericMatrix if2(NumericVector data, int T, int N) {

	int 	NP 			= 5000;
	int 	nPasses 	= 50;
	float 	coolrate 	= 8;

	int 	i_infec 	= I0;

	NumericMatrix paramdata(NP, 6);

	srand(time(NULL));		// Seed PRNG with system time

	float w[NP]; 			// particle weights

	Particle particles[NP]; 	// particle estimates for current step
	Particle particles_old[NP]; // intermediate particle states for resampling

	printf("Initializing particle states\n");

	// initialize particle parameter states (seeding)
	for (int n = 0; n < NP; n++) {

		float R0can, rcan, sigmacan, Iinitcan;

		do {
			R0can = R0true + R0true*randn();
		} while (R0can < 0);
		particles[n].R0 = R0can;

		do {
			rcan = rtrue + rtrue*randn();
		} while (rcan < 0);
		particles[n].r = rcan;

		do {
			sigmacan = merr + merr*randn();
		} while (sigmacan < 0);
		particles[n].sigma = sigmacan;

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

		perturbParticles(particles, N, NP, pass, coolrate);

		// initialize particle system states
		for (int n = 0; n < NP; n++) {

			particles[n].S = particles[n].Sinit;
			particles[n].I = particles[n].Iinit;
			particles[n].R = particles[n].Rinit;

		}

		// between-pass perturbations

		for (int t = 1; t < T; t++) {

			// between-iteration perturbations
			perturbParticles(particles, N, NP, pass, coolrate);

			// generate individual predictions and weight
			for (int n = 0; n < NP; n++) {

				exp_euler_SIR(1.0/10.0, 0.0, 1.0, N, &particles[n]);

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

		}

	}

	ParticleInfo pInfo;
	particleDiagnostics(&pInfo, particles, NP);

	printf("Parameter results (mean | sd)\n");
	printf("-----------------------------\n");
	printf("R0        %f %f\n", pInfo.R0mean, pInfo.R0sd);
	printf("r         %f %f\n", pInfo.rmean, pInfo.rsd);
	printf("sigma     %f %f\n", pInfo.sigmamean, pInfo.sigmasd);
	printf("S_init  %f %f\n", pInfo.Sinitmean, pInfo.Sinitsd);
	printf("I_init    %f %f\n", pInfo.Iinitmean, pInfo.Iinitsd);
	printf("R_init    %f %f\n", pInfo.Rinitmean, pInfo.Rinitsd);

	printf("\n");



	// Get particle results to pass back to R

	for (int n = 0; n < NP; n++) {

		paramdata(n, 0) = particles[n].R0;
		paramdata(n, 1) = particles[n].r;
		paramdata(n, 2) = particles[n].sigma;
		paramdata(n, 3) = particles[n].Sinit;
		paramdata(n, 4) = particles[n].Iinit;
		paramdata(n, 5) = particles[n].Rinit;

	}

	return paramdata;

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
	float B 	= R0 * r / N;

	for(int i = 0; i < num_steps; i++) {
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

}


/*	Particle pertubation function to be run between iterations and passes

	*/
void perturbParticles(Particle * particles, int N, int NP, int passnum, float coolrate) {

	float coolcoef = exp( - (float) passnum / coolrate );

    float spreadR0 		= coolcoef * R0true / 3.0;
    float spreadr 		= coolcoef * rtrue 	/ 3.0;
    float spreadsigma 	= coolcoef * merr 	/ 3.0;
    float spreadIinit 	= coolcoef * 5.0 	/ 3.0;

    float R0can, rcan, sigmacan, Iinitcan;

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
	dst->S 		= src->S;
	dst->I 		= src->I;
	dst->R 		= src->R;
	dst->Sinit  = src->Sinit;
	dst->Iinit  = src->Iinit;
	dst->Rinit  = src->Rinit;

}


/*	Checks to see if particles are collapsed
	This is done by checking if the standard deviations between the particles' parameter
	values are significantly close to one another. Spread threshold may need to be tuned.

	*/
bool isCollapsed(Particle * particles, int NP) {

	bool retVal;

	float R0mean = 0, rmean = 0, sigmamean = 0, Sinitmean = 0, Iinitmean = 0, Rinitmean = 0;

    // means

    for (int n = 0; n < NP; n++) {

    	R0mean 		+= particles[n].R0;
    	rmean 		+= particles[n].r;
    	sigmamean 	+= particles[n].sigma;
    	Sinitmean 	+= particles[n].Sinit;
    	Iinitmean 	+= particles[n].Iinit;
    	Rinitmean 	+= particles[n].Rinit;

    }

    R0mean 		/= NP;
    rmean 		/= NP;
    sigmamean 	/= NP;
    Sinitmean 	/= NP;
    Iinitmean 	/= NP;
    Rinitmean 	/= NP;

	float 	R0sd = 0, rsd = 0, sigmasd = 0, Sinitsd = 0, Iinitsd = 0, Rinitsd = 0;

    for (int n = 0; n < NP; n++) {

    	R0sd 	+= ( particles[n].R0 - R0mean ) * ( particles[n].R0 - R0mean );
    	rsd 	+= ( particles[n].r - rmean ) * ( particles[n].r - rmean );
    	sigmasd += ( particles[n].sigma - sigmamean ) * ( particles[n].sigma - sigmamean );
    	Sinitsd += ( particles[n].Sinit - Sinitmean ) * ( particles[n].Sinit - Sinitmean );
    	Iinitsd += ( particles[n].Iinit - Iinitmean ) * ( particles[n].Iinit - Iinitmean );
    	Rinitsd += ( particles[n].Rinit - Rinitmean ) * ( particles[n].Rinit - Rinitmean );

    }

    R0sd 		/= NP;
    rsd 		/= NP;
    sigmasd 	/= NP;
    Sinitsd 	/= NP;
    Iinitsd 	/= NP;
    Rinitsd 	/= NP;

    if ( (R0sd + rsd + sigmasd) < 1e-5)
    	retVal = true;
    else
    	retVal = false;

    return retVal;

}

void particleDiagnostics(ParticleInfo * partInfo, Particle * particles, int NP) {

	float 	R0mean 		= 0.0,
    		rmean 		= 0.0,
    		sigmamean 	= 0.0,
    		Sinitmean 	= 0.0,
    		Iinitmean 	= 0.0,
    		Rinitmean 	= 0.0;

    // means

    for (int n = 0; n < NP; n++) {

    	R0mean 		+= particles[n].R0;
    	rmean 		+= particles[n].r;
    	sigmamean 	+= particles[n].sigma;
    	Sinitmean 	+= particles[n].Sinit;
    	Iinitmean 	+= particles[n].Iinit;
    	Rinitmean 	+= particles[n].Rinit;

    }

    R0mean 		/= NP;
    rmean 		/= NP;
    sigmamean 	/= NP;
    Sinitmean 	/= NP;
    Iinitmean 	/= NP;
    Rinitmean 	/= NP;

    // standard deviations

    float 	R0sd 	= 0.0,
    		rsd 	= 0.0,
    		sigmasd = 0.0,
    		Sinitsd = 0.0,
    		Iinitsd = 0.0,
    		Rinitsd = 0.0;

    for (int n = 0; n < NP; n++) {

    	R0sd 	+= ( particles[n].R0 - R0mean ) * ( particles[n].R0 - R0mean );
    	rsd 	+= ( particles[n].r - rmean ) * ( particles[n].r - rmean );
    	sigmasd += ( particles[n].sigma - sigmamean ) * ( particles[n].sigma - sigmamean );
    	Sinitsd += ( particles[n].Sinit - Sinitmean ) * ( particles[n].Sinit - Sinitmean );
    	Iinitsd += ( particles[n].Iinit - Iinitmean ) * ( particles[n].Iinit - Iinitmean );
    	Rinitsd += ( particles[n].Rinit - Rinitmean ) * ( particles[n].Rinit - Rinitmean );

    }

    R0sd 		/= NP;
    rsd 		/= NP;
    sigmasd 	/= NP;
    Sinitsd 	/= NP;
    Iinitsd 	/= NP;
    Rinitsd 	/= NP;

    partInfo->R0mean 	= R0mean;
    partInfo->R0sd 		= R0sd;
    partInfo->sigmamean = sigmamean;
    partInfo->sigmasd 	= sigmasd;
    partInfo->rmean 	= rmean;
    partInfo->rsd 		= rsd;
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
