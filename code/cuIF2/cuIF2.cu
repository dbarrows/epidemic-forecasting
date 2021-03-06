/*	Dexter Barrows
	dbarrows.github.io
	McMaster University
	2016

	*/

#include <cuda.h>
#include <iostream>
#include <fstream>
#include <curand.h>
#include <curand_kernel.h>
#include <string>
#include <sstream>
#include <cmath>

#include "timer.h"
#include "rand.h"
#include "readdata.h"

#define NP 			(2*2500) 	// number of particles
#define N 			500.0		// population size
#define R0true 		3.0			// infectiousness
#define rtrue 		0.1			// recovery rate
#define etatrue 	0.5			// real drift attraction strength
#define berrtrue	0.5			// real beta drift noise
#define phitrue 	0.5 		// real connectivity strength
#define merr 		10.0  		// expected measurement error
#define I0 			5.0			// Initial infected individuals
#define PSC 		0.5 		// sensitive parameter perturbation scaling
#define NLOC 		10

#define PI 		3.141592654f

// Wrapper for CUDA calls, from CUDA API
// Modified to also print the error code and string
# define CUDA_CALL(x) do { if ((x) != cudaSuccess ) {					\
	std::cout << " Error at " << __FILE__ << ":" << __LINE__ << std::endl;		\
	std::cout << " Error was " << x << " " << cudaGetErrorString(x) << std::endl;	\
	return EXIT_FAILURE ;}} while (0)									\

typedef struct {
	float R0;
	float r;
	float sigma;
	float eta;
	float berr;
	float phi;
	float S[NLOC];
	float I[NLOC];
	float R[NLOC];
	float B[NLOC];
	float Iinit[NLOC];
	curandState randState; 	// PRNG state
} Particle;

__host__ std::string getHRmemsize (size_t memsize);
__host__ std::string getHRtime (float runtime);

__device__ void exp_euler_SSIR(float h, float t0, float tn, Particle * particle, int * neinum, int * neibmat, int nloc);
__device__ void copyParticle(Particle * dst, Particle * src, int nloc);


/* 	Initialize all PRNG states, get starting state vector using initial distribution
	*/
__global__ void initializeParticles (Particle * particles, int nloc) {

	int id 	= blockIdx.x*blockDim.x + threadIdx.x;	// global thread ID

	if (id < NP) {

		// initialize PRNG state
		curandState state;
		curand_init(id, 0, 0, &state);

		float R0can, rcan, sigmacan, Iinitcan, etacan, berrcan, phican;

		do {
			R0can = R0true + R0true*curand_normal(&state);
		} while (R0can < 0);
		particles[id].R0 = R0can;

		do {
			rcan = rtrue + rtrue*curand_normal(&state);
		} while (rcan < 0);
		particles[id].r = rcan;

		for (int loc = 0; loc < nloc; loc++)
			particles[id].B[loc] = (float) R0can * rcan / N;

		do {
			sigmacan = merr + merr*curand_normal(&state);
		} while (sigmacan < 0);
		particles[id].sigma = sigmacan;

		do {
			etacan = etatrue + PSC*etatrue*curand_normal(&state);
		} while (etacan < 0 || etacan > 1);
		particles[id].eta = etacan;

		do {
			berrcan = berrtrue + PSC*berrtrue*curand_normal(&state);
		} while (berrcan < 0);
		particles[id].berr = berrcan;

		do {
			phican = phitrue + PSC*phitrue*curand_normal(&state);
		} while (phican <= 0 || phican >= 1);
		particles[id].phi = phican;

		for (int loc = 0; loc < nloc; loc++) {
			do {
				Iinitcan = I0 + I0*curand_normal(&state);
			} while (Iinitcan < 0 || N < Iinitcan);
			particles[id].Iinit[loc] = Iinitcan;
		}

		particles[id].randState = state;

	}

}

__global__ void resetStates (Particle * particles, int nloc) {

	int id 	= blockIdx.x*blockDim.x + threadIdx.x;	// global thread ID

	if (id < NP) {

		for (int loc = 0; loc < nloc; loc++) {
			particles[id].S[loc] = N - particles[id].Iinit[loc];
			particles[id].I[loc] = particles[id].Iinit[loc];
			particles[id].R[loc] = 0.0;
		}

	}

}

__global__ void clobberParams (Particle * particles, int nloc) {

	int id 	= blockIdx.x*blockDim.x + threadIdx.x;	// global thread ID

	if (id < NP) {

		particles[id].R0 = R0true;
		particles[id].r = rtrue;
		particles[id].sigma = merr;
		particles[id].eta = etatrue;
		particles[id].berr = berrtrue;
		particles[id].phi = phitrue;

		for (int loc = 0; loc < nloc; loc++) {
			particles[id].Iinit[loc] = I0;
		}

	}

}


/* 	Project particles forward, perturb, and save weight based on data
	int t - time step number (1,...,T)
	*/
__global__ void project (Particle * particles, int * neinum, int * neibmat, int nloc) {

	int id = blockIdx.x*blockDim.x + threadIdx.x;	// global id

	if (id < NP) {
		// project forward
		exp_euler_SSIR(1.0/7.0, 0.0, 1.0, &particles[id], neinum, neibmat, nloc);
	}

}

__global__ void weight(float * data, Particle * particles, double * w, int t, int T, int nloc) {

	int id = blockIdx.x*blockDim.x + threadIdx.x;	// global id

	if (id < NP) {

		float merr_par = particles[id].sigma;

		// Get weight and save
		double w_local = 1.0;
		for (int loc = 0; loc < nloc; loc++) {
			float y_diff = data[loc*T + t] - particles[id].I[loc];
			w_local *= 1.0/(merr_par*sqrt(2.0*PI)) * exp( - y_diff*y_diff / (2.0*merr_par*merr_par) );
		}

		w[id] = w_local;

	}

}

__global__ void stashParticles (Particle * particles, Particle * particles_old, int nloc) {

	int id = blockIdx.x*blockDim.x + threadIdx.x;	// global id
	
	if (id < NP) {
		// COPY PARTICLE
		copyParticle(&particles_old[id], &particles[id], nloc);
	}

}


/* 	The 0th thread will perform cumulative sum on the weights.
	There may be a faster way to do this, will investigate.
	*/
__global__ void cumsumWeights (double * w) {

	int id 	= blockIdx.x*blockDim.x + threadIdx.x;	// global thread ID

	// compute cumulative weights
	if (id == 0) {
		for (int i = 1; i < NP; i++)
			w[i] += w[i-1];
	}

}


/* 	Resample from all particle states within cell
	*/
__global__ void resample (Particle * particles, Particle * particles_old, double * w, int nloc) {

	int id 	= blockIdx.x*blockDim.x + threadIdx.x;

	if (id < NP) {

		// resampling proportional to weights
		double w_r = curand_uniform(&particles[id].randState) * w[NP-1];
		int i = 0;
		while (w_r > w[i]) {
			i++;
		}	

		// i is now the index of the particle to copy from
		copyParticle(&particles[id], &particles_old[i], nloc);

	}

}

// launch this with probably just nloc threads... block structure/size probably not important
__global__ void reduceStates (Particle * particles, float * countmeans, int t, int T, int nloc) {

	int id 	= blockIdx.x*blockDim.x + threadIdx.x;

	if (id < nloc) {

		int loc = id;

		double countmean_local = 0.0;
		for (int n = 0; n < NP; n++) {
			countmean_local += particles[n].I[loc] / NP;
		}

		countmeans[loc*T + t] = (float) countmean_local;

	}

}

__global__ void perturbParticles(Particle * particles, int nloc, int passnum, double coolrate) {

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

    int id 	= blockIdx.x*blockDim.x + threadIdx.x;

    if (id < NP) {

		do {
			R0can = particles[id].R0 + spreadR0*curand_normal(&particles[id].randState);
		} while (R0can < 0);
		particles[id].R0 = R0can;

		do {
			rcan = particles[id].r + spreadr*curand_normal(&particles[id].randState);
		} while (rcan < 0);
		particles[id].r = rcan;

		do {
			sigmacan = particles[id].sigma + spreadsigma*curand_normal(&particles[id].randState);
		} while (sigmacan < 0);
		particles[id].sigma = sigmacan;

		do {
			etacan = particles[id].eta + PSC*spreadeta*curand_normal(&particles[id].randState);
		} while (etacan < 0 || etacan > 1);
		particles[id].eta = etacan;

		do {
			berrcan = particles[id].berr + PSC*spreadberr*curand_normal(&particles[id].randState);
		} while (berrcan < 0);
		particles[id].berr = berrcan;

		do {
			phican = particles[id].phi + PSC*spreadphi*curand_normal(&particles[id].randState);
		} while (phican <= 0 || phican >= 1);
		particles[id].phi = phican;

		for (int loc = 0; loc < nloc; loc++) {
	    	do {
	    		Iinitcan = particles[id].Iinit[loc] + spreadIinit*curand_normal(&particles[id].randState);
	    	} while (Iinitcan < 0 || Iinitcan > 500);
	    	particles[id].Iinit[loc] = Iinitcan;
	    }

	}

}


int main (int argc, char *argv[]) {


	int T, nloc;

	double restime;
	struct timeval tdr0, tdr1, tdrMaster;

	// Parse arguments **********************************************

	if (argc < 4) {
		std::cout << "Not enough arguments" << std::endl;
		return 0;
	}

	std::string arg1(argv[1]); 	// infection counts
	std::string arg2(argv[2]);	// neighbour counts
	std::string arg3(argv[3]);	// neighbour indices
	std::string arg4(argv[4]); 	// outfile: params + runtime

	std::cout << "Arguments:" << std::endl;
	std::cout << "Infection data: 	 " << arg1 << std::endl;
	std::cout << "Neighbour counts:  " << arg2 << std::endl;
	std::cout << "Neighbour indices: " << arg3 << std::endl;
	std::cout << "Outfile            " << arg4 << std::endl;

	// **************************************************************


	// Read count data **********************************************

	std::cout << "Getting count data" << std::endl;
	float * data = getDataFloat(arg1, &T, &nloc);
	size_t datasize = nloc*T*sizeof(float);

	// **************************************************************

	// Read neinum matrix data **************************************

	std::cout << "Getting neighbour count data" << std::endl;
	int * neinum = getDataInt(arg2, NULL, NULL);
	size_t neinumsize = nloc * sizeof(int);

	// **************************************************************

	// Read neibmat matrix data *************************************

	std::cout << "Getting neighbour count data" << std::endl;
	int * neibmat = getDataInt(arg3, NULL, NULL);
	size_t neibmatsize = nloc * nloc * sizeof(int);

	// **************************************************************

	// *****************************************************************************************************

    // start timing
	gettimeofday (&tdr0, NULL);

	// CUDA data ****************************************************

	std::cout << "Allocating device storage" << std::endl;

	float 		* d_data;			// device copy of data
	Particle 	* particles;		// particles
	Particle 	* particles_old; 	// intermediate particle states
	double 		* w;				// weights
	int         * d_neinum; 		// device copy of adjacency matrix
	int 		* d_neibmat; 		// device copy of neighbour counts matrix
	float 		* countmeans; 		// host copy of reduced infection count means from last pass
	float 		* d_countmeans; 	// device copy of reduced infection count means from last pass

	CUDA_CALL( cudaMalloc( (void**) &d_data 		, datasize )			);
	CUDA_CALL( cudaMalloc( (void**) &particles 		, NP*sizeof(Particle)) 	);
	CUDA_CALL( cudaMalloc( (void**) &particles_old 	, NP*sizeof(Particle)) 	);
	CUDA_CALL( cudaMalloc( (void**) &w 				, NP*sizeof(double)) 	);
	CUDA_CALL( cudaMalloc( (void**) &d_neinum 		, neinumsize) 			);
	CUDA_CALL( cudaMalloc( (void**) &d_neibmat 		, neibmatsize) 			);
	CUDA_CALL( cudaMalloc( (void**) &d_countmeans 	, nloc*T*sizeof(float)) );


	gettimeofday (&tdr1, NULL);
    timeval_subtract (&restime, &tdr1, &tdr0);

    std::cout << "\t" << getHRtime(restime) << std::endl;

	size_t avail, total;
	cudaMemGetInfo( &avail, &total );
	size_t used = total - avail;

	std::cout << "\t[" << getHRmemsize(used) << "] used of [" << getHRmemsize(total) << "]" <<std::endl;

	std::cout << "Copying data to device" << std::endl;

	gettimeofday (&tdr0, NULL);

	CUDA_CALL( cudaMemcpy(d_data	, data 		, datasize		, cudaMemcpyHostToDevice)	);
	CUDA_CALL( cudaMemcpy(d_neinum	, neinum 	, neinumsize	, cudaMemcpyHostToDevice)	);
	CUDA_CALL( cudaMemcpy(d_neibmat , neibmat 	, neibmatsize	, cudaMemcpyHostToDevice)	);

	gettimeofday (&tdr1, NULL);
    timeval_subtract (&restime, &tdr1, &tdr0);

    std::cout << "\t" << getHRtime(restime) << std::endl;

	// **************************************************************



	// Initialize particles *****************************************

	std::cout << "Initializing particles" << std::endl;

	//gettimeofday (&tdr0, NULL);

	int nThreads 	= 32;
	int nBlocks 	= ceil( (float) NP / nThreads);

	initializeParticles <<< nBlocks, nThreads >>> (particles, nloc);
	CUDA_CALL( cudaGetLastError() );
	CUDA_CALL( cudaDeviceSynchronize() );

	initializeParticles <<< nBlocks, nThreads >>> (particles_old, nloc);
	CUDA_CALL( cudaGetLastError() );
	CUDA_CALL( cudaDeviceSynchronize() );

	//gettimeofday (&tdr1, NULL);
    //timeval_subtract (&restime, &tdr1, &tdr0);
    //std::cout << "\t" << getHRtime(restime) << std::endl;

    cudaMemGetInfo( &avail, &total );
	used = total - avail;
	std::cout << "\t[" << getHRmemsize(used) << "] used of [" << getHRmemsize(total) << "]" <<std::endl;

	// **************************************************************

	// Starting filtering *******************************************

	for (int pass = 0; pass < 50; pass++) {

		nThreads 	= 32;
		nBlocks 	= ceil( (float) NP / nThreads);

		resetStates <<< nBlocks, nThreads >>> (particles, nloc);
		CUDA_CALL( cudaGetLastError() );
		CUDA_CALL( cudaDeviceSynchronize() );

		nThreads = 1;
		nBlocks  = 10;

		if (pass == 49) {
			reduceStates <<< nBlocks, nThreads >>> (particles, d_countmeans, 0, T, nloc);
			CUDA_CALL( cudaGetLastError() );
			CUDA_CALL( cudaDeviceSynchronize() );
		}

		int Tlim = T;

		for (int t = 1; t < Tlim; t++) {

			// Projection ************************************************

			nThreads 	= 32;
			nBlocks 	= ceil( (float) NP / nThreads);

			project <<< nBlocks, nThreads >>> (particles, d_neinum, d_neibmat, nloc);
			CUDA_CALL( cudaGetLastError() );
			CUDA_CALL( cudaDeviceSynchronize() );

		    // Weighting *************************************************

			nThreads 	= 32;
			nBlocks 	= ceil( (float) NP / nThreads);

			weight <<< nBlocks, nThreads >>>(d_data, particles, w, t, T, nloc);
			CUDA_CALL( cudaGetLastError() );
			CUDA_CALL( cudaDeviceSynchronize() );

		    // Cumulative sum ********************************************

			nThreads 	= 1;
			nBlocks 	= 1;

			cumsumWeights <<< nBlocks, nThreads >>> (w);
			CUDA_CALL( cudaGetLastError() );
			CUDA_CALL( cudaDeviceSynchronize() );

		    // Save particles for resampling from *************************

		    nThreads 	= 32;
			nBlocks 	= ceil( (float) NP / nThreads);

			stashParticles <<< nBlocks, nThreads >>> (particles, particles_old, nloc); 
			CUDA_CALL( cudaGetLastError() );
			CUDA_CALL( cudaDeviceSynchronize() );


		    // Resampling *************************************************

			nThreads 	= 32;
			nBlocks 	= ceil( (float) NP/ nThreads);

			resample <<< nBlocks, nThreads >>> (particles, particles_old, w, nloc);
			CUDA_CALL( cudaGetLastError() );
			CUDA_CALL( cudaDeviceSynchronize() );

		    // Reduction **************************************************

		    if (pass == 49) {

		    	nThreads = 1;
		    	nBlocks  = 10;

		    	reduceStates <<< nBlocks, nThreads >>> (particles, d_countmeans, t, T, nloc);
		    	CUDA_CALL( cudaGetLastError() );
				CUDA_CALL( cudaDeviceSynchronize() );

			}

		    // Perturb particles ******************************************

		    nThreads 	= 32;
			nBlocks 	= ceil( (float) NP/ nThreads);

		    perturbParticles <<< nBlocks, nThreads >>> (particles, nloc, pass, 0.975);
		    CUDA_CALL( cudaGetLastError() );
			CUDA_CALL( cudaDeviceSynchronize() );


		} // end time

	} // end pass

	std::cout.precision(10);

	countmeans = (float*) malloc (nloc*T*sizeof(float));
	cudaMemcpy(countmeans, d_countmeans, nloc*T*sizeof(float), cudaMemcpyDeviceToHost);

	// stop master timer and print

	gettimeofday (&tdrMaster, NULL);
	timeval_subtract(&restime, &tdrMaster, &tdr0);
	std::cout << "Time: " << getHRtime(restime) << std::endl;
	std::cout << "Rawtime: " << restime << std::endl;

	// Write results out

	std::string filename = arg4;

	std::cout << "Writing results to file '" << filename << "' ..." << std::endl;

	std::ofstream outfile;
	outfile.open(filename.c_str());

	for(int loc = 0; loc < nloc; loc++) {
		for (int t = 0; t < T; t++) {
			outfile << countmeans[loc*T + t] << " ";
		}
		outfile << std::endl;
	}

	outfile.close();

	cudaFree(d_data);
	cudaFree(particles);
	cudaFree(particles_old);
	cudaFree(w);
	cudaFree(d_neinum);
	cudaFree(d_neibmat);
	cudaFree(d_countmeans);

	exit (EXIT_SUCCESS);

}


/*	Use the Explicit Euler integration scheme to integrate SIR model forward in time
	float h 	- time step size
	float t0 	- start time
	float tn 	- stop time
	float * y 	- current system state; a three-component vector representing [S I R], susceptible-infected-recovered
	*/
__device__ void exp_euler_SSIR(float h, float t0, float tn, Particle * particle, int * neinum, int * neibmat, int nloc) {

	int num_steps = floor( (tn-t0) / h );

	float * S = particle->S;
	float * I = particle->I;
	float * R = particle->R;
	float * B = particle->B;

	// create last state vectors
	float * S_last = (float*) malloc (nloc*sizeof(float));
	float * I_last = (float*) malloc (nloc*sizeof(float));
	float * R_last = (float*) malloc (nloc*sizeof(float));
	float * B_last = (float*) malloc (nloc*sizeof(float));

	float R0 	= particle->R0;
	float r 	= particle->r;
	float B0 	= R0 * r / N;
	float eta 	= particle->eta;
	float berr  = particle->berr;
	float phi   = particle->phi;

	for(int t = 0; t < num_steps; t++) {

		for (int loc = 0; loc < nloc; loc++) {
			S_last[loc] = S[loc];
			I_last[loc] = I[loc];
			R_last[loc] = R[loc];
			B_last[loc] = B[loc];
		}

		for (int loc = 0; loc < nloc; loc++) {

			B[loc] = exp( log(B_last[loc]) + eta*(log(B0) - log(B_last[loc])) + berr*curand_normal(&(particle->randState)) );

			int n = neinum[loc];
        	float sphi = 1.0 - phi*( (float) n/(n+1.0) );
        	float ophi = phi/(n+1.0);

        	float nBIsum = 0.0;
        	for (int j = 0; j < n; j++)
        		nBIsum += B_last[neibmat[nloc*loc + j]-1] * I_last[neibmat[nloc*loc + j]-1];

        	float BSI = S_last[loc]*( sphi*B_last[loc]*I_last[loc] + ophi*nBIsum );
        	float rI  = r*I_last[loc];

			// get derivatives
			float dS = - BSI;
			float dI = BSI - rI;
			float dR = rI;

			// step forward by h
			S[loc] += h*dS;
			I[loc] += h*dI;
			R[loc] += h*dR;

		}

	}

	free(S_last);
	free(I_last);
	free(R_last);
	free(B_last);

}

/*	Convinience function for particle resampling process
	*/
__device__ void copyParticle(Particle * dst, Particle * src, int nloc) {

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

/*	Convert memory size in bytes to human-readable format
	*/
std::string getHRmemsize (size_t memsize) {

	std::stringstream ss;
	std::string valstring;

	int kb = 1024;
	int mb = kb*1024;
	int gb = mb*1024;
	
	if (memsize <= kb)
		ss << memsize << " B";
	else if (memsize > kb && memsize <= mb)
		ss << (float) memsize/ kb << " KB";
	else if (memsize > mb && memsize <= gb)
		ss << (float) memsize/ mb << " MB";
	else
		ss << (float) memsize/ gb << " GB";

	valstring = ss.str();
	
	return valstring;

}


/*	Convert time in seconds to human readable format
	*/
std::string getHRtime (float runtime) {

	std::stringstream ss;
	std::string valstring;

	int mt = 60;
	int ht = mt*60;
	int dt = ht*24;
	
	if (runtime <= mt)
		ss << runtime << " s";
	else if (runtime > mt && runtime <= ht)
		ss << runtime/mt << " m";
	else if (runtime > ht && runtime <= dt)
		ss << runtime/dt << " h";
	else
		ss << runtime/ht << " d";

	valstring = ss.str();
	
	return valstring;

}