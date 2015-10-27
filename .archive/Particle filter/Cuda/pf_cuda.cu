/*	Author: Dexter Barrows
	Github: dbarrows.github.io
	*/

/*	Runs a particle filter on synthetic noisy data and attempts to
	reconstruct underlying true state at each time step. Note that
	this program uses gnuplot to plot the data, so an x11
	environment must be present. Also the multiplier of 1024 in the
	definition of NP below should be set to a multiple of the number
	of multiprocessors of your GPU for optimal results.

	Also, the accompanying "pf.plg" file contains the instructions
	gnuplot will use. It must be present in the same directory as
	the executable generated by compiling this file.

	Compile with:

	nvcc -arch=sm_20 -O2 pf_cuda.cu timer.cpp rand.cpp -o pf_cuda.x

	*/

#include <cuda.h>
#include <iostream>
#include <fstream>
#include <curand.h>
#include <curand_kernel.h>
#include <string>
#include <sstream>
#include "timer.h"
#include "rand.h"

#define NP 		98*1024 		// number of particles, should optimally be a multiple of 1024*NMP,
							// where NMP is the number of multiprocessors on the cuda device
#define T 		100			// time to simulate over
#define R0 		3.0			// infectiousness
#define r 		1e-1		// recovery rate
#define N 		500 		// population size
#define B 		R0*r/N		// transmission factor
#define merr 	20  		// measurement error

#define PI 		3.141592654f
#define RMAX    4294967296

using namespace std;

// Wrapper for CUDA calls, from CUDA API
// Modified to also print the error code and string
# define CUDA_CALL(x) do { if ((x) != cudaSuccess ) {					\
	cout << " Error at " << __FILE__ << ":" << __LINE__ << endl;		\
	cout << " Error was " << x << " " << cudaGetErrorString(x) << endl;	\
	return EXIT_FAILURE ;}} while (0)									\

__device__ float 		y_est[NP*T]; 		// each particle's estimate at each step, to be reduced at the end
__device__ float 		y_save[T];			// reduced estimate for each time step
__device__ float 		d_y_noise[T]; 		// device copy of noisy data
__device__ float 		Y[3*NP]; 			// each particle's state, should be NP*L, where L is the length of the state vector
__device__ float 		Y_temp[3*NP]; 		// temp storage between first and second resampling stages
__device__ float 		w[NP];				// each particle's weight
__device__ curandState 	curand_states[NP]; 	// each particle's PRNG state


__device__ __host__ void exp_euler_SIR (float h, float t0, float tn, float * y);
string getHRmemsize (size_t memsize);
string getHRtime (double runtime);


/* 	Initialize all PRNG states, get starting state vector using initial distribution
	*/
__global__ void initializeParticles (float i_infec) {

	// global id
	int id 	= blockIdx.x*blockDim.x + threadIdx.x;

	// initialize PRNG state
	curandState state;
	curand_init(id, 0, 0, &state);

	float i_infec_pert = i_infec + merr*curand_normal(&state);
	if (i_infec_pert < 0)
		i_infec_pert = 0;
	float i_sus_pert = N - i_infec_pert;

	// save initial infected states
	y_est[id] = i_infec_pert;

	// save SIR state
	Y_temp[3*id] 	= i_sus_pert;
	Y_temp[3*id+1] 	= i_infec_pert;
	Y_temp[3*id+2] 	= 0;

	// save PRNG state
	curand_states[id] = state;

}


/* 	Project particles forward, perturb, and save weight based on data
	int t - time step number (1,...,T)
	*/
__global__ void project (int t) {

	// global id
	int id 	= blockIdx.x*blockDim.x + threadIdx.x;

	// retrieve PRNG state
	curandState state = curand_states[id];

	float y[3];

	// retrieve SIR state
	y[0] = Y_temp[3*id];
	y[1] = Y_temp[3*id+1];
	y[2] = Y_temp[3*id+2];

	exp_euler_SIR(1.0/100, 0.0, 1.0, y); 							// project particle forward, could take fewer steps
	float y_par_noise = y[1] + (float) merr*curand_normal(&state);	// perturb with expected measurement noise
	if (y_par_noise < 0)											// make sure we don't go negative
		y_par_noise = 0;
	float y_diff = d_y_noise[t] - y_par_noise;

	// get weight and save
	w[id] = 1.0/(merr*sqrt(2.0*PI)) * expf( - y_diff*y_diff / (2.0*merr*merr) );

	// save SIR state
	Y[3*id]   = y[0];
	Y[3*id+1] = y[1];
	Y[3*id+2] = y[2];

	// save PRNG state
	curand_states[id] = state;

}


/* 	The 0th thread will perform cumulative sum on the weights.
	There may be a faster way to do this, will investigate.
	*/
__global__ void cumsumWeights () {

	// global id
	int id 	= blockIdx.x*blockDim.x + threadIdx.x;

	if (id == 0) {
		for (int n = 1; n < NP; n++) {
			w[n] += w[n-1];
		}
	}

}


/* 	Resample from all particle states
	int t - time step number (1,...,T)
	*/
__global__ void resample (int t) {

	// global id
	int id 	= blockIdx.x*blockDim.x + threadIdx.x;

	// retrieve PRNG state
	curandState state = curand_states[id];

	// resampling proportional to weights
	float w_r = curand_uniform(&state) * w[NP-1];
	int i = 0;
	while (w_r > w[i]) {
		i++;
	}

	// index i was the index to resample from
	Y_temp[3*id] 	= Y[3*i];
	Y_temp[3*id+1] 	= Y[3*i+1];
	Y_temp[3*id+2] 	= Y[3*i+2];

	// save all true state estimates for later reduction
	y_est[NP*t+id] = Y[3*i+1];

	// save PRNG state
	curand_states[id] = state;

}


/* 	Reduce all the estimates made at each time step into a sequence of mean estimates
	Should by launched with T threads
	*/
__global__ void reduce () {

	// global id
	int id 	= blockIdx.x*blockDim.x + threadIdx.x;

	if (id < T) {

		int t0 = id*NP;
		float sum = 0.0;

		for (int n = 0; n < NP; n++)
			sum += y_est[t0+n];

		y_save[id] = sum / (NP);

	}

}


int main (int argc, char *argv[]) {

	float i_infec = 5;

	float y[3] = {N - i_infec, i_infec, 0}; 	// 495 susceptible, 5 infected, 0 recovered

	unsigned int rootseed = rand();

	cout << "-----------------" << endl;
	cout << "SYSTEM PARAMETERS" << endl;
	cout << endl;
	cout << "R0  " << R0 << endl;
	cout << "r   " << r << endl;
	cout << "B   " << B << endl;
	cout << "-----------------" << endl;

	printf("Running with %d particles\n", NP);

	double restime;
	struct timeval  tdr0, tdr1;

	float y_true[T];		// true number of infected peeps
	float y_noise[T];		// true number of infected peeps with observation noise

	// Generate our synthetic true trajectory and noisy observation data
	// Note that we could instead use real data, but synthetic data is more flexible for the purposes of this implementation
	y_true[0] = y[1];
	y_noise[0] = y[1] + merr*randn();
	if (y_noise[0] < 0)
		y_noise[0] = 0;
	for (int i = 1; i < T; i++) {
		exp_euler_SIR( 1.0/5, 0.0, 1.0, y);
		y_true[i] = y[1];
		y_noise[i] = y[1] + merr*randn();
		if (y_noise[i] < 0)
			y_noise[i] = 0;
	}

	size_t memsize_y = T * sizeof(float);

	// copy noisy data to device
	CUDA_CALL( cudaMemcpyToSymbol(d_y_noise, y_noise, memsize_y, 0, cudaMemcpyHostToDevice) );
	CUDA_CALL( cudaDeviceSynchronize() );

	cout << "Launching kernels ..."  << endl;

	gettimeofday (&tdr0, NULL);

	int BLOCK_DIM = 1024; // should be the maximum number of threads per multiprocessor
	int NUM_BLOCKS = NP / BLOCK_DIM + (NP % BLOCK_DIM);

	initializeParticles <<< NUM_BLOCKS, BLOCK_DIM >>> (i_infec);
	CUDA_CALL( cudaGetLastError() );
	CUDA_CALL( cudaDeviceSynchronize() );

	for (int t = 1; t < T; t++) {

		project <<< NUM_BLOCKS, BLOCK_DIM >>> (t);
		CUDA_CALL( cudaGetLastError() );
		CUDA_CALL( cudaDeviceSynchronize() );

		cumsumWeights <<< NUM_BLOCKS, BLOCK_DIM >>> ();
		CUDA_CALL( cudaGetLastError() );
		CUDA_CALL( cudaDeviceSynchronize() );

		resample <<< NUM_BLOCKS, BLOCK_DIM >>> (t);
		CUDA_CALL( cudaGetLastError() );
		CUDA_CALL( cudaDeviceSynchronize() );

	}

	int R_BLOCKS = T / BLOCK_DIM + (T % BLOCK_DIM);

	size_t avail, total;
	cudaMemGetInfo( &avail, &total );
	size_t used = total - avail;
	

	reduce <<< R_BLOCKS, BLOCK_DIM >>> ();
	CUDA_CALL( cudaGetLastError() );
	CUDA_CALL( cudaDeviceSynchronize() );

	gettimeofday (&tdr1, NULL);
    timeval_subtract (&restime, &tdr1, &tdr0);

    cout << "---------------------------" << endl;
    cout << "CUDA STATS" << endl;
    cout << endl;
    cout << "Runtime          " << getHRtime(restime) << endl;
    cout << "Device mem used  " << getHRmemsize(used) << endl;
    cout << "---------------------------" << endl;

	float h_y_save[T];
	CUDA_CALL( cudaMemcpyFromSymbol(h_y_save, y_save, T*sizeof(float), 0, cudaMemcpyDeviceToHost) );
	CUDA_CALL( cudaDeviceSynchronize() );

	string filename = "pf.dat";

	cout << "Writing results to file '" << filename << "' ..." << endl;

	ofstream outfile;
	outfile.open(filename.c_str());
	
	for (int t = 0; t < T; t++)
		outfile << t << " " << y_true[t] << " " << y_noise[t] << " " << h_y_save[t] << endl;  

	outfile.close();

	cout << "Plotting using gnuplot ..." << endl;
	cout << "Press ENTER to close plot and exit" << endl;

	string syscall("gnuplot -e \"filename='");
	syscall += filename;
	syscall += "'\" pf.plg";

	system( syscall.c_str() );

}


/*	Use the Explicit Euler integration scheme to integrate SIR model forward in time
	float h 	- time step size
	float t0 	- start time
	float tn 	- stop time
	float * y 	- current system state; a three-component vector representing [S I R], susceptible-infected-recovered
	*/
__device__ __host__ void exp_euler_SIR (float h, float t0, float tn, float * y) {

	int num_steps = floor( (tn-t0) / h );

	float S = y[0];
	float I = y[1];
	float R = y[2];

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

	y[0] = S;
	y[1] = I;
	y[2] = R;

}


/*	convert memory size in bytes to human-readable format
	*/
string getHRmemsize (size_t memsize) {

	stringstream ss;
	string valstring;

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


string getHRtime (double runtime) {

	stringstream ss;
	string valstring;

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