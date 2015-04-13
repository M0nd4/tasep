#include "cu_lattice.hpp"

#include <stdio.h>
#include <stdlib.h>

// for one single polysom  only
__global__ static 
void computePolysome (double* deviceRates, int length)
{
	
}


void runSinglePolysome (double* rates, int length)
{
	double *deviceRates;
    cudaMalloc(&deviceRates, length * sizeof(double));
    
    cudaMemcpy(&deviceRates, rates,
               length * sizeof(double), cudaMemcpyHostToDevice);

    
    computePolysome
        <<< 1, 1, 0 >>>
        (deviceRates, length);
}
