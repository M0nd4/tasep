#include "cu_lattice.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <thrust/device_vector.h>

#include <curand.h>
#include <curand_kernel.h>

#include <iomanip>

struct Ribosome {
    double time;
    int pos;
};


struct Codon {
    double time;
    double rate;
    bool occupied;
    double accumtime;
};


__global__ static 
void setupRand ( curandState * state, unsigned long seed )
{
    int id = threadIdx.x;
    curand_init ( seed, id, 0, &state[id] );
} 


__device__ inline static
double sampleTime (double rate, curandState* globalState)
{
    curandState localState = globalState[threadIdx.x];
    double u = curand_uniform_double( &localState );
    globalState[threadIdx.x] = localState;
    if (u == 0) u += 0.1; // should almost never happen
    return -log(u) / rate;  // let rate == 0 throw
}


//__host__ __device__ inline static
__global__ static 
void updatePolysome (Codon* codons, Ribosome* ribosomes, int length, curandState* globalState )
{
    // copy ribosom data from global memory to registers
    int riboId = threadIdx.x;
    int nextId = (riboId + blockDim.x - 1) % blockDim.x;
    Ribosome ribo     = ribosomes[riboId];
    Ribosome riboNext = ribosomes[nextId];
    __syncthreads();

    // get the position of the next and previous ribosomes
    bool nextIsFar = (riboNext.pos - ribo.pos + length) % length > 1 || (ribo.pos != 0 && riboNext.pos == 0);
    if (!nextIsFar) return;

    // sample the period
    double dt = sampleTime(codons[ribo.pos].rate, globalState);

    // flip two codons and update the ribosome position
    int newpos = (ribo.pos + length + 1) % length;
    codons[newpos].occupied = true;
    codons[ribo.pos].occupied = false;
    ribosomes[riboId].pos = newpos;

    // wait for the following codon to clear (if necessary), then add dt
    double t = max(ribo.time, codons[newpos].time) + dt;
    // find the occupied time
    codons[ribo.pos].accumtime += (t - codons[ribo.pos].time);
    // update its own time and the time of the codon
    ribosomes[riboId].time = t;
    codons[ribo.pos].time = t;   // for the next ribosome
    codons[newpos].time = t;     // for computing time of occupancy by next ribosome

    // take care of the border
    if (ribo.pos == length - 1)
    {
        ribosomes[riboId].time = 0;
        codons[0].time = 0;
        codons[0].occupied = false;
    }
}

/*
__global__ static 
void computePolysome (double* rates, bool* cover, int length, int padding)
{
}
*/

using namespace std;


vector<double> runSinglePolysome (const vector<double>& rates, double initRate)
{
    int padding = 1;

    // pad the vector
    int length = rates.size();
    int lengthPadded = length + padding;

    cout << "length: " << length << endl;
    cout << "padding: " << padding << endl;

    // init codons
    thrust::device_vector<Codon> codonsVector (lengthPadded);
    for (int i = 0; i != length; ++i)
    {
        Codon codon; codon.rate = rates[i]; codon.time = 0; codon.occupied = false; codon.accumtime = 0;
        codonsVector[i+padding] = codon;
    }
    Codon codon0; codon0.rate = initRate; codon0.time = 0; codon0.occupied = false; codon0.accumtime = 0;
    codonsVector.front() = codon0;
    Codon codon1 = codonsVector[1];
    codon1.occupied = true;
    codonsVector[1] = codon1;
    Codon* deviceCodons = thrust::raw_pointer_cast( &codonsVector[0] );

    // init ribosomes
    //const int RiboWidth = 10;
    const int numRibosomes = 10;//((lengthPadded - 1) / 32 / RiboWidth + 1) * 32;
    cout << "numRibosomes: " << numRibosomes << endl;
    thrust::device_vector<Ribosome> ribosomesVector (numRibosomes);
    Ribosome ribo00; ribo00.pos = 0; ribo00.time = 0;
    Ribosome ribo10; ribo10.pos = 1; ribo10.time = 0;
    for (int i = 0; i != numRibosomes; ++i) ribosomesVector[i] = ribo00;
    ribosomesVector[0] = ribo10;
    Ribosome* deviceRibosomes = thrust::raw_pointer_cast( &ribosomesVector[0] );

    // setup seeds
    curandState* deviceStates;
    cudaMalloc ( &deviceStates, numRibosomes * sizeof(curandState) );
    setupRand <<< 1, numRibosomes >>> ( deviceStates, time(NULL) );    
        
    //cout << "iteration   codon.occupied   codon.time   ribosome.pos   ribosome.time" << endl;
    cout << "iteration   codon.occupied   codon.accumtime" << endl;
    for (int it = 0; it != 30; ++it)
    {
        cout << "iteration " << it << ": " << flush;
        for (int i = 0; i != lengthPadded; ++i)
        {
            Codon codon = codonsVector[i];
            cout << (codon.occupied ? '*' : '.');
        }
        cout << "    ";
        /*
        {
            Codon codon = codonsVector[i];
            cout << setprecision(3) << codon.time << " ";
        }
        cout << "    ";
        *//*
        for (int i = 0; i != numRibosomes; ++i)
        {
            Ribosome ribosome = ribosomesVector[i];
            cout << ribosome.pos << " ";
        }
        cout << "    ";
        *//*
        for (int i = 0; i != numRibosomes; ++i)
        {
            Ribosome ribosome = ribosomesVector[i];
            cout << setprecision(3) << ribosome.time << " ";
        }
        */
        cout << "    ";
        for (int i = 0; i != lengthPadded; ++i)
        {
            Codon codon = codonsVector[i];
            cout << setprecision(3) << codon.accumtime << " ";
        }
        
        cout << endl;

        updatePolysome 
            <<< 1, numRibosomes, 0 >>> 
            (deviceCodons, deviceRibosomes, lengthPadded, deviceStates);

    }

    const int epoch = 1;
    vector<double> probs (length);
    for (int i = 0; i != probs.size(); ++i)
    {
        Codon codon = codonsVector[i+padding];
        probs[i] = codon.accumtime / epoch;
    }

    return probs;
}
