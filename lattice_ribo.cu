#include "cu_lattice.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <thrust/device_vector.h>

#include <curand.h>
#include <curand_kernel.h>

#include <iomanip>

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
void updatePolysome (Codon* codons, Ribosome* ribosomes, int length, double epoch, curandState* globalState)
{
    // copy ribosom data from global memory to registers
    int riboId = threadIdx.x;
    int nextId = (riboId + blockDim.x - 1) % blockDim.x;
    Ribosome ribo     = ribosomes[riboId];
    Ribosome nextribo = ribosomes[nextId];

    // copy codon data from global memory to registers
    int pos = ribo.pos;
    int nextpos = (pos + length + 1) % length;
    Codon codon = codons[pos];
    Codon nextcodon = codons[nextpos];
    __syncthreads();

    // update current time with time of the next codon
    double t0 = max(codon.time, nextcodon.time);
    ribosomes[riboId].time = t0;
    codons[pos].accumtime += t0 - codon.time; 
    codons[pos].time = t0;
    
    // do not jump if can not
    bool nextIsFar = (nextribo.pos - pos + length) % length > 1 || (pos != 0 && nextribo.pos == 0); 
    if (!nextIsFar) return;

   // sample the period
    double dt = sampleTime(codon.rate, globalState);
    dt = min(dt, epoch - t0);     // when about to finish
    double t = t0 + dt;

    // update times of the ribo and of the codon
    codons[pos].accumtime += dt;
    ribosomes[riboId].time = t;
    codons[pos].time = t;         // for the next ribosome
    codons[nextpos].time = t;     // for computing time of occupancy by next ribosome

    // finish simulation for this ribo when time reaches the epoch
    if (t >= epoch) return;

    // flip two codons and update the ribosome position
    codons[nextpos].occupied = true;
    codons[pos].occupied = false;
    ribosomes[riboId].pos = nextpos;

     // take care of the border
    if (pos == length - 1)
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


vector<double> runSinglePolysome (const vector<double>& rates, double initRate, double epoch, int verbose)
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

    // set up seeds
    curandState* deviceStates;
    cudaMalloc ( &deviceStates, numRibosomes * sizeof(curandState) );
    setupRand <<< 1, numRibosomes >>> ( deviceStates, time(NULL) );    
    
    const int MaxIter = 100;
    for (int it = 0; it != MaxIter; ++it)
    {
        // stop condition
        int countNonactive = 0;
        for (int i = 0; i != numRibosomes; ++i)
        {
            Ribosome ribo = ribosomesVector[i];
            if (ribo.time >= epoch || ribo.pos == 0)
                ++countNonactive;
        }
        int countActive = ribosomesVector.size() - countNonactive;
        if (countActive == 0) break;

        if (verbose)
        {
        cout << setw(2) <<  it << "  &  " << countActive << "  &   " << flush;
        for (int i = 0; i != lengthPadded; ++i)
        {
            Codon codon = codonsVector[i];
            cout << (codon.occupied ? '*' : '.');
        }
        cout << "  &  ";
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
        for (int i = padding; i != lengthPadded; ++i)
        {
            Codon codon = codonsVector[i];
            cout << setprecision(2) << setw(2) << codon.accumtime << " ";
        }
        cout << " \\\\" << endl;
        }

        updatePolysome 
            <<< 1, numRibosomes, 0 >>> 
            (deviceCodons, deviceRibosomes, lengthPadded, epoch, deviceStates);

    }

    vector<double> probs (length);
    for (int i = 0; i != probs.size(); ++i)
    {
        Codon codon = codonsVector[i+padding];
        probs[i] = codon.accumtime / epoch;
    }

    return probs;
}
