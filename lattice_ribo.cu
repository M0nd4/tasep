#include "lattice.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

#include <thrust/device_vector.h>

#include <curand.h>
#include <curand_kernel.h>

#include <iomanip>

// FIXME: for many blocks probably need to create many states
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


__device__ static
//__global__ static 
void updatePolysome (Codon* codons, Ribosome* ribosomes, int length, double epoch, curandState* globalState)
{
    // TODO: ribosomes AoS to SoA

    // copy ribosom data from global memory to registers
    int riboId = threadIdx.x;
    int nextId = (riboId + blockDim.x - 1) % blockDim.x;
    Ribosome ribo     = ribosomes[riboId];
    Ribosome nextribo = ribosomes[nextId];

    // TODO: many ribosomes will be in pos == 0. Check early

    int pos = ribo.pos;
    int nextpos = (pos + length + (RiboWidth - RiboKeyCodon)) % length;
    
    // copy codon data from global memory to registers
    Codon codon = codons[pos];
    Codon nextcodon = codons[nextpos];

    // after copy
    __syncthreads();

    // range of covered codons, the range follows the convention [a, b)
    int beginCoveredPos = max (pos - RiboKeyCodon, 0);
    int endCoveredPos   = min (pos - RiboKeyCodon + RiboWidth, length);
    
    // update current time with time of the next codon
    double t0 = max(codon.time, nextcodon.time);
    ribosomes[riboId].time = t0;
    codons[pos].accumtime += t0 - codon.time; 
    for (int i = beginCoveredPos; i != endCoveredPos; ++i) codons[i].time = t0;

    // do not jump if can not
    bool nextIsFar = (nextribo.pos - pos > RiboWidth) || (pos != 0 && nextribo.pos == 0); 
    if (!nextIsFar) return;

    // sample the period
    double dt = sampleTime(codon.rate, globalState);
    dt = min(dt, epoch - t0);     // when about to finish
    double t = t0 + dt;

    // update times of the ribo and of the codon
    codons[pos].accumtime += dt;
    ribosomes[riboId].time = t;
    for (int i = beginCoveredPos; i != endCoveredPos; ++i) codons[i].time = t;
    codons[nextpos].time = t;     // for computing time of occupancy by next ribosome

    // finish simulation for this ribo when time reaches the epoch
    if (t >= epoch) return;

    // flip two codons and update the ribosome position
    int jumppos = (pos + 1) % length;
    codons[jumppos].occupied = true;
    codons[pos].occupied = false;
    ribosomes[riboId].pos = jumppos;

     // take care of the border
    if (pos == length - 1)
        ribosomes[riboId].time = 0;

    __syncthreads();
    if (threadIdx.x == 0)
    {
        codons[0].time = 0;
        codons[0].occupied = false;
    }
}


__device__ static
int countActiveRibos (Ribosome* ribosomes, double epoch)
{
    // TODO: rewrite with reduce
    int countNonactive = 0;
    for (int i = 0; i != blockDim.x; ++i)
    {
        Ribosome ribo = ribosomes[i];
        if (ribo.time >= epoch || ribo.pos == 0)
            ++countNonactive;
    }
    return blockDim.x - countNonactive;
}


// pass info forward and backwards
struct In {
    int maxIterMult;
    double epoch;
    int frontpadding;
    int iters4display;
};

struct Out {
    int iter;
    double* prob;
    char* occupancy;
    int* activeRibos;
};


__global__ static 
void computePolysome (Codon** codonsPtr, Ribosome** ribosomesPtr, int* lengthsPtr, 
                      In* inPtr, Out* outPtr, curandState* globalState)
{
    __shared__ int activeRibos;
   
    // each block has its own arrays (numRibosomes is the same in every block)
    Codon*    codons = codonsPtr[blockIdx.x];
    Ribosome* ribosomes = ribosomesPtr[blockIdx.x];
    int       length = lengthsPtr[blockIdx.x];
    In        in = inPtr[blockIdx.x];
    Out       out = outPtr[blockIdx.x];

    for (out.iter = 0; out.iter != in.maxIterMult * length; ++out.iter)
    {
        // stop condition
        if (threadIdx.x == 0)
            activeRibos = countActiveRibos (ribosomes, in.epoch);
        __syncthreads();
        if (activeRibos == 0) break;

        // write occupancy into specially pre-allocated memory 
        if (out.iter < in.iters4display)
        {
            if (threadIdx.x == 0) out.activeRibos[out.iter] = activeRibos;
            for (int i = threadIdx.x; i < length; i += blockDim.x)
                out.occupancy[out.iter * length + i] = codons[i].occupied;
            __syncthreads();
        }

        updatePolysome (codons, ribosomes, length, in.epoch, globalState);
    }

    // calculate resulting probability
    __syncthreads();
    
    for (int i = threadIdx.x; i < length; i += blockDim.x)
        out.prob[i] = codons[i+in.frontpadding].accumtime / in.epoch;
    
    if (threadIdx.x == 0)
        outPtr[blockIdx.x] = out;
}



using namespace std;


Codon* initCodons (const vector<double>& rates)
{
    int length = rates.size();
    vector<Codon> hostCodons (length);
    for (int i = 0; i != length; ++i)
    {
        Codon codon; codon.rate = rates[i]; codon.time = 0; codon.occupied = (i == 1); codon.accumtime = 0;
        hostCodons[i] = codon;
    }
    Codon* deviceCodons;
    cudaMalloc(&deviceCodons, length*sizeof(Codon));
    cudaMemcpy(deviceCodons, &hostCodons[0], length*sizeof(Codon), cudaMemcpyHostToDevice);
    return deviceCodons;
}

Ribosome* initRibosomes (int numRibosomes)
{
    vector<Ribosome> hostRibosomes (numRibosomes);
    Ribosome ribo00; ribo00.pos = 0; ribo00.time = 0;
    Ribosome ribo10; ribo10.pos = 1; ribo10.time = 0;
    for (int i = 0; i != numRibosomes; ++i) hostRibosomes[i] = ribo00;
    hostRibosomes[0] = ribo10;
    Ribosome* deviceRibosomes;
    cudaMalloc(&deviceRibosomes, numRibosomes*sizeof(Ribosome));
    cudaMemcpy(deviceRibosomes, &hostRibosomes[0], numRibosomes*sizeof(Ribosome), cudaMemcpyHostToDevice);
    return deviceRibosomes;
}

void initDebug (In& in, Out& out, int length, int iters4display)
{
    int space4occupancy = iters4display*length*sizeof(char);
    char* deviceOccupancy;
    int* deviceActiveRibos;
    if (iters4display) 
    {
        cudaMalloc (&deviceOccupancy, space4occupancy);
        cudaMemset (deviceOccupancy, 0, space4occupancy);
        cudaMalloc (&deviceActiveRibos, iters4display*sizeof(int));
        cudaMemset (deviceActiveRibos, 0, iters4display*sizeof(int));
    }
    in.iters4display = iters4display;
    out.occupancy = deviceOccupancy;
    out.activeRibos = deviceActiveRibos;                        
}

void printDebug (const In& in, const Out& out, int length)
{
    if (in.iters4display)
    {
        int space4occupancy = in.iters4display*length*sizeof(char);
        vector<char> vectorOccupancy (in.iters4display*length);
        cudaMemcpy (&vectorOccupancy[0], out.occupancy, space4occupancy, cudaMemcpyDeviceToHost);
        vector<int> vectorActiveRibos (in.iters4display);
        cudaMemcpy (&vectorActiveRibos[0], out.activeRibos, in.iters4display*sizeof(int), cudaMemcpyDeviceToHost);
        for (int iter = 0; iter != min(in.iters4display, out.iter); ++iter)
        {
            cout << setw(3) << iter << "  &  " << setw(3) << vectorActiveRibos[iter] << "  &  ";
            for (int i = 0; i != length; ++i)
                cout << (vectorOccupancy[iter * length + i] ? '*' : '.');
            cout << endl;
        }
        cudaFree (out.occupancy);
        cudaFree (out.activeRibos);
    }
}



void runSinglePolysome (const vector<double>& rates, double epoch, 
                        vector<double>& probs, int verbose)
{
    int length = rates.size();
    cout << "length: " << length << endl;

    // init
    const int numRibosomes = min (1024, ((length - 1) / 32 / RiboWidth + 1) * 32);
    Ribosome* deviceRibosomes = initRibosomes(numRibosomes);
    Codon*    deviceCodons    = initCodons(rates);

    // set up seeds
    curandState* deviceStates;
    cudaMalloc ( &deviceStates, numRibosomes * sizeof(curandState) );
    setupRand <<< 1, numRibosomes >>> ( deviceStates, time(NULL) );

    // pass constants
    In in; in.epoch = epoch; in.maxIterMult = 100; in.frontpadding = 1;

    // info to return
    double* deviceProb;
    cudaMalloc (&deviceProb, length*sizeof(double));
    Out out; out.prob = deviceProb;

    // debugging/visualization info to return
    int iters4display = verbose >= 2 ? 200 : 0;
    initDebug (in, out, length, iters4display);

    // copy the in/out structs to device
    thrust::device_vector<Codon*> codonsPtr       (1, deviceCodons);
    thrust::device_vector<Ribosome*> ribosomesPtr (1, deviceRibosomes);
    thrust::device_vector<int> lengthPtr          (1, length);
    thrust::device_vector<In> inPtr               (1, in);
    thrust::device_vector<Out> outPtr             (1, out);

    if (verbose) cout << "in: " << in.epoch << " " << in.maxIterMult << endl;

    computePolysome <<< 1, numRibosomes >>> (thrust::raw_pointer_cast( codonsPtr.data() ), 
                                             thrust::raw_pointer_cast( ribosomesPtr.data() ), 
                                             thrust::raw_pointer_cast( lengthPtr.data() ), 
                                             thrust::raw_pointer_cast( inPtr.data() ),
                                             thrust::raw_pointer_cast( outPtr.data() ),
                                             deviceStates);

    out = outPtr[0];
    cout << "finished in " << out.iter << " iterations" << endl;
    if (out.iter >= in.maxIterMult * length)
        cerr << "warning: reached the maximum number of iterations" << endl;

    // debugging/visualization info
    printDebug (in, out, length);

    // write result
    probs.resize(length);
    cudaMemcpy (&probs[0], deviceProb, length*sizeof(double), cudaMemcpyDeviceToHost);

    // clean up
    cudaFree (deviceProb);
    cudaFree (deviceCodons);
    cudaFree (deviceRibosomes);
    cudaFree (deviceStates);
}



// sort input array and return permutation indices. Can be done with lambdas in Cuda 7.
// TODO: do sort in a kernel
/*
struct LengthComparatorByIndex
{
    LengthComparatorByIndex (const vector<vector<double> >& data) : m_data(data) { }
    bool operator()(int left, int right) const { return m_data[left].size() < m_data[right].size(); }
    const vector< vector<double> > & m_data;
};
vector<size_t> orderedLength (const vector< vector<double> >& values) {
    vector<size_t> indices (values.size());
    for (int i = 0; i != values.size(); ++i) indices[i] = i;
    sort( indices.begin(), indices.end(), LengthComparatorByIndex(values));
    return indices;
}
*/

bool lengthCompare (const vector<double>& left, const vector<double>& right)
{
    return left.size() < right.size();
}

void runMultiplePolysomes (const vector< vector<double> > rates, double epoch,
                           vector< vector<double> >& probs, int verbose)
{
    int numRNAs = rates.size();
    probs.resize(numRNAs);

    // sort rates vectors based on length
    //vector<size_t> indices = orderedLength (rates);

    const int MaxIterMult = 100;
    
    int maxLength = max_element(rates.begin(), rates.end(), lengthCompare)->size();
    int numRibosomes = min (1024, ((maxLength - 1) / 32 / RiboWidth + 1) * 32);
    if (verbose)
        cout << "epoch: " << epoch << ", MaxIterMult: " << MaxIterMult << ", numRibos: " << numRibosomes << endl;

    // set up seeds
    curandState* deviceStates;
    cudaMalloc ( &deviceStates, numRibosomes * sizeof(curandState) );
    setupRand <<< 1, numRibosomes >>> ( deviceStates, time(NULL) );

    // one element per RNA
    thrust::device_vector<Codon*> codonsPtr       (numRNAs);
    thrust::device_vector<Ribosome*> ribosomesPtr (numRNAs);
    thrust::device_vector<int> lengthPtr          (numRNAs);
    thrust::device_vector<In> inPtr               (numRNAs);
    thrust::device_vector<Out> outPtr             (numRNAs);

    // prepare inputs
    for (int rna = 0; rna != numRNAs; ++rna)
    {
        int length = rates[rna].size();

        // init
        Ribosome* deviceRibosomes = initRibosomes(numRibosomes);
        Codon*    deviceCodons    = initCodons(rates[rna]);

        // pass constants
        In in; in.epoch = epoch; in.maxIterMult = MaxIterMult; in.frontpadding = 1;

        // info to return
        double* deviceProb;
        cudaMalloc (&deviceProb, length*sizeof(double));
        Out out; out.prob = deviceProb;

        // debugging/visualization info to return
        int iters4display = verbose >= 2 ? 200 : 0;
        initDebug (in, out, length, iters4display);

        // copy the in/out structs to device
        codonsPtr[rna]    = deviceCodons;
        ribosomesPtr[rna] = deviceRibosomes;
        lengthPtr[rna]    = length;
        inPtr[rna]        = in;
        outPtr[rna]       = out;

        if (verbose > 1) cout << "rna: " << rna << ", length: " << length << endl;
    }

    computePolysome <<< numRNAs, numRibosomes >>> (thrust::raw_pointer_cast( codonsPtr.data() ), 
                                                   thrust::raw_pointer_cast( ribosomesPtr.data() ), 
                                                   thrust::raw_pointer_cast( lengthPtr.data() ), 
                                                   thrust::raw_pointer_cast( inPtr.data() ),
                                                   thrust::raw_pointer_cast( outPtr.data() ),
                                                   deviceStates);

    // process outputs
    for (int rna = 0; rna != numRNAs; ++rna)
    {
        int length = rates[rna].size();
        Codon* deviceCodons = codonsPtr[rna];
        Ribosome* deviceRibosomes = ribosomesPtr[rna];
        Out out = outPtr[rna];
        In in = inPtr[rna];
        double* deviceProb = out.prob;

        if (verbose)
            cout << "rna: " << rna << ", finished in " << out.iter << " iterations" << endl;
        if (out.iter >= in.maxIterMult * length)
            cerr << "warning: reached the maximum number of iterations" << endl;

        // debugging/visualization info
        printDebug (in, out, length);

        // write result
        probs[rna].resize(length);
        cudaMemcpy (&probs[rna][0], deviceProb, length*sizeof(double), cudaMemcpyDeviceToHost);

        // clean up
        cudaFree (deviceProb); 
        cudaFree (deviceCodons);
        cudaFree (deviceRibosomes);
    }

    cudaFree (deviceStates);
}

