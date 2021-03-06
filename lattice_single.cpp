#include <deque>
#include <iostream>
#include <random>
#include <tuple>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <assert.h>
#include <sys/time.h>

#include "lattice.hpp"


using namespace std;


mt19937_64 _rg (std::random_device{}());
uniform_real_distribution<double> _rand (uniform_real_distribution<double>(0,1));


double sampleTime (double rate)
{
    // exponential  distribution
    double u(0);
    while (u==0) u = _rand(_rg);
    double dt = -std::log(u) / rate;
    assert(!std::isinf(dt));
    return dt;
}


tuple<int, double> dynamicMonteCarlo (const vector<Codon>& codons, const deque<Ribosome>& ribosomes)
{
    // fill in a vector of rates, element per ribosome
    size_t N = ribosomes.size();
    vector<double> ribosome_rates (N+1, 0);
    for (int i = 0; i != N; ++i)
        if (i == ribosomes.size()-1 || ribosomes[i+1].pos - ribosomes[i].pos > RiboWidth)
            ribosome_rates[i] = codons[ribosomes[i].pos].rate;
        else
            ribosome_rates[i] = 0;

    // cumulative sum
    vector<double> cum_interval (N+1);
    partial_sum (ribosome_rates.begin(), ribosome_rates.end(), cum_interval.begin());

    // sample (choose) a ribosome to move
    double rand_num = _rand(_rg) * cum_interval[N];
    auto lower_it = lower_bound(cum_interval.begin(), cum_interval.end(), rand_num);
    size_t ribo_id = distance(cum_interval.begin(), lower_it);

    // sample (choose) time of jump
    double dt = sampleTime (cum_interval.back());

    return make_tuple (ribo_id, dt);
}


void step (vector<Codon>& codons, deque<Ribosome>& ribosomes, double& t, double epoch)
{
    // choose ribosome and time for the jump
    int ribo_id;
    double dt;
    tie(ribo_id, dt) = dynamicMonteCarlo (codons, ribosomes);

    // correct the time
    dt = min(dt, epoch - t);

    // increase cumulative occupied time for each occupied codon
    for (const Ribosome& ribo : ribosomes)
        codons[ribo.pos].accumtime += dt;
    
    t += dt;

    // perform jump
    Ribosome& ribo = ribosomes[ribo_id];
    if (ribo.pos == codons.size()-1)
    {
        codons[ribo.pos].occupied = false;
        ribosomes.pop_back();
    }
    else
    {
        codons[ribo.pos].occupied = false;
        codons[ribo.pos+1].occupied = true;
        ++ribo.pos;
    }

    // add new ribosome if necessary
    if (ribo_id == 0)
        ribosomes.push_front( (Ribosome) { .time = t, .pos = 0 } );
}


void runSinglePolysome (const vector<double>& rates, double epoch, 
                        vector<double>& probs, int verbose)
{
    int length = rates.size();
    if (verbose > 1) cout << "length: " << length << endl;

    // init codons
    vector<Codon> codons (length);
    for (int i = 0; i != length; ++i)
        codons[i] = { .time = 0, .rate = rates[i], .occupied = false, .accumtime = 0 };

    // init ribosomes
    deque<Ribosome> ribosomes;
    ribosomes.push_front( (Ribosome) { .time = 0, .pos = 0 } );

    double t = 0;

    const int MaxIter = MaxIterMult * codons.size() * codons.size();
    int it = 0;
    for (it = 0; it != MaxIter; ++it)
    {
        // stop condition
        if (t >= epoch) break;

        if (verbose > 2)
        {
            cout << setw(2) <<  it << "   &   " << ribosomes.size()-1 << "   &    " << flush;
            // occupied
            for (int i = 0; i != length; ++i)
                cout << (codons[i].occupied ? '*' : '.');
            cout << "   &   ";
            // accumtime
            for (int i = 0; i != length; ++i)
                cout << setprecision(2) << setw(5) << codons[i].accumtime << " ";
            cout << " \\\\" << endl;
        }

        step (codons, ribosomes, t, epoch);
    }

    if (verbose > 1) cout << "finished in " << it << " iterations" << endl;
    if (it == MaxIter)
        cerr << "warning: reached the maximum number of iterations" << endl;

    // write result
    probs.resize (length);
    for (int i = 0; i != probs.size(); ++i)
        probs[i] = codons[i].accumtime / epoch;
}


void runMultiplePolysomes (const vector< vector<double> > rates, double epoch,
                           vector< vector<double> >& probs, int verbose)
{
    // time profile
    struct timeval tv1, tv2;
    gettimeofday(&tv1, NULL);

    int numRNAs = rates.size();
    probs.resize(numRNAs);

    for (int rna = 0; rna != numRNAs; ++rna)
    {
        int length = rates[rna].size();

        // init codons
        vector<Codon> codons (length);
        for (int i = 0; i != length; ++i)
            codons[i] = { .time = 0, .rate = rates[rna][i], .occupied = false, .accumtime = 0 };

        // init ribosomes
        deque<Ribosome> ribosomes;
        ribosomes.push_front( (Ribosome) { .time = 0, .pos = 0 } );

        double t = 0;

        const int MaxIter = MaxIterMult * codons.size() * codons.size();
        int it = 0;
        for (it = 0; it != MaxIter; ++it)
        {
            // stop condition
            if (t >= epoch) break;

            if (verbose > 2)
            {
                cout << setw(2) <<  it << "   &   " << ribosomes.size()-1 << "   &    " << flush;
                // occupied
                for (int i = 0; i != length; ++i)
                    cout << (codons[i].occupied ? '*' : '.');
                cout << "   &   ";
                // accumtime
                for (int i = 0; i != length; ++i)
                    cout << setprecision(1) << setw(3) << codons[i].accumtime << " ";
                cout << " \\\\" << endl;
            }

            step (codons, ribosomes, t, epoch);
        }

        if (verbose > 1)
            cout << "rna: " << rna << ", length: " << length << ", finished in " << it << " iterations" << endl;
        if (it == MaxIter)
            cerr << "warning: rna " << rna << " reached the maximum number of iterations" << endl;

        // write result
        probs[rna].resize (length);
        for (int i = 0; i != probs[rna].size(); ++i)
            probs[rna][i] = codons[i].accumtime / epoch;
    }

    gettimeofday(&tv2, NULL);
    printf("single-CPU-time-sec: %f\n",
           (double) (tv2.tv_usec - tv1.tv_usec) / (double)1000000 +
           (double) (tv2.tv_sec - tv1.tv_sec));
}

