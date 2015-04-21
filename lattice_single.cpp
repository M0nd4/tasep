#include <deque>
#include <iostream>
#include <random>
#include <tuple>
#include <iomanip>
#include <assert.h>

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
        if (ribosomes[i].pos == codons.size() - 1 || not codons[ribosomes[i].pos+1].occupied)
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


void run (vector<Codon>& codons, deque<Ribosome>& ribosomes, double epoch, int verbose = 0)
{
}


void runSinglePolysome (const vector<double>& rates, double epoch, 
                        vector<double>& probs, int verbose)
{
    int padding = 1;

    // pad the vector
    int lengthPadded = rates.size();
    int length = rates.size() - 1;

    cout << "length: " << length << endl;
    cout << "padding: " << padding << endl;

    // init codons
    vector<Codon> codons (lengthPadded);
    for (int i = 0; i != lengthPadded; ++i)
        codons[i] = { .time = 0, .rate = rates[i], .occupied = false, .accumtime = 0 };

    // init ribosomes
    deque<Ribosome> ribosomes;
    ribosomes.push_front( (Ribosome) { .time = 0, .pos = 0 } );

    double t = 0;

    const int MaxIter = 10000 * codons.size();
    int it = 0;
    for (it = 0; it != MaxIter; ++it)
    {
        // stop condition
        if (t >= epoch) break;

        if (verbose >= 2)
        {
            cout << setw(2) <<  it << "   &   " << ribosomes.size()-1 << "   &    " << flush;
            // occupied
            for (int i = 0; i != lengthPadded; ++i)
                cout << (codons[i].occupied ? '*' : '.');
            cout << "   &   ";
            // accumtime
            for (int i = padding; i != lengthPadded; ++i)
                cout << setprecision(2) << setw(5) << codons[i].accumtime << " ";
            cout << " \\\\" << endl;
        }

        step (codons, ribosomes, t, epoch);
    }

    cout << "finished in " << it << " iterations" << endl;
    if (it == MaxIter)
        cerr << "warning: reached the maximum number of iterations" << endl;

    // write result
    probs.resize (length);
    for (int i = 0; i != probs.size(); ++i)
        probs[i] = codons[i].accumtime / epoch;
}
