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
    {
        assert (ribosomes[i].pos != codons.size() - 1); // virtual last codon
        if (not codons[ribosomes[i].pos+1].occupied)
            ribosome_rates[i] = codons[ribosomes[i].pos].rate;
        else
            ribosome_rates[i] = 0;
    }

    // cumulative sum
    vector<double> cum_interval (N+1);
    partial_sum (ribosome_rates.begin(), ribosome_rates.end(), cum_interval.begin());

    // all ribosomes are stalled
    if (cum_interval.back() == 0)
        return make_tuple (-1, -1);

    // sample (choose) a ribosome to move
    double rand_num = _rand(_rg) * cum_interval[N];
    auto lower_it = lower_bound(cum_interval.begin(), cum_interval.end(), rand_num);
    size_t ribo_id = distance(cum_interval.begin(), lower_it);

    // sample (choose) time of jump
    double dt = sampleTime (cum_interval.back());

    return make_tuple (ribo_id, dt);
}


void step (vector< vector<Codon> >& codons_parts, 
           vector< deque<Ribosome> >& ribosomes_parts, 
           vector< double >& t_parts, 
           double epoch, bool approximate, int verbose = 0)
{
    int parts = t_parts.size();
    assert (parts == ribosomes_parts.size());
    assert (parts == codons_parts.size());

    // the (1) goal of the last virtual codon is to reflect the next part
    for (int part = 0; part != parts - 1; ++part)
        codons_parts[part].back().occupied = codons_parts[part+1].front().occupied;

    for (int part = 0; part != parts; ++part)
    {
        if (verbose >= 2)
            cout << (part==0 ? "parts: [" : "[") << part << (part==parts-1 ? "]\n" : "] ") << flush;

        vector<Codon>& codons (codons_parts[part]);
        deque<Ribosome>& ribosomes (ribosomes_parts[part]);

        if (ribosomes.empty())
        { 
            // set the time from the previous parts if necessary
            assert (part > 0);
            t_parts[part] = max(t_parts[part], t_parts[part-1]);
            continue;
        }

        // skip iteration if reached the epoch
        if (not approximate and t_parts[part] >= epoch)
            continue;

        // skip iteration in exact scheme if times are not in sequential order
        if (not approximate and part > 0 and part < parts-1 and t_parts[part] > t_parts[part+1])
            continue;

        // the last virtual codon must be empty
        assert (ribosomes.back().pos != codons.size() - 1);
        // ribosomes must be in sequential order
        for (int i = 0; i < ribosomes.size()-1; ++i)
            assert (ribosomes[i+1].pos > ribosomes[i].pos);

        // choose ribosome and time for the jump
        int ribo_id;
        double dt;
        tie(ribo_id, dt) = dynamicMonteCarlo (codons, ribosomes);

        // skip iteration and adjust time if all ribosomes are stalled
        if (ribo_id == -1)
        {
            assert (part < parts-1);
            t_parts[part] = t_parts[part+1];
            continue;
        }

        // reset after the (1) goal of the last virtual codon
        codons.back().occupied = false;

        // correct the time
        dt = min(dt, epoch - t_parts[part]);

        // increase cumulative occupied time for each occupied codon
        for (const Ribosome& ribo : ribosomes)
            codons[ribo.pos].accumtime += dt;
        
        t_parts[part] += dt;
        //cout << "t: " << t_parts[i] << " ";

        // perform jump
        assert (ribosomes[ribo_id].pos != codons.size()-1);
        codons[ribosomes[ribo_id].pos].occupied = false;
        codons[ribosomes[ribo_id].pos+1].occupied = true;
        ++ribosomes[ribo_id].pos;
    }

    for (int part = 0; part != parts; ++part)
    {
        vector<Codon>& codons (codons_parts[part]);
        deque<Ribosome>& ribosomes (ribosomes_parts[part]);

        // move ribosome from last virtual codon
        // ribosomes move front -> back
        if (not ribosomes.empty() && ribosomes.back().pos == codons.size() - 1)
        {
            // out of RNA
            if (part == parts-1)
            {
                assert (not ribosomes.empty());
                ribosomes.pop_back();

                codons.back().occupied = false;
            }
            // to the first codon of the following part
            else
            {
                ribosomes_parts[part+1].push_front( (Ribosome) { .pos = 0 } );

                assert (not ribosomes.empty());
                ribosomes.pop_back();

                codons.back().occupied = false;
                codons_parts[part+1].front().occupied = true;

                // this hack uses codons.back().accumtime as temp storage of dt
                double adj_dt = codons.back().accumtime + t_parts[part + 1] - t_parts[part];
                codons_parts[part+1].front().accumtime += adj_dt;
                codons.back().accumtime = 0;
            }
        }
    }

    // add new ribosome if necessary
    assert (not ribosomes_parts[0].empty());
    if (ribosomes_parts[0].front().pos != 0)
        ribosomes_parts[0].emplace_front( (Ribosome) { .pos = 0 } );
}


std::vector<double> runSinglePolysome (const std::vector<double>& rates, double initRate, 
                                       double epoch, int verbose)
{
    /* initialize random seed: */
    srand (time(NULL));

    cout << "size: " << rates.size() << endl;

    const int numparts = 32;   // TODO: move to arguments (requires changes all files)
    const bool approximate = false;

    // split codons into parts
    vector< vector<Codon> >    codons_parts (numparts);
    vector< deque<Ribosome> >  ribosomes_parts (numparts);
    vector< double >           t_parts (numparts, 0);

    Codon codon = (Codon) { .time = 0, .rate = 0, .occupied = false, .accumtime = 0 };

    for (int i = 0; i != numparts; ++i)
    {
        int pos1 = rates.size() * i     / numparts;
        int pos2 = rates.size() * (i+1) / numparts;
        // logic on padding with virtual codons here
        int length = i == 0 ? (pos2 - pos1 + 2) : (pos2 - pos1 + 1);

        // init codons
        vector<Codon> codons (length, codon);
        int padfront = (i == 0 ? 1 : 0);  // logic on padding the front of the first part
        for (int j = 0; j != pos2 - pos1; ++j)
            codons[j + padfront].rate = rates[pos1+j];

        codons_parts[i] = codons;

        if (verbose)
            cout << "part " << i <<  ", length: " << length << endl;
    }

    // init the very first codon
    codons_parts[0][0].rate = initRate;

    // init the very first ribosome
    ribosomes_parts[0].push_front( (Ribosome) { .time = 0, .pos = 0 } );

    double t = 0;

    if (verbose >= 2)
    {
        cout << "rates:" << endl; 
        for (int i = 0; i != numparts; ++i)
            for (int j = 0; j != codons_parts[i].size(); ++j)
                cout << codons_parts[i][j].rate << (j == codons_parts[i].size() - 1 ? "\n" : " ");
    }

    const int MaxIter = 100 * epoch * rates.size();
    for (int it = 0; it != MaxIter; ++it)
    {
        // stop condition
        double t_mean = accumulate( t_parts.begin(), t_parts.end(), 0.0 ) / t_parts.size();
        if (t_mean >= epoch) break;

        if (verbose)
        {
            cout << setw(3) << it << "  &  " << flush;
            // occupied
            for (int i = 0; i != numparts; ++i)
                for (int j = (i == 0 ? 1 : 0); j != codons_parts[i].size() - 1; ++j)
                    cout << (codons_parts[i][j].occupied ? '*' : '.') 
                         << (j == codons_parts[i].size() - 2 ? " " : "");
            cout << "  &  " << setprecision(2);
            for (int i = 0; i != numparts; ++i)
                cout << t_parts[i] << " ";
            cout << " \\\\" << endl;
        }

        step (codons_parts, ribosomes_parts, t_parts, epoch, approximate, verbose);
    }

    // final result
    vector<double> probs;
    probs.reserve(rates.size());
    for (int i = 0; i != numparts; ++i)
        for (int j = (i == 0 ? 1 : 0); j != codons_parts[i].size() - 1; ++j)
            probs.push_back (codons_parts[i][j].accumtime / epoch);
    assert (probs.size() == rates.size());

    return probs;
}
