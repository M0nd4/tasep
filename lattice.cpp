#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iterator>
#include "lattice.hpp"
#include "utils.hpp"
#include "math_utils.hpp"

mt19937_64 polysome::_rg(SEED);


double sample_time (uniform_real_distribution<double>& rand, 
                    mt19937_64& rg, 
                    double rate)
{
  // exponential  distribution
  double u(0);
  while (u==0) u = rand(rg);
  double dt = -std::log(u) / rate;
  assert(!std::isinf(dt));
  return dt;
}


ostream& operator<<(ostream& os, const Codon& c) {
  os<<"["<<c.occupied<<" "<<c.tpre<<" "<<c.tsum<<"]";
  return os;
}


polysome::polysome(const vector<double> *rate_vec): 
  _rate_vec  (rate_vec), 
  _mRNA_len  (rate_vec->size()-1),
  _Acover    (_mRNA_len+1, Codon{false, 0, 0}),
  _Rcover    (_mRNA_len, Codon{false, 0, 0}),
  _Aprob     (_mRNA_len, 0),
  _Aprob_pre (_mRNA_len, 0),
  _Rprob     (_mRNA_len, 0),
  _Rprob_pre (_mRNA_len, 0),
  _rand      (uniform_real_distribution<double>(0,1))
{ }


void polysome::flip_codon(size_t i)
{
  // update Asite occupancy
  if (_Acover[i].occupied)
    _Acover[i].tsum += t-_Acover[i].tpre;
  _Acover[i].occupied = !_Acover[i].occupied;
  _Acover[i].tpre = t;
}


void polysome::update_Rcover(size_t i)
{
  // update ribosome occupancy
  int ilow(i-Ribosome::Asite), ihigh(i+Ribosome::ribosome_len-Ribosome::Asite);
  // only the left most and right most position
  // covered by the ribosome need to be updated
  // update left most position
  if (ilow >= 0) {
    if ( not _Rcover[ilow].occupied ) {
      cerr<<"left most codon turned unoccupied before expected! i: "<<i<<endl;
      cout<<"mlen: "<<size()<<endl;
      exit(1);
    }
    else {
      // occupancy 1 to 0: accumulate occupancy time
      _Rcover[ilow].occupied = false;
      _Rcover[ilow].tsum += t - _Rcover[ilow].tpre;
    }
  }
  // update right most position
  if (ihigh < size()) {
    if (_Rcover[ihigh].occupied) {
      cerr<<"right most codon turned occupied before expected! i: "<<i<<endl;
      exit(1);
    }
    else {
      // occupancy 0 to 1: update time stamp
      _Rcover[ihigh].occupied = true;
      _Rcover[ihigh].tpre = t;
    }
  }
  // left corner case
  // set all occupancies in the middle
  if (i == 0) {
    for (size_t j=0; j<ihigh; ++j) {
      if ( _Rcover[j].occupied ) {
      	cerr<<"initiating ribosome before site available! j: "<<j<<endl;
      	exit(1);
      }
      else {
      	// occupancy 0 to 1: update time stamp
      	_Rcover[j].occupied = true;
      	_Rcover[j].tpre = t;
      }
    }
    return;
  }
  // right corner case
  // clear all occupancies in the middle
  // i = size() instead of size()-1 is the point where
  // the ribosome move out of the last codon location
  else if (i == size()) {
    for (size_t j=ilow+1; j!=size(); ++j) {
      if ( not _Rcover[j].occupied ) {
	      cerr<<"terminating ribosome before possible! j: "<<j<<endl;
	      exit(1);
      }
      else {
	      // occupancy 1 to 0: accumulate occupancy time
	      _Rcover[j].occupied = false;
	      _Rcover[j].tsum += t - _Rcover[j].tpre;
      }
    }
    return;
  }
}


void polysome::initiate()
{
  _ribosomes.emplace_back(Ribosome{0, _rate_vec->at(1), !_exists_tagged });
  _exists_tagged = true;
  flip_codon(0);
  update_Rcover(0);
}


void polysome::move(size_t ribo_id)
{
  _ribosomes[ribo_id].pos++;
  flip_codon(_ribosomes[ribo_id].pos);
  update_Rcover(_ribosomes[ribo_id].pos);
  flip_codon(_ribosomes[ribo_id].pos-1);
}


void polysome::terminate()
{
  if (_ribosomes.front().is_tagged) {
    update_steadiness();
    _exists_tagged = false;
  }
  _ribosomes.pop_front();
  _terminated_cnt++;
}


void polysome::update_ribosome_rate(size_t ribo_id) 
{
  // check if a ribosome can't move
  bool stalled;
  if (not ribo_id) 
    stalled = false;
  else
  {
    size_t pos_diff = _ribosomes[ribo_id-1].pos - _ribosomes[ribo_id].pos;
    stalled = (pos_diff <= Ribosome::ribosome_len);
  }

  _ribosomes.at(ribo_id).rate = (stalled ? 0 : _rate_vec->at(_ribosomes[ribo_id].pos+1));
}


size_t polysome::jump_event()
{
  // fill in a vector of rates, element per ribosome
  size_t N = _ribosomes.size();
  vector<double> ribosome_rates (N+1, 0);
  for (int i = 0; i != N; ++i)
    ribosome_rates[i] = _ribosomes[i].rate;
  if (_ribosomes.back().pos >= Ribosome::ribosome_len)
    // initiation is possible if first codon is not occupied
    ribosome_rates.back() = (*_rate_vec)[0];

  // cumulative sum
  vector<double> cum_interval (N+1);
  partial_sum (ribosome_rates.begin(), ribosome_rates.end(), cum_interval.begin());

  // sample (choose) a ribosome to move
  double rand_event = _rand(_rg) * cum_interval[N];
  auto lower_it = lower_bound(cum_interval.begin(), cum_interval.end(), rand_event);
  size_t ribo_id = distance(cum_interval.begin(), lower_it);

  // sample (choose) time of jump
  t += sample_time (_rand, _rg, cum_interval.back());
  return ribo_id;
}


void polysome::update()
{
  // step 1: choose a jump event
  // step 2: update ribosome position
  // step 3: update occupancy time
  // step 4: update rate in case of stall
  // initial initiation

  // empty
  if (is_empty()) {
    t += sample_time (_rand, _rg, (*_rate_vec)[0]);
    initiate();
    return;
  }

  size_t ribo_id = jump_event();

  // initiation
  if (ribo_id == _ribosomes.size()) {
    initiate();
    update_ribosome_rate(ribo_id);
    return;
  }

  move(ribo_id);

  // termination
  if ( _ribosomes[ribo_id].pos == size() ) {
  	terminate();
  	// front ribosome popped out
  	// new front guy's rate should be updated
  	// to unstall new guy in case it was stalled by the previous guy
  	if (not is_empty())
      update_ribosome_rate(0);
  }
  // non-termination
  else {
  	update_ribosome_rate(ribo_id);
  	if (ribo_id != _ribosomes.size()-1)
	    update_ribosome_rate(ribo_id+1);
  }
}


void polysome::compute_Aprofile()
{
  for (size_t i=0; i!=size(); ++i)
    _Aprob[i] = (_Acover[i].tsum + _Acover[i].occupied * (t - _Acover[i].tpre))/t;
}


void polysome::compute_Rprofile()
{
  for (size_t i=0; i!=size(); ++i)
    _Rprob[i] = (_Rcover[i].tsum + _Rcover[i].occupied * (t - _Rcover[i].tpre))/t;
}


void polysome::update_Aprofile()
{
  (_Aprob_pre).swap(_Aprob);
  compute_Aprofile();
}


void polysome::update_Rprofile()
{
  (_Rprob_pre).swap(_Rprob);
  compute_Rprofile();
}


void polysome::update_steadiness (double eps0)
{
  // normalize threshold on length and density of RNA
  compute_Aprofile();
  double eps = median(_Aprob) * eps0;

  update_Rprofile();
  update_Aprofile();
  double dpdt = euclidean_dist(_Aprob, _Aprob_pre)/(t - _tpre);
  if (dpdt<0) {
    cout<<"profile: ";
    for (auto pi: _Aprob)
      cout<<pi<<" ";
    cout<<"\nprofile_pre: ";
    for (auto pi: _Aprob_pre)
      cout<<pi<<" ";
    cout<<"\nt: "<<t<<" tpre: "<<_tpre<<endl;
    exit(1);
  }
  _tpre = t;
  _steady = dpdt < eps;
}


void polysome::run()
{
  if ( Ribosome::ribosome_len < Ribosome::Asite ) {
    cerr<<"A-site not within ribosome! A: "<<Ribosome::Asite<<" ribolen: "<<Ribosome::ribosome_len<<endl;
    exit(1);
  }
  /*
  // burn in
  while (_terminated_cnt <= 200) {
    update();
    iteration++;
  }
  // reset states
  _steady = false;
  iteration=0;
  */
  // collect data again
  double threshold(0);
  while (not (_steady and iteration>size()*100)) {
    update();
    iteration++;
  }
}

