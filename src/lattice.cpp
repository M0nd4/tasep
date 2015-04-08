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

ostream& operator<<(ostream& os, const codon_state& c) {
  os<<"["<<c.occupied<<" "<<c.tpre<<" "<<c.tsum<<"]";
  return os;
}

// ribosome lenghth: 10 codons, A-site location: 6th codon (zero-indexed)
// params from Steve Skiena et al.'s paper
polysome::polysome(const vector<double> *rate_vec): t(0), iteration(0), _tagged(false), _should_check(false), _steady(false), _pep_cnt(0), _ribosome_len(10), _Asite(6), _event('i'), _event_id(0), _rate_vec(rate_vec), _tpre(0)
{
  _rand = uniform_real_distribution<double>(0,1);
  _mRNA_len = rate_vec->size()-1;
  _Acover = vector<codon_state>(_mRNA_len+1, codon_state{false, 0, 0});
  _Rcover = vector<codon_state>(_mRNA_len, codon_state{false, 0, 0});
  _Aprob = vector<double>(_mRNA_len, 0);
  _Aprob_pre = vector<double>(_mRNA_len, 0);
  _Rprob = vector<double>(_mRNA_len, 0);
  _Rprob_pre = vector<double>(_mRNA_len, 0);
}

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
  int ilow(i-_Asite), ihigh(i+_ribosome_len-_Asite);
  // only the left most and right most position
  // covered by the ribosome need to be updated
  // update left most position
  if (ilow >= 0) {
    if ( not _Rcover[ilow].occupied ) {
      cerr<<"left most codon turned unoccupied before expected! i: "<<i<<endl;
      cout<<"mlen: "<<_mRNA_len<<endl;
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
  _ribosome.emplace_back(particle{0, _rate_vec->at(1), !_tagged });
  if (not _tagged) _tagged = true;
  flip_codon(0);
  update_Rcover(0);
  _event = 'i';
}

void polysome::move(size_t ribo_id)
{
  _ribosome[ribo_id].pos++;
  flip_codon(_ribosome[ribo_id].pos);
  update_Rcover(_ribosome[ribo_id].pos);
  flip_codon(_ribosome[ribo_id].pos-1);
}

void polysome::terminate()
{
  if (_ribosome.front().sample_tagged) {
    _should_check = true;
    _tagged = false;
  }
  _ribosome.pop_front();
  _pep_cnt++;
  _event = 't';
}

size_t polysome::jump_event()
{
  vector<double> cum_interval(_ribosome.size()+1, 0);
  // compute cumulative intervals of all jumping events
  // fill out the first one
  cum_interval[0] = _ribosome[0].rate;
  // fill out the middle part
  size_t n = cum_interval.size();
  for (size_t i=0; i!=n-2; ++i) {
    cum_interval[i+1] = cum_interval[i] + _ribosome[i+1].rate;
  }
  // fill out the last one
  cum_interval.back() = cum_interval[n-2];
  // initiation possible if first cocon not occupied
  if (_ribosome.back().pos >= _ribosome_len)
    cum_interval.back() += _rate_vec->at(0);
  // sample jump event
  double rand_event = _rand(_rg)*cum_interval[n-1];
  auto lower_it = lower_bound(cum_interval.begin(), cum_interval.end(), rand_event);
  size_t event_id = distance(cum_interval.begin(), lower_it);
  // sample jump time
  double u(0);
  while (u==0) u = _rand(_rg);
  double dt = -std::log(u)/cum_interval.back();
  assert(!std::isinf(dt));
  t += dt;
  // cout<<"rate_vec: ";
  // for (auto r: _ribosome)
  //   cout<<r.rate<<" ";
  // cout<<endl;
  // cout<<"cum_vec: ";
  // for (auto c: cum_interval)
  //   cout<<c<<" ";
  // cout<<endl;
  // cout<<"rand1: "<<rand_event<<" jump event: "<<event_id<<" t "<<t<<" dt "<<dt<<" ";
  return event_id;
}

void polysome::update()
{
  // step 1: choose a jump event
  // step 2: update ribosome position
  // step 3: update occupancy time
  // step 4: update rate in case of stall
  // initial initiation
  if (is_empty()) {
    double dt = -std::log(_rand(_rg))/_rate_vec->at(0);
    t += dt;
    initiate();
    _event_id = 0;
    //cout<<"initial initiation "<<endl;
  } //if empty
  // sample jump events
  else {
    size_t ribo_id = jump_event();
    _event_id = ribo_id;
    // initiation
    if (ribo_id == _ribosome.size()) {
      initiate();
      update_transition_rate(ribo_id);
      //cout<<"initiate "<<ribo_id<<endl;
    } //if initiation
    else {
      move(ribo_id);
      // termintation
      if ( _ribosome[ribo_id].pos == _mRNA_len ) {
	terminate();
	// front ribosome popped out
	// new front guy's rate should be updated
	// to unstall new guy in case it was stalled by the previous guy
	if (not is_empty())
	  update_transition_rate(0);
	//cout<<"terminate "<<ribo_id<<endl;
      }
      // elongation
      else {
	_event = 'e';
	update_transition_rate(ribo_id);
	if (ribo_id != _ribosome.size()-1)
	  update_transition_rate(ribo_id+1);
	//cout<<"elongate "<<ribo_id<<endl;
      }// else non-termination
    }// else non-initiation
  }// else non-empty
}

void polysome::compute_profile(const char type)
{
  vector<double> *profile;
  vector<codon_state> *codon;
  if (type=='A') {
    profile = &_Aprob;
    codon = &_Acover;
  }
  else if (type=='R') {
    profile = &_Rprob;
    codon = &_Rcover;
  }
  else {
    cerr<<"profile type "<<type<<" not supported! fail to compute polysome::compute_profile()!"<<endl;
    exit(1);
  }
  for (size_t i=0; i!=size(); ++i)
    (*profile)[i] =  ((*codon)[i].tsum + (*codon)[i].occupied * (t - (*codon)[i].tpre))/t;
}

void polysome::update_profile(const char type)
{
  vector<double> *p, *ppre;
  if (type=='A') {
    p = &_Aprob;
    ppre = &_Aprob_pre;
  }
  else if (type=='R') {
    p = &_Rprob;
    ppre = &_Rprob_pre;
  }
  else {
    cerr<<"profile type "<<type<<" not supported! fail to compute polysome::update_profile()!"<<endl;
    exit(1);

  }
  (*ppre).swap(*p);
  compute_profile(type);
}

bool polysome::check_steady(double eps)
{
  update_profile('R');
  update_profile('A');
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
  _should_check = false;
  _tpre = t;
  return dpdt<eps;
}

double polysome::set_dpdt_threshold(double eps)
{
  compute_profile('A');
  double scale(median(_Aprob));
  return scale*eps;
}


bool polysome::done_burn_in()
{
  // reference: The effect of tRNA levels on decoding times of mRNA codon
  return _pep_cnt > 200;
}


void polysome::start_sample()
{
  _steady = false;
  iteration=0;
  // t = 0;
  // _tpre = 0;
  // for (auto& c: _Acover) {
  //   c.tpre = 0;
  //   c.tsum = 0;
  // }
  // for (auto& c: _Rcover) {
  //   c.tpre = 0;
  //   c.tsum = 0;
  // }
}


void polysome::run()
{
  if ( _ribosome_len < _Asite ) {
    cerr<<"A-site not within ribosome! A: "<<_Asite<<" ribolen: "<<_ribosome_len<<endl;
    exit(1);
  }
  // burn in
  while (not done_burn_in()) {
    update();
    iteration++;
  }
  // reset states
  start_sample();
  // collect data again
  double threshold(0);
  while(not done_sample()) {
    update();
    iteration++;
    if (_should_check) {
      threshold = set_dpdt_threshold(0.05);
      _steady = check_steady(threshold);
    }
  }
}

