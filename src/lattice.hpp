#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <deque>
#include <vector>
#include <random>

using namespace std;
//------forward declaration------//
struct particle;
struct codon_state;
class polysome;
struct tasep_state;
class tasep_validator;

struct particle {
  size_t pos;
  double rate;
  bool sample_tagged; 
  // tag a ribosome
  // update profile estimation after this ribosome termintates
  // only one tagged ribosome in the queue at a time
};

// for storing cumulative ribosome/ribosome Asite occupancy time
struct codon_state {
  bool occupied;
  double tpre;
  double tsum;
};
ostream& operator<<(ostream& os, const codon_state& c); 

class polysome {
public:
  double t;
  double dt;
  int iteration;
  polysome(const vector<double> *rate_vec);
  size_t size() const { return _mRNA_len; } 
  void set_ribowidth(int w) { _ribosome_len=w; }
  void set_riboAsite(int A) { _Asite=A; }
  void store_state (tasep_state& state) const;
  const vector<double>& get_Aprob() const { return _Aprob; }
  const vector<double>& get_Rprob() const { return _Rprob; }
  void run();
  double compute_translation_rate() const { return _pep_cnt/t; }
  // intermediate analysis during developing
  void test_run();
  void profile_recorder(int steps, const string& fname);
  void translation_time_recorder();
  void check_run();
private:
  static mt19937_64 _rg;
  bool _tagged;
  bool _should_check;
  bool _steady;
  int _pep_cnt;
  int _ribosome_len;
  int _Asite;
  size_t _mRNA_len;
  uniform_real_distribution<double> _rand;
  // a list of bound ribosome
  // front always has the furthest location
  deque<particle> _ribosome; 
  // vector of ribosome Asite occupancy accumulative time on codons
  // mRNA_len+1 elements to make flip_codon function easier to implement
  // (flip without monitoring pos being the end of the mRNA)
  vector<codon_state> _Acover; 
  // vector of ribosome occupancy accumulative time on codons
  vector<codon_state> _Rcover;
  char _event;
  size_t _event_id;
  const vector<double> *_rate_vec;
  double _tpre; // previous timestamp of a teminated tagged ribosome
  vector<double> _Aprob;
  vector<double> _Aprob_pre;
  vector<double> _Rprob;
  vector<double> _Rprob_pre;
  void update();
  size_t jump_event();
  bool is_empty() const { return _ribosome.empty(); }
  void flip_codon(size_t pos);
  void update_Rcover(size_t pos);
  void initiate();
  void move(size_t ribo_id);
  bool should_stall(size_t ribo_id) const { return (ribo_id) ? (_ribosome.at(ribo_id-1).pos - _ribosome.at(ribo_id).pos <= _ribosome_len) : false ; }
  void update_transition_rate(size_t event_id) {_ribosome.at(event_id).rate = should_stall(event_id) ? 0 : _rate_vec->at(_ribosome[event_id].pos+1); }
  void terminate();
  void compute_profile(const char type = 'A');
  void update_profile(const char type ='A');
  bool is_steady() { return true; }
  double set_dpdt_threshold(double eps);
  bool check_steady(double eps = 1e-3);
  bool done_burn_in();
  void start_sample();
  bool done_sample() { return _steady and iteration>size()*100; }
};

struct tasep_state {
  void print_state(bool simplified = true) const;
  vector<bool> get_occupancy() const;
  bool check_dup_posi() const;
  bool check_long_toccupied(const polysome& p, double threshold) const;
  vector<codon_state> Acover;
  vector<codon_state> Rcover;
  deque<particle> ribosome;
  double t;
  double dt;
  char event;
  size_t event_id;
};

class tasep_validator {
public:
  bool check(const polysome& p);
  void show_states() const;
  bool check_single_move() const;
  tasep_state state_pre;
  tasep_state state_current;
};

bool check_steady_flow(const vector<double>& rate, const vector<double>& profile, double eps=0.1);

#endif
