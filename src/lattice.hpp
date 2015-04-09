#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <deque>
#include <vector>
#include <random>

using namespace std;


struct Ribosome {
  // ribosome lenghth: 10 codons, A-site location: 6th codon (zero-indexed)
  // params from Steve Skiena et al.'s paper
  const static int ribosome_len = 10;
  const static int Asite = 6;
  size_t pos;
  double rate;
  bool is_tagged; 
  // tag a ribosome
  // update profile estimation after this ribosome termintates
  // only one tagged ribosome in the queue at a time
};


// for storing cumulative ribosome/ribosome Asite occupancy time
struct Codon {
  bool occupied;
  double tpre;         // time when latest flip happened
  double tsum;         // ammount of time it's been occupied so far
};
ostream& operator<<(ostream& os, const Codon& c); 


class Polysome {
public:
  double t = 0;        // current time in sec
  int iteration = 0;
  Polysome(const vector<double> *rate_vec);
  size_t size() const { return _mRNA_len; } 
  const vector<double>& get_Aprob() const { return _Aprob; }
  const vector<double>& get_Rprob() const { return _Rprob; }
  void run();
  double compute_translation_rate() const { return _terminated_cnt / t; }
  
private:
  static mt19937_64 _rg;
  uniform_real_distribution<double> _rand;

  size_t _mRNA_len;
  bool is_empty() const { return _ribosomes.empty(); }

  bool _exists_tagged = false;  // a tagged ribosome is on the RNA
  bool _steady = false;         // flag set when the flow is stabilized
  int  _terminated_cnt = 0;     // number of terminated ribosomes

  // ribosomes jump in in the back of the deque na dmove towards the front
  deque<Ribosome> _ribosomes; 

  // vector of ribosome Asite occupancy accumulative time on codons
  // mRNA_len+1 elements to make flip_codon function easier to implement
  // (flip without monitoring pos being the end of the mRNA)
  vector<Codon> _Acover; 
  // vector of ribosome occupancy accumulative time on codons
  vector<Codon> _Rcover;
  const vector<double> *_rate_vec;
  double _tpre = 0; // previous timestamp of a terminated tagged ribosome
  vector<double> _Aprob, _Aprob_pre;
  vector<double> _Rprob, _Rprob_pre;

  void update();
  size_t jump_event();
  void flip_codon(size_t pos);
  void update_Rcover(size_t pos);

  void initiate();             // a new ribosome jump in
  void move(size_t ribo_id);   // exisiting ribosome moves
  void terminate();            // the last ribosome jumps off

  void update_ribosome_rate(size_t ribo_id);
  void compute_Aprofile();
  void compute_Rprofile();
  void update_Aprofile();
  void update_Rprofile();

  void update_steadiness (double eps = 1e-3);   // update _steady
};


#endif
