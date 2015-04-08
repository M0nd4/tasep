#include <fstream>
#include <vector>
#include <memory>
#include "lattice.hpp" 
#include "utils.hpp"

using namespace std;

shared_ptr<vector<double>> make_rate_vec(double rfast, double rslow, size_t rlen);

int main()
{
  shared_ptr<vector<double>> rate_vec=make_rate_vec(10, 0.1, 100);
  (*rate_vec)[0] = 10;
  polysome p(&*rate_vec);
  p.set_ribowidth(10);
  p.run();
  cout<<"A: "<<p.get_Aprob()<<endl;
  cout<<"R: "<<p.get_Rprob()<<endl;
  return 0;
}

shared_ptr<vector<double>> make_rate_vec(double rfast, double rslow, size_t rlen)
{
  shared_ptr<vector<double>> rate_vec(new vector<double>(rlen,rfast));
  size_t mid(rlen/2);
  (*rate_vec)[mid] = rslow;
  return rate_vec;
}

