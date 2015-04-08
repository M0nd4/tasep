#include <fstream>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
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
  // print to a string
  ostringstream ss_result;
  ss_result <<"A: "<<p.get_Aprob()<<endl;
  ss_result <<"R: "<<p.get_Rprob()<<endl;
  // load ground truth
  ifstream ifs ("test/truth1.txt");
  stringstream ss_truth;
  ss_truth << ifs.rdbuf();
  // compare the current result and the saved file. 
  // TODO: I'm expecting problems with cross-platform line endings differences
  if (ss_result.str() == ss_truth.str())
    cout << "passed" << endl;
  else
    cout << "failed" << endl;    

  return 0;
}

shared_ptr<vector<double>> make_rate_vec(double rfast, double rslow, size_t rlen)
{
  shared_ptr<vector<double>> rate_vec(new vector<double>(rlen,rfast));
  size_t mid(rlen/2);
  (*rate_vec)[mid] = rslow;
  return rate_vec;
}

