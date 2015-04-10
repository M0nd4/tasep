#include <fstream>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include "lattice.hpp" 
#include "utils.hpp"

#include <tclap/CmdLine.h>

using namespace std;
using namespace TCLAP;

shared_ptr<vector<double>> make_rate_vec(double rfast, double rslow, size_t rlen);

int main(int argc, char **argv)
{
  // init the parser
  CmdLine cmd ("run experiment, and compare to stored file");
  // describe arguments
  ValueArg<string> cmdTruthPath ("", "truth", "path of ground truth", true, "", "string", cmd);
  ValueArg<string> cmdOutPath ("", "out", "path of output", false, "/dev/null", "string", cmd);
  MultiSwitchArg   cmdVerbose ("v", "verbose", "verbosity level", cmd);
  // parse arguments
  cmd.parse(argc, argv);
  string truthPath = cmdTruthPath.getValue();
  string outPath   = cmdOutPath.getValue();
  int    verbose   = cmdVerbose.getValue();

  shared_ptr<vector<double>> rate_vec=make_rate_vec(10, 0.1, 100);
  (*rate_vec)[0] = 10;
  polysome p(&*rate_vec);
  p.run();
  // print to a string
  ostringstream ss_result;
  ss_result <<"A: "<<p.get_Aprob()<<endl;
  ss_result <<"R: "<<p.get_Rprob()<<endl;
  // load ground truth
  ifstream ifs (truthPath);
  stringstream ss_truth;
  ss_truth << ifs.rdbuf();
  // compare the current result and the saved file. 
  // TODO: I'm expecting problems with cross-platform line endings differences
  if (ss_result.str() == ss_truth.str())
    cout << "passed" << endl;
  else
    cout << "failed" << endl;    
  // write output to file 
  ofstream ofs (outPath);
  ofs << ss_result.str() << flush;

  double* rates = &(*rate_vec)[0];

  return 0;
}

shared_ptr<vector<double>> make_rate_vec(double rfast, double rslow, size_t rlen)
{
  shared_ptr<vector<double>> rate_vec(new vector<double>(rlen,rfast));
  size_t mid(rlen/2);
  (*rate_vec)[mid] = rslow;
  return rate_vec;
}

