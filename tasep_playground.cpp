#include <fstream>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <tclap/CmdLine.h>

#include "lattice.hpp" 

using namespace std;
using namespace TCLAP;


vector<double> loadRates (const string& ratesPath)
{
    vector<double> rates;
    ifstream ifs (ratesPath.c_str());
    istringstream iss;
    string s;
    double d;
    ifs >> s;
    while (ifs)
    {
        iss.str(s);
        iss >> d;
        rates.push_back(d);
        iss.clear();
        ifs >> s;
    }
    ifs.close();
    return rates;
}


int main(int argc, char **argv)
{
    // init the parser
    CmdLine cmd ("run experiment, and compare to stored file");
    // describe arguments
    ValueArg<string> cmdRatesPath ("", "rates", "path of rates", true, "", "string", cmd);
    ValueArg<string> cmdTruthPath ("", "truth", "path of ground truth", false, "", "string", cmd);
    ValueArg<string> cmdOutPath ("", "out", "path of output", false, "/dev/null", "string", cmd);
    ValueArg<double> cmdEpoch ("", "epoch", "time to finish simulation", false, 10, "double", cmd);
    MultiSwitchArg   cmdVerbose ("v", "verbose", "verbosity level", cmd);
    // parse arguments
    cmd.parse(argc, argv);
    string ratesPath = cmdRatesPath.getValue();
    string truthPath = cmdTruthPath.getValue();
    string outPath   = cmdOutPath.getValue();
    double epoch     = cmdEpoch.getValue();
    int    verbose   = cmdVerbose.getValue();

    vector<double> rates = loadRates (ratesPath);
    if (verbose >= 2)
    {
        cout << "rates: ";
        for (int i = 0; i != rates.size(); ++i)
            cout << setprecision(3) << rates[i] << " ";
        cout << endl;
    }

    vector<double> probs;
    runSinglePolysome (rates, epoch, probs, verbose);

    // write results
    ofstream ofs (outPath.c_str());
    ostringstream oss;
    for (int i = 0; i != probs.size(); ++i)
    {
        ofs << probs[i] << endl;
        oss << probs[i] << endl;
    }
    ofs.close();

    // load ground truth
    if (truthPath != "")
    {
        ifstream ifs (truthPath);
        stringstream ss_truth;
        ss_truth << ifs.rdbuf();
        // compare the current result and the saved file. 
        // TODO: I'm expecting problems with cross-platform line endings differences
        if (oss.str() == ss_truth.str())
            cout << "passed" << endl;
        else
            cout << "failed" << endl;    
        // write output to file 
        ofstream ofs (outPath);
        ofs << oss.str() << flush;
    }

    return 0;
}



