#include <fstream>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "cu_lattice.hpp" 
#include "utils.hpp"

#include <tclap/CmdLine.h>

using namespace std;
using namespace TCLAP;


vector<double> loadRates (const string& ratesPath)
{
    vector<double> rates;
    ifstream ifs (ratesPath.c_str());
    istringstream iss;
    string s;
    double d;
    while (ifs)
    {
        ifs >> s;
        iss.str(s);
        iss >> d;
        rates.push_back(d);
        iss.clear();
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
    ValueArg<string> cmdOutPath ("", "out", "path of output", false, "/dev/null", "string", cmd);
    MultiSwitchArg   cmdVerbose ("v", "verbose", "verbosity level", cmd);
    // parse arguments
    cmd.parse(argc, argv);
    string ratesPath = cmdRatesPath.getValue();
    string outPath   = cmdOutPath.getValue();
    int    verbose   = cmdVerbose.getValue();

    vector<double> rate_vec = loadRates (ratesPath);
    if (verbose >= 2)
    {
        cout << "rates: ";
        for (int i = 0; i != rate_vec.size(); ++i)
            cout << setprecision(3) << rate_vec[i] << " ";
        cout << endl;
    }

    runSinglePolysome (*rate_vec, 1);
    return 0;
}


