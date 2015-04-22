#include <fstream>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <assert.h>

#include <tclap/CmdLine.h>

#include "lattice.hpp" 

using namespace std;
using namespace TCLAP;


// each line in the ratesPath file is a vector of rates
// number of lines -> number of vectors
// 
vector< vector<double> > loadRates (const string& ratesPath)
{
    vector< vector<double> > rates;

    ifstream ifs (ratesPath.c_str());
    string line;
    getline (ifs, line);
    while (ifs)
    {
        istringstream issline (line);
        rates.push_back( vector<double>() );

        // process a line in the file
        istringstream iss;
        string s;
        double d;
        issline >> s;
        while (issline)
        {
            iss.str(s);
            iss >> d;
            rates.back().push_back(d);
            iss.clear();
            issline >> s;
        }

        getline (ifs, line);
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
    ValueArg<double> cmdEpoch ("", "epoch", "time to finish simulation", false, 10, "double", cmd);
    MultiSwitchArg   cmdVerbose ("v", "verbose", "verbosity level", cmd);
    SwitchArg        cmdUseMultiple ("m", "multiple", "run multiple polysomes", cmd);
    // parse arguments
    cmd.parse(argc, argv);
    string ratesPath = cmdRatesPath.getValue();
    string outPath   = cmdOutPath.getValue();
    double epoch     = cmdEpoch.getValue();
    int    verbose   = cmdVerbose.getValue();
    bool   multiple  = cmdUseMultiple.getValue();

    // read and show rates
    vector< vector<double> > rates = loadRates (ratesPath);
    if (rates.empty())
    {
        cerr << "rates are empty" << endl;
        return 1;
    }
    if (verbose > 2)
    {
        cout << "rates: " << endl;
        for (int rna = 0; rna != rates.size(); ++rna)
            for (int i = 0; i != rates[rna].size(); ++i)
                cout << setprecision(3) << rates[rna][i] << (i == rates[rna].size()-1 ? "\n" : " ");
        cout << flush;
    }

    // write results
    ofstream ofs (outPath.c_str());
    if (multiple)
    {
        cout << "running runMultiplePolysomes" << endl;
        vector< vector<double> > probs;
        runMultiplePolysomes (rates, epoch, probs, verbose);
        assert (rates.size() == probs.size());
        for (int rna = 0; rna != probs.size(); ++rna)
            for (int i = 0; i != probs[rna].size(); ++i)
                ofs << probs[rna][i] << (i == probs[rna].size()-1 ? "\n" : " ") << flush; 
    }
    else
    {
        cout << "running runSinglePolysome" << endl;
        vector<double> probs;
        runSinglePolysome (rates[0], epoch, probs, verbose);
        for (int i = 0; i != probs.size(); ++i)
            ofs << probs[i] << " " << flush;
    }
    ofs.close();
    
    return 0;
}



