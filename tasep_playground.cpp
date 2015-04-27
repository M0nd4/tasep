#include <fstream>
#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <assert.h>

#include <tclap/CmdLine.h>

#include "lattice.hpp" 
#include <sys/time.h>

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
    CmdLine cmd ("run experiment, and save result to file");
    // describe arguments
    ValueArg<string> cmdRatesPath ("", "rates", "path of rates", true, "", "string", cmd);
    ValueArg<string> cmdOutPath ("", "out", "path of output", false, "/dev/null", "string", cmd);
    ValueArg<double> cmdEpoch ("", "epoch", "time to finish simulation", false, 10, "double", cmd);
    SwitchArg        cmdUseMultiple ("m", "multiple", "run multiple polysomes", cmd);
    MultiSwitchArg   cmdVerbose ("v", "verbose", "verbosity level: "
                                 "0: none, 1: O(const), 2: O(RNAs), 3: O(RNAs*epoch)", cmd);
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
        
        
        struct timeval tv1, tv2;
        gettimeofday(&tv1, NULL);     
        runMultiplePolysomes (rates, epoch, probs, verbose);
        gettimeofday(&tv2, NULL);
        printf("-----TOTAL MULTIPLE-----: Time taken in execution = %f seconds\n",
              (double) (tv2.tv_usec - tv1.tv_usec) / (double)1000000 +
              (double) (tv2.tv_sec - tv1.tv_sec));
              
              
        assert (rates.size() == probs.size());
        for (int rna = 0; rna != probs.size(); ++rna)
            for (int i = 0; i != probs[rna].size(); ++i)
                ofs << probs[rna][i] << (i == probs[rna].size()-1 ? "\n" : " ") << flush; 
    }
    else
    {
        cout << "running runSinglePolysome" << endl;
        vector<double> probs;
        
        
        struct timeval tv1, tv2;
        gettimeofday(&tv1, NULL);
        runSinglePolysome (rates[0], epoch, probs, verbose);
        gettimeofday(&tv2, NULL);
        printf("-----TOTAL SINGLE-----: Time taken in execution = %f seconds\n",
              (double) (tv2.tv_usec - tv1.tv_usec) / (double)1000000 +
              (double) (tv2.tv_sec - tv1.tv_sec));
              
              
        for (int i = 0; i != probs.size(); ++i)
            ofs << probs[i] << " " << flush;
    }
    ofs.close();
    
    return 0;
}



