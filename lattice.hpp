#include <vector>


struct Ribosome {
    double time;
    int pos;
};


struct Codon {
    double time;  // used when time is assyncronized
    double rate;
    bool occupied;
    double accumtime;
};



void runSinglePolysome (const std::vector<double>& rates, double epoch, 
                        std::vector<double>& prob,  int verbose = 0);

void runMultiplePolysome (const std::vector< std::vector<double> > rates, double epoch,
                          std::vector< std::vector<double> > prob, int verbose = 0);
