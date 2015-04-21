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



// return probability
std::vector<double> runSinglePolysome (const std::vector<double>& rates, double epoch, int verbose = 0);

