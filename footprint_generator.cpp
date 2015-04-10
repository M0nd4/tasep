/**
	This file is part of the Ribomap suite
	-- an automatic pipeline for quantifying isoform-level
	ribosome profiles from ribosome profiling data.


	Copyright 2015 Hao Wang

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

	    http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
**/



#include <string>
#include <ios>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <numeric>
#include <cstdint>

#include "utils.hpp"
#include "reference_info_builder.hpp"
#include "ribomap_profiler.hpp"
#include "transcript_model.hpp" 

using namespace std;

int main(int argc, char** argv)
{
  if (argc != 9) {
    cout<< "Usage: ./footprint_generator transcript_fa cds_range sailfish_result cutoff erate_fn ilow ihigh nproc"<<endl;
    exit(1);
  }

  // cout.setf(ios::scientific);
  // cout.precision(3);
  
  // command parser
  char* ref_fa=argv[1];
  char* cds_range=argv[2]; 
  char* sf_fn=argv[3];
  char* erate_fn=argv[5];
  double ilow(std::stof(argv[6])), ihigh(std::stof(argv[7])), cutoff(std::stof(argv[4]));
  int nproc(std::stoi(argv[8]));
  //cds range
  cout<<"getting transcript info...\n";
  transcript_info tinfo(ref_fa, cds_range);
  //profile
  cout<<"constructing profile class...\n";
  ribo_profile rprofile(tinfo, sf_fn, "sailfish", cutoff);
  cout<<"number of transcripts in profile class: "<<rprofile.number_of_transcripts()<<endl;
  //model
  cout<<"initializing tasep model parameters..."<<endl;
  ribo_model ribo_rate(ref_fa, rprofile, tinfo);
  cout<<"tRNA abundance as elongation rate"<<endl;
  ribo_rate.elongation_rate_from_file(erate_fn);
  cout<<"initiation rate range: "<<ilow<<" "<<ihigh<<endl;
  ribo_rate.random_initiation_rate(ilow,ihigh,SEED);
  cout<<"compute model profile..."<<endl;
  time_point start, end;
  start = chrono::system_clock::now();
  ribo_rate.update_profile(rprofile,false,nproc);
  end = chrono::system_clock::now();
  time_period t = end-start;
  cout<<"solving the entire system uses: "<<t.count()<<" secs"<<endl;
  return 0;
}
