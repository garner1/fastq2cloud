// This is my first C++ program

// g++ -std=c++11 first.cpp -o first1 
// ./first1 reads.txt ~/Work/dataset/fastq2cloud/test1K/MCmodel/

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <stdlib.h>
#include <map>
#include <vector>
#include <iterator>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <math.h>
#include <boost/range/adaptor/sliced.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/assign.hpp>

using namespace std;

using std::cout;
using std::endl;
using std::ifstream;

using namespace boost::accumulators;

float std_func(vector<float> v, int start, int end)
{
  double sum = accumulate(v.begin()+start,v.begin()+end,0.0);
  double mean = sum / ((end-start)*1.0);
  double sq_sum = inner_product(v.begin()+start,v.begin()+end, v.begin()+start, 0.0);
  double stdev = sqrt(sq_sum / ((end-start)*1.0) - mean * mean);

  return stdev;
}

float median_func(vector<float> newvec)
{
  float median;
  int size = newvec.size();
  sort(newvec.begin(), newvec.end());
  if (size % 2 == 0)
    {
      median = (newvec[size / 2 - 1] + newvec[size / 2]) / 2.0;
    }
  else
    median = newvec[size / 2]; 

  return median;
}

float mean_func(vector<float> v, int start, int end)
{
  double sum = accumulate(v.begin()+start,v.begin()+end,0.0);
  double mean = sum / ((end-start)*1.0);

  return mean;
}

vector<int> stop_finder( vector<float> y, int dim, int lag, float threshold, float influence )
{
  vector<int> signals(y.size());      // init the stop-words vector
  vector<float> filteredY(y.size());  // init filtering vector
  vector<float> avgFilter(y.size());  // init average vector
  vector<float> stdFilter(y.size());  // init std vector

  filteredY = y;			   // Init filteredY
  avgFilter[dim+lag] = mean_func(y,dim,dim+lag); // Init first value
  stdFilter[dim+lag] = std_func(y,dim,dim+lag); // Init first value

  for ( unsigned i = dim+lag+1; i < y.size(); i++ ) // loop over read
    {
      if ( ( y[i] - avgFilter[i-1] ) > threshold*stdFilter[i-1] && y[i-1] < y[i] )
	{
	  signals[i] = 1;
	  filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]; // update local data points
	}
      else                                            // if not stop word
	{
	  signals[i] = 0;                               
	  filteredY[i] = y[i];
	}
      avgFilter[i] = mean_func(filteredY,i-lag,i); // update the local mean
      stdFilter[i] = std_func(filteredY,i-lag,i);  // update the local std
    }
  return signals;
}

int main ( int argc, char* argv[] )
{ 
  string path;
  ifstream model6;
  ifstream model8;
  ifstream model10;

  path = string(argv[2])+"merged_6.jf.csv";
  model6.open( path ); // the MC model
  if (!model6.good()) 
    return 1;                         // exit if file not found

  path = string(argv[2])+"merged_8.jf.csv";
  model8.open( path ); // the MC model
  if (!model8.good()) 
    return 1;                         // exit if file not found

  path = string(argv[2])+"merged_10.jf.csv";
  model10.open( path ); // the MC model
  if (!model10.good()) 
    return 1;                         // exit if file not found

  map<string, float> m6,m8,m10;
  map<string, float>::iterator p6,p8,p10;
  string read,kmer;
  float value;

  while (!model6.eof())	// loop over model lines
    {
      model6 >> kmer >> value;
      m6[kmer] = value;
    }

  while (!model8.eof())	// loop over model lines
    {
      model8 >> kmer >> value;
      m8[kmer] = value;
    }

  while (!model10.eof())	// loop over model lines
    {
      model10 >> kmer >> value;
      m10[kmer] = value;
    }

  ifstream reads;
  reads.open( argv[1] );              // We assume argv[1] is the file with reads
  if (!reads.good()) 
    return 1;                         // exit if file not found

  vector<int> contextsize,signal;
  vector<float> medianvec,selection;
  float number, info, median;
  int k,location;


  while (getline(reads, read))	// loop over reads
    {
      medianvec = {};
      vector<float> profile6(read.size(),0.0),profile8(read.size(),0.0),profile10(read.size(),0.0);
      
      k = 6;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p6 = m6.find(kmer);
      	  if (p6 != m6.end())	// kmer is present as a key
      	    {
      	      info = p6->second;  
      	      profile6[i+k-1] += info;
      	    }
      	}
      median = median_func(profile6); // evaluate the median
      medianvec.push_back( median ); 

      k = 8;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p8 = m8.find(kmer);
      	  if (p8 != m8.end())	// kmer is present as a key
      	    {
      	      info = p8->second;  
      	      profile8[i+k-1] += info;
      	    }
      	}
      median = median_func(profile8); // evaluate the median
      medianvec.push_back( median ); 

      k = 10;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p10 = m10.find(kmer);
      	  if (p10 != m10.end())	// kmer is present as a key
      	    {
      	      info = p10->second;  
      	      profile10[i+k-1] += info;
      	    }
      	}
      median = median_func(profile10); // evaluate the median
      medianvec.push_back( median ); 

      // Find the optimal kmer size
      auto smallest = min_element( begin(medianvec), end(medianvec) );
      location = distance( begin(medianvec), smallest );
      k = 2*location + 6;

      // Select the optimal profile for the read
      if ( k == 6 )
      	selection = profile6;
      if ( k == 8 )
      	selection = profile8;
      if ( k == 10 )
      	selection = profile10;

      vector<float> y = selection;
      int context = k-1;
      int lag = 10;
      float threshold = 1.0;
      float influence =0.5;      
      signal = stop_finder( y, context, lag, threshold, influence );
      for( unsigned i=0; i < read.size() ; i++ )
      	{
      	  if (i == 0)
      	    cout << "['";
      	  if ( signal[i] == 0 )
      	    cout << read[i];
      	  else
      	    cout << "','";
      	  if (i == read.size() - 1 )
      	    cout << "']";
      	}
      cout << endl;
    } // end of loop over reads
} 


