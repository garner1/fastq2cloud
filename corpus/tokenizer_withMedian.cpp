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
  ifstream model10;
  ifstream model12;
  ifstream model14;
  path = string(argv[2])+"merged_10.jf.csv";
  model10.open( path ); // the MC model
  if (!model10.good()) 
    return 1;                         // exit if file not found
  path = string(argv[2])+"merged_12.jf.csv";
  model12.open( path ); // the MC model
  if (!model12.good()) 
    return 1;                         // exit if file not found
  path = string(argv[2])+"merged_14.jf.csv";
  model14.open( path ); // the MC model
  if (!model14.good()) 
    return 1;                         // exit if file not found

  map<string, float> m10,m12,m14;
  map<string, float>::iterator p10,p12,p14;
  string read,kmer;
  float value;
  while (!model10.eof())	// loop over model lines
    {
      model10 >> kmer >> value;
      m10[kmer] = value;
    }
  while (!model12.eof())	// loop over model lines
    {
      model12 >> kmer >> value;
      m12[kmer] = value;
    }
  while (!model14.eof())	// loop over model lines
    {
      model14 >> kmer >> value;
      m14[kmer] = value;
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
      vector<float> profile4(read.size(),0.0),profile6(read.size(),0.0),profile8(read.size(),0.0),profile10(read.size(),0.0),profile12(read.size(),0.0),profile14(read.size(),0.0);
      
      cout << read << endl;

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

      k = 12;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p12 = m12.find(kmer);
      	  if (p12 != m12.end())	// kmer is present as a key
      	    {
      	      info = p12->second;  
      	      profile12[i+k-1] += info;
      	    }
      	}
      median = median_func(profile12); // evaluate the median
      medianvec.push_back( median ); 

      k = 14;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p14 = m14.find(kmer);
      	  if (p14 != m14.end())	// kmer is present as a key
      	    {
      	      info = p14->second;  
      	      profile14[i+k-1] += info;
      	    }
      	}
      median = median_func(profile14); // evaluate the median
      medianvec.push_back( median ); 

      // Find the optimal kmer size
      auto smallest = min_element( begin(medianvec), end(medianvec) );
      location = distance( begin(medianvec), smallest );
      k = 2*location + 10;

      // Select the optimal profile for the read
      if ( k == 10 )
      	selection = profile10;
      if ( k == 12 )
      	selection = profile12;
      if ( k == 14 )
      	selection = profile14;

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




// cout << "The optimal kmer is fot k= " << k << endl;
// Show the map content
// cout << m.find(key)->second;
// for(map<string, float>::iterator it = m.begin(); it != m.end(); ++it)
//   {
//     cout << it->first << " kmer " << it->second << " value" << endl;
//   }
// // Show the information profile vector
// for(unsigned int i = 0; i < profile.size(); i++)
// 	{
// 	  cout << profile[i] << endl;
// 	}
