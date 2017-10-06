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

float mean_func(vector<float> v, int start, int end)
{
  double sum = accumulate(v.begin()+start,v.begin()+end,0.0);
  double mean = sum / ((end-start)*1.0);

  return mean;
}

float std_func(vector<float> v, int start, int end)
{
  double sum = accumulate(v.begin()+start,v.begin()+end,0.0);
  double mean = sum / ((end-start)*1.0);
  double sq_sum = inner_product(v.begin()+start,v.begin()+end, v.begin()+start, 0.0);
  double stdev = sqrt(sq_sum / ((end-start)*1.0) - mean * mean);

  return stdev;
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

int main ( int argc, char* argv[] )
{ 
  string path;
  ifstream model4;
  ifstream model6;
  ifstream model8;
  ifstream model10;
  ifstream model12;
  ifstream model14;
  path = string(argv[2])+"_4transitionMatrix_fromFastq_4to1.csv";
  model4.open( path ); // the MC model
  if (!model4.good()) 
    return 1;                         // exit if file not found
  path = string(argv[2])+"_6transitionMatrix_fromFastq_6to1.csv";
  model6.open( path ); // the MC model
  if (!model6.good()) 
    return 1;                         // exit if file not found
  path = string(argv[2])+"_8transitionMatrix_fromFastq_8to1.csv";
  model8.open( path ); // the MC model
  if (!model8.good()) 
    return 1;                         // exit if file not found
  path = string(argv[2])+"_10transitionMatrix_fromFastq_10to1.csv";
  model10.open( path ); // the MC model
  if (!model10.good()) 
    return 1;                         // exit if file not found
  path = string(argv[2])+"_12transitionMatrix_fromFastq_12to1.csv";
  model12.open( path ); // the MC model
  if (!model12.good()) 
    return 1;                         // exit if file not found
  path = string(argv[2])+"_14transitionMatrix_fromFastq_14to1.csv";
  model14.open( path ); // the MC model
  if (!model14.good()) 
    return 1;                         // exit if file not found

  map<string, float> m4,m6,m8,m10,m12,m14;
  map<string, float>::iterator p4,p6,p8,p10,p12,p14;
  string read,kmer;
  float value;
  while (!model4.eof())	// loop over model lines
    {
      model4 >> kmer >> value;
      m4[kmer] = value;
    }
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
  vector<float> sumvec,selection;
  float number, info, sum;
  int k,location;


  while (getline(reads, read))	// loop over reads
    {
      sumvec = {};
      vector<float> profile4(read.size(),0.0),profile6(read.size(),0.0),profile8(read.size(),0.0),profile10(read.size(),0.0),profile12(read.size(),0.0),profile14(read.size(),0.0);

      k = 4;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1) // loop over read
      	{
      	  kmer = read.substr(i, k);
      	  p4 = m4.find(kmer);
      	  if (p4 != m4.end())
      	    {
      	      info = p4->second;  
      	      profile4[i+k-1] += info; // and the end of the kmer, for bi-directionality
      	    }
      	}
      sum = accumulate( profile4.begin()+(k-1), profile4.end(), 0.0 );
      sumvec.push_back( sum/(1.0*(read.size()-k+1)) ); // evaluate the numb of bits per effective read length times the numb of bit of a single kmer

      k = 6;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1) // loop over read
      	{
      	  kmer = read.substr(i, k);
      	  p6 = m6.find(kmer);
      	  if (p6 != m6.end())
      	    {
      	      info = p6->second;  
      	      profile6[i+k-1] += info; // and the end of the kmer, for bi-directionality
      	    }
      	}
      sum = accumulate( profile6.begin()+(k-1), profile6.end(), 0.0 );
      sumvec.push_back( sum/(1.0*(read.size()-k+1)) ); // evaluate the numb of bits per effective read length times the numb of bit of a single kmer

      k = 8;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1) // loop over read
      	{
      	  kmer = read.substr(i, k);
      	  p8 = m8.find(kmer);
      	  if (p8 != m8.end())
      	    {
      	      info = p8->second;  
      	      profile8[i+k-1] += info; // and the end of the kmer, for bi-directionality
      	    }
      	}
      sum = accumulate( profile8.begin()+(k-1), profile8.end(), 0.0 );
      sumvec.push_back( sum/(1.0*(read.size()-k+1)) ); // evaluate the numb of bits per effective read length times the numb of bit of a single kmer

      k = 10;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p10 = m10.find(kmer);
      	  if (p10 != m10.end())
      	    {
      	      info = p10->second;  
      	      profile10[i+k-1] += info;
      	    }
      	}
      sum = accumulate( profile10.begin()+(k-1), profile10.end(), 0.0 );
      sumvec.push_back( sum/(1.0*(read.size()-k+1)) ); // evaluate the numb of bits per effective read length times the numb of bit of a single kmer

      k = 12;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p12 = m12.find(kmer);
      	  if (p12 != m12.end())
      	    {
      	      info = p12->second;  
      	      profile12[i+k-1] += info;
      	    }
      	}
      sum = accumulate( profile12.begin()+(k-1), profile12.end(), 0.0 );
      sumvec.push_back( sum/(1.0*(read.size()-k+1)) ); // evaluate the numb of bits per effective read length times the numb of bit of a single kmer

      k = 14;
      for (unsigned i = 0; i < read.length()-(k-1); i += 1)
      	{
      	  kmer = read.substr(i, k);
      	  p14 = m14.find(kmer);
      	  if (p14 != m14.end())
      	    {
      	      info = p14->second;  
      	      profile14[i+k-1] += info;
      	    }
      	}
      sum = accumulate( profile14.begin()+(k-1), profile14.end(), 0.0 );
      sumvec.push_back( sum/(1.0*(read.size()-k+1)) ); // evaluate the numb of bits per effective read length times the numb of bit of a single kmer

      // Find the optimal kmer size
      auto smallest = min_element(begin(sumvec), end(sumvec));
      location = distance(begin(sumvec), smallest);
      k = 2*location + 4;

      // Select the optimal profile for the read
      if ( k == 4 )
      	selection = profile4;
      if ( k == 6 )
      	selection = profile6;
      if ( k == 8 )
      	selection = profile8;
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
    }// end of loop over reads
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
