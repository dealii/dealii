//----------------------------  function.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  function.cc  ---------------------------


#include <base/memory_consumption.h>


namespace MemoryConsumption 
{
  unsigned int memory_consumption (const std::vector<std::string> &v)
  {
    unsigned int mem = sizeof(v);
    for (unsigned int i=0; i<v.size(); ++i)
      mem += memory_consumption(v[i]);
    return mem;
  };
  


  unsigned int memory_consumption (const bool) 
  {
    return sizeof(bool);
  };
  
  
  
  unsigned int memory_consumption (const char)
  {
    return sizeof(char);
  };
  


  unsigned int memory_consumption (const short int) 
  {
    return sizeof(short int);
  };
  


  unsigned int memory_consumption (const short unsigned int) 
  {
    return sizeof(short unsigned int);
  };



  unsigned int memory_consumption (const int) 
  {
    return sizeof(int);
  };
  


  unsigned int memory_consumption (const unsigned int) 
  {
    return sizeof(unsigned int);
  };



  unsigned int memory_consumption (const float)
  {
    return sizeof(float);
  };



  unsigned int memory_consumption (const double)
  {
    return sizeof(double);
  };


  
  unsigned int memory_consumption (const std::string &s)
  {
    return sizeof(s) + s.length();
  };  



  template <typename T>
  unsigned int memory_consumption (const typename std::vector<T> &v)
  {
    unsigned int mem = sizeof(std::vector<T>);
    const unsigned int n = v.size();
    for (unsigned int i=0; i<n; ++i)
      mem += memory_consumption(v[i]);
    mem += (v.capacity() - n)*sizeof(T);
    return mem;
  };



  template <typename T, int N>
  unsigned int memory_consumption (const T (&v)[N])
  {
    unsigned int mem = 0;
    for (unsigned int i=0; i<N; ++i)
      mem += memory_consumption(v[i]);
    return mem;
  };



  unsigned int memory_consumption (const std::vector<bool> &v)
  {
    return v.capacity() / 8 + sizeof(v);
  };


  
  unsigned int memory_consumption (const std::vector<int> &v)
  {
    return (v.capacity() * sizeof(int) +
	    sizeof(v));
  };
    
    

  unsigned int memory_consumption (const std::vector<double> &v)
  {
    return (v.capacity() * sizeof(double) +
	    sizeof(v));
  };
    
    

  unsigned int memory_consumption (const std::vector<float> &v)
  {
    return (v.capacity() * sizeof(float) +
	    sizeof(v));
  };
    
    
	
  unsigned int memory_consumption (const std::vector<char> &v)
  {
    return (v.capacity() * sizeof(char) +
	    sizeof(v));
  };
    

    
  unsigned int memory_consumption (const std::vector<unsigned char> &v)
  {
    return (v.capacity() * sizeof(unsigned char) +
	    sizeof(v));
  };


    
  template <typename T>
  unsigned int memory_consumption (const typename std::vector<T *> &v)
  {
    return (v.capacity() * sizeof(T *) +
	    sizeof(v));
  };
    

				    
  template <typename A, typename B>
  unsigned int memory_consumption (const typename std::pair<A,B> &p)
  {
    return (memory_consumption(p.first) +
	    memory_consumption(p.second));
  };

  
		
  template <typename T>
  unsigned int
  memory_consumption (const T * const)
  {
    return sizeof(T*);
  };


		
  template <typename T>
  unsigned int
  memory_consumption (T * const)
  {
    return sizeof(T*);
  };

  
	
  unsigned int
  memory_consumption (void * const)
  {
    return sizeof(void*);
  };
    
    
	
  template <typename T>
  unsigned int
  memory_consumption (const T &t)
  {
    return t.memory_consumption();
  };
};
