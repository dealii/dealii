// $Id$

// This file was once part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier
// Revised and extended by Wolfgang Bangerth, 1998, 1999

#include <lac/block_vector.h>
#include <cmath>
#include <algorithm>


template <int n_blocks, typename Number>
static inline Number sqr (const Number x)
{
  return x*x;
};



template <int n_blocks, typename Number>
BlockVector<n_blocks,Number>::BlockVector ()
{
  for (unsigned int i=0;i<=n_blocks;++i)
    start[n_blocks] = 0;
}



template <int n_blocks, typename Number>
BlockVector<n_blocks,Number>::BlockVector (const vector<unsigned int>& n)
{
  Assert(n.size()==n_blocks, ExcDimensionMismatch(n.size(), n_blocks));

  reinit (n, false);
}



template <int n_blocks, typename Number>
BlockVector<n_blocks,Number>::BlockVector (const BlockVector<n_blocks,Number>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      start[i] = v.start[i];
      components[i] = v.components[i];
    }
  start[n_blocks] = v.start[n_blocks];
}



// see the .h file for why this function was disabled
//
// template <int n_blocks, typename Number>
// template <typename OtherNumber>
// BlockVector<n_blocks,Number>::Vector (const BlockVector<n_blocks, OtherNumber>& v) :
// 		dim(v.size()),
// 		maxdim(v.size()),
// 		val(0)
// {
//   if (dim)
//     {
//       val = new Number[maxdim];
//       Assert (val != 0, ExcOutOfMemory());
//       copy (v.begin(), v.end(), begin());
//     }
// }



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::reinit (const vector<unsigned int>& n,
					    const bool fast)
{
  Assert (n.size() == n_blocks,
	  ExcDimensionsDontMatch(n.size(), n_blocks));
  
  unsigned int j=0;
  start[0] = j;
  for (unsigned int i=0; i<n_blocks; ++i)
    {
      components[i].reinit(n[i], fast);
      j += n[i];
      start[i+1] = j;
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::reinit (const BlockVector<n_blocks,Number>& v,
					    const bool fast)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].reinit(v.components[i], fast);
      start[i] = v.start[i];
    }
  start[n_blocks] = v.start[n_blocks];
}




template <int n_blocks, typename Number>
BlockVector<n_blocks,Number>::~BlockVector ()
{
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::clear ()
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].clear();
    }
}



template <int n_blocks, typename Number>
bool BlockVector<n_blocks,Number>::all_zero () const
{
  bool result = true;
  for (unsigned int i=0;i<n_blocks;++i)
    {
      result = result && components[i].all_zero();
    }
  return result;
}



template <int n_blocks, typename Number>
Number BlockVector<n_blocks,Number>::operator * (const BlockVector<n_blocks,Number>& v) const
{
  Number sum = 0.;
  for (unsigned int i=0;i<n_blocks;++i)
    {
      sum += components[i]*v.components[i];
    }
  return sum;
}



template <int n_blocks, typename Number>
Number BlockVector<n_blocks,Number>::norm_sqr () const
{
  Number sum = 0.;
  for (unsigned int i=0;i<n_blocks;++i)
    {
      sum += components[i].norm_sqr();
    }
  return sum;
};



template <int n_blocks, typename Number>
Number BlockVector<n_blocks,Number>::mean_value () const
{
  Number sum = 0.;
  for (unsigned int i=0;i<n_blocks;++i)
    {
      sum += components[i].mean_value() * components[i].size();
    }
  return sum/size();
}



template <int n_blocks, typename Number>
Number BlockVector<n_blocks,Number>::l1_norm () const
{
  Number sum = 0.;
  for (unsigned int i=0;i<n_blocks;++i)
    {
      sum += components[i].l1_norm();
    }
  return sum;
}



template <int n_blocks, typename Number>
Number BlockVector<n_blocks,Number>::l2_norm () const
{
  return sqrt(norm_sqr());
}



template <int n_blocks, typename Number>
Number BlockVector<n_blocks,Number>::linfty_norm () const
{
  Number sum = 0.;
  for (unsigned int i=0;i<n_blocks;++i)
    {
      Number newval = components[i].linfty_norm();
      if (sum<newval)
	sum = newval;
    }
  return sum;
}
  




template <int n_blocks, typename Number>
BlockVector<n_blocks,Number>& BlockVector<n_blocks,Number>::operator += (const BlockVector<n_blocks,Number>& v)
{
  add (v);
  return *this;
}



template <int n_blocks, typename Number>
BlockVector<n_blocks,Number>& BlockVector<n_blocks,Number>::operator -= (const BlockVector<n_blocks,Number>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
       components[i] -= v.components[i];
    }
  return *this;
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::add (const Number v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].add(v);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::add (const BlockVector<n_blocks,Number>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].add(v.components[i]);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::add (const Number a,
					 const BlockVector<n_blocks,Number>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].add(a, v.components[i]);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::add (const Number a,
					 const BlockVector<n_blocks,Number>& v,
					 const Number b,
					 const BlockVector<n_blocks,Number>& w)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].add(a, v.components[i], b, w.components[i]);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::sadd (const Number x,
					  const BlockVector<n_blocks,Number>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].sadd(x, v.components[i]);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::sadd (const Number x, const Number a,
					  const BlockVector<n_blocks,Number>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].sadd(x, a, v.components[i]);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::sadd (const Number x, const Number a,
					  const BlockVector<n_blocks,Number>& v,
					  const Number b,
					  const BlockVector<n_blocks,Number>& w)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].sadd(x, a, v.components[i], b, w.components[i]);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::sadd (const Number x, const Number a,
					  const BlockVector<n_blocks,Number>& v,
					  const Number b,
					  const BlockVector<n_blocks,Number>& w,
					  const Number c,
					  const BlockVector<n_blocks,Number>& y)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].sadd(x, a, v.components[i],
			 b, w.components[i], c, y.components[i]);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::scale (const Number factor)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].scale(factor);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::equ (const Number a,
					 const BlockVector<n_blocks,Number>& v,
					 const Number b,
					 const BlockVector<n_blocks,Number>& w)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].equ( a, v.components[i], b, w.components[i]);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::equ (const Number a,
					 const BlockVector<n_blocks,Number>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].equ( a, v.components[i]);
    }
}


template <int n_blocks, typename Number>
BlockVector<n_blocks,Number>& BlockVector<n_blocks,Number>::operator = (const Number s)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i] = s;
    }
  return *this;
}



template <int n_blocks, typename Number>
BlockVector<n_blocks,Number>&
BlockVector<n_blocks,Number>::operator = (const BlockVector<n_blocks,Number>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i] = v.components[i];
    }
  return *this;
}



template <int n_blocks, typename Number>
template<typename Number2>
BlockVector<n_blocks,Number>&
BlockVector<n_blocks,Number>::operator = (const BlockVector<n_blocks, Number2>& v)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i] = v.components[i];
    }
  return *this;
}


template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::print (ostream &out,
			    unsigned int precision,
			    bool scientific,
			    bool across) const
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      if (across)
	out << 'C' << i << ':';
      else
	out << "Component " << i << endl;
      components[i].print(out, precision, scientific, across);
    }
}



template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::block_write (ostream &out) const
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].block_write(out);
    }
}


template <int n_blocks, typename Number>
void BlockVector<n_blocks,Number>::block_read (istream &in)
{
  for (unsigned int i=0;i<n_blocks;++i)
    {
      components[i].block_read(in);
    }  
}
