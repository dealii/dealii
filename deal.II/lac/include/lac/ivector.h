// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_ivector_h
#define __lac_ivector_h
#ifndef __base_types_h
#include <base/types.h>
#endif
//#ifndef __base_errors_h
//#include <deal/errors.h>
//#endif

/*
CLASS
   iVector
   */
class iVector
{
  friend class dFMatrix;
protected:
  //////////
  int dim, maxdim;
  //////////
  int *val;
public:
  //////////
  iVector();
  //////////
  iVector(const iVector& v);
  //////////
  iVector(int n);
  //////////
  ~iVector();
  //////////
  void reinit(int n, int fast = 0);
  //////////
  void reinit(const iVector&, int fast = 0);

  //////////
  int n() const; // Abfrage der Dimension

  //////////
  int operator()(int i) const; //read-only Zugriff
  //////////
  int& operator()(int i); //Zugriff auf die Komponenten

  //////////
  iVector& operator=(int i);
  //////////
  iVector& operator=(const iVector& v);

  //////////
  void add(const iVector&);
  //////////
  void add(int, const iVector&);

  // Zuweisung

  //////////
  void equ(int, const iVector&);
};

inline int iVector::n() const
{
  return dim;
}

inline int iVector::operator() (int i) const
{
  THROW2( (i<0) || (i>=dim), IntError(IntError::Range,i));
  return val[i];
}

inline int& iVector::operator() (int i)
{
  THROW2( (i<0) || (i>=dim), IntError(IntError::Range,i));
  return val[i];
}
#endif
