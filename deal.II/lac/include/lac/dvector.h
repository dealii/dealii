// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_dvector_h
#define __lac_dvector_h
#include <stdio.h>
#ifndef __base_types_h
#include <base/types.h>
#endif
//#ifndef __base_errors_h
//#include <base/errors.h>
//#endif
#ifndef __lac_vectorbase_h
#include <lac/vectorbase.h>
#endif

/*
CLASS
   dVector
   */
class dVector : public VectorBase
{
  friend class dFMatrix;
protected:
  int dim, maxdim;
  double *val;
public:
  dVector();
  dVector(const dVector& v);
  dVector(int n);
  ~dVector();
  void reinit(int n, int fast = 0);
  void reinit(const dVector&, int fast = 0);

  int n() const; // Abfrage der Dimension

  double operator()(int i) const; // read-only Zugriff
  double& operator()(int i); //Zugriff auf die Komponenten

  double operator*(const dVector& v) const; //Skalarprodukt

  dVector& operator=(double s);
  dVector& operator=(const dVector& v);

  // GROUP: Addition

  void add(const double);
  void add(const dVector&);
  void add(double, const dVector&);
  void add(double, const dVector&, double, const dVector&);

  void sadd(double, const dVector&);
  void sadd(double, double, const dVector&);
  void sadd(double, double, const dVector&, double, const dVector&);
  void sadd(double, double, const dVector&, double, const dVector&, 
	    double, const dVector&);

  void equ(double);
  void equ(double, const dVector&);
  void equ(double, const dVector&, double, const dVector&);

  void czero(int);
  void cequ(int, const VectorBase&, double, int);
  void cequ(int, const VectorBase&, double, int, double, int);
  void cequ(int, const VectorBase&, double, int, double, int, double, int, double, int);

  void cadd(int, const VectorBase&, double, int);
  void cadd(int, const VectorBase&, double, int, double, int);
  void cadd(int, const VectorBase&, double, int, double, int, double, int, double, int);

  virtual const char* name() const;
  //
  // Output of the vector in user-defined format.
  //
  void print(FILE* fp, const char* format = 0) const;
  void print(const char* format = 0) const;
};

inline int dVector::n() const
{
  return dim;
}

inline double dVector::operator() (int i) const
{
  THROW2( (i<0) || (i>=dim), IntError(IntError::Range,i));
  return val[i];
}

inline double& dVector::operator() (int i)
{
  THROW2( (i<0) || (i>=dim), IntError(IntError::Range,i));
  return val[i];
}

#endif
