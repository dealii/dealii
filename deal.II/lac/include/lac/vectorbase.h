// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_vectorbase_h
#define __lac_vectorbase_h

/*
CLASS
   VectorBase
   */
class VectorBase
{
  
  public:
  virtual ~VectorBase() {}

  // Komponentenweises Eintragen

  virtual void czero(int) = 0;
  virtual void cequ(int, const VectorBase&, double, int) = 0;
  virtual void cequ(int, const VectorBase&, double, int, double, int) = 0;
  virtual void cequ(int, const VectorBase&, double, int, double, int, double, int, double, int) = 0;

  // Komponentenweise Addition

  virtual void cadd(int, const VectorBase&, double, int) = 0;
  virtual void cadd(int, const VectorBase&, double, int, double, int) = 0;
  virtual void cadd(int, const VectorBase&, double, int, double, int, double, int, double, int) = 0;

  virtual const char* name() const = 0;
};

#endif
