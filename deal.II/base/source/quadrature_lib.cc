/* $Id$ */

#include <fe/quadrature_lib.h>



QGauss2<1>::QGauss2 () :
		Quadrature<1> (2)
{
  static const double xpts[] = { 0.288675135, 0.71132486 };
  static const double wts[]  = { 0.5, 0.5 };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};



QGauss2x4<1>::QGauss2x4 () :
		Quadrature<1> (8)
{
  static const double G0=0.930568156,
		      G1=0.669990522,
		      G2=0.330009478,
		      G3=0.069431844;
  static const double W0=0.173927423,
		      W1=0.326072577;
  
  static const double xpts[] = { 0.5*G0, 0.5*G1, 0.5*G2, 0.5*G3,
				 0.5*G0+0.5, 0.5*G1+0.5, 0.5*G2+0.5, 0.5*G3+0.5 };
  static const double wts[]  = { 0.5*W0, 0.5*W1, 0.5*W1, 0.5*W0,
				 0.5*W0, 0.5*W1, 0.5*W1, 0.5*W0 };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};



QGauss4<1>::QGauss4 () :
		Quadrature<1> (4)
{
  static const double G0=0.930568156,
		      G1=0.669990522,
		      G2=0.330009478,
		      G3=0.069431844;
  static const double W0=0.173927423,
		      W1=0.326072577;
  
  static const double xpts[] = { G0, G1, G2, G3 };
  static const double wts[]  = { W0, W1, W1, W0 };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};



QGauss8<1>::QGauss8 () :
		Quadrature<1> (8)
{
  static const double G0=0.0198550717512321,
		      G1=0.1016667612931866,
		      G2=0.2372337950418355,
		      G3=0.4082826787521749,
		      G4=0.5917173212478251,
		      G5=0.7627662049581646,
		      G6=0.8983332387068134,
		      G7=0.9801449282487679;
  static const double W0=0.0506142681451880,
		      W1=0.1111905172266870,
		      W2=0.1568533229389435,
		      W3=0.1813418916891810;
  
  static const double xpts[] = { G0, G1, G2, G3, G4, G5, G6, G7 };
  static const double wts[]  = { W0, W1, W2, W3, W3, W2, W1, W0 };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};



QMidpoint<1>::QMidpoint () :
		Quadrature<1>(1)
{
  quadrature_points[0] = 0.5;
  weights[0] = 1.0;
};



QSimpson<1>::QSimpson () :
		Quadrature<1> (3)
{
  static const double xpts[] = { 0.0, 0.5, 1.0 };
  static const double wts[]  = { 1./6., 2./3., 1./6. };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};



QTrapez<1>::QTrapez () :
		Quadrature<1> (2)
{
  static const double xpts[] = { 0.0, 1.0 };
  static const double wts[]  = { 0.5, 0.5 };

  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<1>(xpts[i]);
      weights[i] = wts[i];
    };
};




QGauss4<2>::QGauss4 () :
		Quadrature<2> (4) {
  static const double xpts[] = { 0.211324865,   0.788675135, 
				 0.211324865,   0.788675135  };
  static const double ypts[] = { 0.211324865,   0.211324865,
				 0.788675135,   0.788675135  };
  static const double wts[]  = { 1./4., 1./4., 1./4., 1./4.  };
  
  for (unsigned int i=0; i<n_quadrature_points; ++i) 
    {
      quadrature_points[i] = Point<2>(xpts[i], ypts[i]);
      weights[i] = wts[i];
    };
  
};
