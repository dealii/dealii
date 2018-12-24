# //-----------------------------------------------------------
# //
# //    Copyright (C) 2018 by the deal.II authors
# //
# //    This file is part of the deal.II library.
# //
# //    The deal.II library is free software; you can use it, redistribute
# //    it, and/or modify it under the terms of the GNU Lesser General
# //    Public License as published by the Free Software Foundation; either
# //    version 2.1 of the License, or (at your option) any later version.
# //    The full text of the license can be found in the file LICENSE.md at
# //    the top level directory of deal.II.
# //
# //---------------------------------------------------------------

# A companion to bfgs_05.cc which minimizes the Rosenbrok function using numpy

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import rosen
from scipy.optimize import rosen_der

# see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html

# dimension of Rosenbrok function
dim = 20

x0 = np.zeros(dim+1)
one = np.ones(dim+1)
location = np.zeros(dim+1)
for i in range(dim+1):
    location[i] = (1.*i)/dim


def v_rosen(theta):
    return rosen(theta - location + one)


def g_rosen(theta):
    return rosen_der(theta - location + one)


x, min_val, info = fmin_l_bfgs_b(func=v_rosen,x0=x0,fprime=g_rosen,m=3, factr=10)
dx = x - location

print "{0} iterations".format(info['nit'])
print "function value: {0}".format(min_val)
print "linf_norm =     {0}".format(np.linalg.norm(dx,ord=np.inf))
print "Gradient noorm: {0}".format(np.linalg.norm(g_rosen(x)))
print "function calls: {0}".format(info['funcalls'])
print "Solution:"
print x
