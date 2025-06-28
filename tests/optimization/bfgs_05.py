## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2018 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# A companion to bfgs_05.cc which minimizes the Rosenbrok function using numpy

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import rosen
from scipy.optimize import rosen_der

# see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html

# dimension of Rosenbrok function
dim = 20

x0 = np.zeros(dim + 1)
one = np.ones(dim + 1)
location = np.zeros(dim + 1)
for i in range(dim + 1):
    location[i] = (1.0 * i) / dim


def v_rosen(theta):
    return rosen(theta - location + one)


def g_rosen(theta):
    return rosen_der(theta - location + one)


x, min_val, info = fmin_l_bfgs_b(func=v_rosen, x0=x0, fprime=g_rosen, m=3, factr=10)
dx = x - location

print("{0} iterations".format(info["nit"]))
print("function value: {0}".format(min_val))
print("linf_norm =     {0}".format(np.linalg.norm(dx, ord=np.inf)))
print("Gradient noorm: {0}".format(np.linalg.norm(g_rosen(x))))
print("function calls: {0}".format(info["funcalls"]))
print("Solution:")
print(x)
