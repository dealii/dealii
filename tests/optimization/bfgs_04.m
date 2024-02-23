%  ------------------------------------------------------------------------
%
%  SPDX-License-Identifier: LGPL-2.1-or-later
%  Copyright (C) 2018 by the deal.II authors
%
%  This file is part of the deal.II library.
%
%  Part of the source code is dual licensed under Apache-2.0 WITH
%  LLVM-exception OR LGPL-2.1-or-later. Detailed license information
%  governing the source code and code contributions can be found in
%  LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
%
%  ------------------------------------------------------------------------

% a companion to bfgs_04.cc which minimizes the same function using Octave.

% bfgs_04.cc matches the Octave. Strangely enough "stepsize" in Octave differs starting from
% iteration 4, however the value, gradient and increment in value (three columns) match perfectly.

% This file uses bfgsmin() from optim package: https://octave.sourceforge.io/optim/
% https://octave.sourceforge.io/optim/package_doc/bfgsmin.html#bfgsmin
pkg load optim

% dimension of Laplace
global N=20;

% 1D Laplace with zero Dirichlet BC on both sides
function M = get_m()
  global N;
  for i=1:N
    M(i,i)=2;
    if i>1
      M(i,i-1)=-1;
    end
    if i<N
      M(i,i+1)=-1;
    end
  end
endfunction

% RHS
function b = get_b()
  global N;
  b(1:N,1)=1;
endfunction

global A=get_m();
global b=get_b();

x0(1:N,1)=1;

% objective function with analytic gradient:
% f(x) = 0.5 x*Ax - x*b
% g(x) = Ax - b
function [obj_value, gradient] = func(x)
  global b;
  global A;
  Ax = A * x;
  obj_value = 0.5 * (x' * Ax) - x'*b;
  gradient  = Ax - b;
endfunction

mMax    = 4;      % maximum number of stored residuals
itmax   = 100;    % maximum allowable number of iterations
ftol    = 1e-12;  % function change tolerance
xtol    = 1e-6;   % parameter change tolerance
gtol    = 1e-5;   % gradient tolerance

verb    = 2;      % verbosity [0,3]
%
% "used analytic gradient" prints 3 columns:
% x_k  g_k  dx = alpha*p
%

control = {itmax;verb;1;1;mMax;ftol;xtol;gtol};
[x, obj_value, convergence, iters] = bfgsmin("func", {x0}, control);

fprintf("BFGS solution:\n");
for i=1:N
  fprintf("%d ", x(i))
end
fprintf("\n");

xsol = A\b;

fprintf("Exact solution:\n");
for i=1:N
  fprintf("%d ", xsol(i))
end
fprintf("\n");
