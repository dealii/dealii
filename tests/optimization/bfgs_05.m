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

% a companion to bfgs_05.cc which minimizes the Rosenbrok function using Octave.

% This file uses bfgsmin() from optim package: https://octave.sourceforge.io/optim/
% https://octave.sourceforge.io/optim/package_doc/bfgsmin.html#bfgsmin
%
% take Example 3 from bfgsmin_example.m
pkg load optim

function [obj_value, gradient] = objective(theta, location)
	x = theta - location + ones(rows(theta),1); % move minimizer to "location"
	[obj_value, gradient] = rosenbrock(x);
endfunction

% problem parameters
dim = 20; % dimension of Rosenbrock function
theta0 = zeros(dim+1,1);  % starting values
location = (0:dim)/dim;   % true values
location = location';

% solver parameers
mMax    = 3;      % maximum number of stored residuals
itmax   = 120;    % maximum allowable number of iterations
ftol    = 1;      % function change tolerance
xtol    = 1;      % parameter change tolerance
gtol    = 1e-5;   % gradient tolerance

verb    = 2;      % verbosity [0,3]

control = {itmax;verb;1;1;mMax;ftol;xtol;gtol};

[theta, obj_value, convergence] = bfgsmin("objective", {theta0, location}, control);

linf_norm = norm(theta-location, 'inf')
