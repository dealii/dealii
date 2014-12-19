// ---------------------------------------------------------------------
//
// Copyright (C) 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include <deal.II/base/time_stepping.h>
#include <deal.II/lac/vector.h>

// test Runge-Kutta methods
Vector<double> f1(double const t, Vector<double> const &y) 
{ 
  Vector<double> values(y);
  for (unsigned int i=0; i<values.size(); ++i)
    values[i] = 1.0;

  return values;
}

Vector<double> f2(double const t, Vector<double> const &y) 
{ 
  Vector<double> values(y);
  for (unsigned int i=0; i<values.size(); ++i)
    values[i] = 2.0*t;

  return values; 
}

Vector<double> f3(double const t, Vector<double> const &y) 
{ 
  Vector<double> values(y);
  for (unsigned int i=0; i<values.size(); ++i)
    values[i] = 3.0*t*t;

  return values; 
}

Vector<double> f4(double const t, Vector<double> const &y) 
{ 
  Vector<double> values(y);
  for (unsigned int i=0; i<values.size(); ++i)
    values[i] = 4.0*t*t*t;
  
  return values; 
}

Vector<double> f5(double const t, Vector<double> const &y) 
{
  Vector<double> values(y);
  for (unsigned int i=0; i<values.size(); ++i)
    values[i] = 5.0*t*t*t*t;
  
  return values; 
}

Vector<double> id_minus_tau_J_inv1(double const t, double const tau, Vector<double> const &y) 
{
  return y;
}

Vector<double> id_minus_tau_J_inv2(double const t, double const tau, Vector<double> const &y) 
{
  return y;
}

Vector<double> id_minus_tau_J_inv3(double const t, double const tau, Vector<double> const &y) 
{
  return y;
}

Vector<double> id_minus_tau_J_inv4(double const t, double const tau, Vector<double> const &y) 
{
  return y;
}

Vector<double> id_minus_tau_J_inv5(double const t, double const tau, Vector<double> const &y) 
{
  return y;
}

double my1(double const t) 
{ 
  return t; 
}

double my2(double const t) 
{ 
  return t*t; 
}

double my3(double const t) 
{ 
  return t*t*t; 
}

double my4(double const t) 
{ 
  return t*t*t*t; 
}

double my5(double const t) 
{ 
  return t*t*t*t*t; 
}

void test(TimeStepping::RungeKutta<Vector<double> > &solver,
    std_cxx11::function<Vector<double> (double const, Vector<double> const &)> f,
    std_cxx11::function<Vector<double> (double const, double const, Vector<double> const &)> id_minus_tau_J_inv,
    std_cxx11::function<double (double const)> my)
{
  unsigned int n_time_steps = 1;
  unsigned int size = 1;
  double initial_time = 0.0, final_time = 1.0;
  double time_step = (final_time-initial_time)/static_cast<double> (n_time_steps);
  double time = initial_time;
  Vector<double> solution(size);
  Vector<double> exact_solution(size);
  for (unsigned int i=0; i<size; ++i)
  {
    solution[i] = my(initial_time);
    exact_solution[i] = my(final_time);
  }

  for (unsigned int i=0; i<n_time_steps; ++i)
    time = solver.evolve_one_time_step(f,id_minus_tau_J_inv,time,time_step,solution);

  Vector<double> error(exact_solution);
  error.sadd(1.0,-1.0,solution);
  double error_norm = error.l2_norm();
  deallog << error_norm <<std::endl;
}

void test2(TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > &solver,
    std_cxx11::function<Vector<double> (double const, Vector<double> const &)> f,
    std_cxx11::function<Vector<double> (double const, double const, Vector<double> const &)> id_minus_tau_J_inv,
    std_cxx11::function<double (double const)> my)
{
  double initial_time = 0.0, final_time = 1.0;
  double time_step = 1.0;
  unsigned int size = 1;
  unsigned int n_time_steps = 0;
  double time = initial_time;
  Vector<double> solution(size);
  Vector<double> exact_solution(size);
  for (unsigned int i=0; i<size; ++i)
  {
    solution[i] = my(initial_time);
    exact_solution[i] = my(final_time);
  }


  while (time<final_time)
  {
    if (time+time_step>final_time)
      time_step = final_time - time;
    time = solver.evolve_one_time_step(f,id_minus_tau_J_inv,time,time_step,solution);
    time_step = solver.get_status().delta_t_guess;
  }

  Vector<double> error(exact_solution);
  error.sadd(1.0,-1.0,solution);
  double error_norm = error.l2_norm();
  deallog<<error_norm<<std::endl;
}

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog<<"Forward Euler"<<std::endl;
  TimeStepping::ExplicitRungeKutta<Vector<double> > fe(TimeStepping::FORWARD_EULER);
  test(fe,f1,id_minus_tau_J_inv1,my1);

  deallog<<"Runge-Kutta third order"<<std::endl;
  TimeStepping::ExplicitRungeKutta<Vector<double> > rk3(TimeStepping::RK_THIRD_ORDER);
  test(rk3,f3,id_minus_tau_J_inv3,my3);

  deallog<<"Runge-Kutta fourth order"<<std::endl;
  TimeStepping::ExplicitRungeKutta<Vector<double> > rk4(TimeStepping::RK_CLASSIC_FOURTH_ORDER);
  test(rk4,f4,id_minus_tau_J_inv4,my4);

  deallog<<"Backward Euler"<<std::endl;
  TimeStepping::ImplicitRungeKutta<Vector<double> > be(TimeStepping::BACKWARD_EULER);
  test(be,f1,id_minus_tau_J_inv1,my1);

  deallog<<"Implicit midpoint"<<std::endl;
  TimeStepping::ImplicitRungeKutta<Vector<double> > im(TimeStepping::IMPLICIT_MIDPOINT);
  test(im,f2,id_minus_tau_J_inv2,my2);

  deallog<<"Crank-Nicolson"<<std::endl;
  TimeStepping::ImplicitRungeKutta<Vector<double> > cn(TimeStepping::CRANK_NICOLSON);
  test(cn,f2,id_minus_tau_J_inv2,my2);

  deallog<<"SDIRK"<<std::endl;
  TimeStepping::ImplicitRungeKutta<Vector<double> > sdirk(TimeStepping::SDIRK_TWO_STAGES);
  test(sdirk,f2,id_minus_tau_J_inv2,my2);

  deallog<<"Heun-Euler"<<std::endl;
  TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > he(TimeStepping::HEUN_EULER);
  test2(he,f2,id_minus_tau_J_inv2,my2);

  deallog<<"Bogacki-Shampine"<<std::endl;
  TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > bs(TimeStepping::BOGACKI_SHAMPINE);
  test2(bs,f3,id_minus_tau_J_inv3,my3);
  bs.free_memory();

  deallog<<"DOPRI"<<std::endl;
  TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > dopri(TimeStepping::DOPRI);
  test2(dopri,f5,id_minus_tau_J_inv5,my5);
  dopri.free_memory();

  deallog<<"Fehlberg"<<std::endl;
  TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > fehlberg(TimeStepping::FEHLBERG);
  test2(fehlberg,f5,id_minus_tau_J_inv5,my5);

  deallog<<"Cash-Karp"<<std::endl;
  TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > ck(TimeStepping::CASH_KARP);
  test2(ck,f5,id_minus_tau_J_inv5,my5);

  return 0;
}
