/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
/*!
  \file   ConjugateGradient.cpp
  \brief  

  The Conjugate Gradient class is a concrete vertex mover
  which performs conjugate gradient minimizaiton.

  \author Michael Brewer
  \date   2002-06-19
*/

#include "ConjugateGradient.hpp"
#include <math.h>
#include "MsqDebug.hpp"
#include "MsqTimer.hpp"
//#include "MsqFreeVertexIndexIterator.hpp"

namespace MESQUITE_NS {

extern int get_parallel_rank();
extern int get_parallel_size();

std::string ConjugateGradient::get_name() const 
  { return "ConjugateGradient"; }
  
PatchSet* ConjugateGradient::get_patch_set()
  { return PatchSetUser::get_patch_set(); }

ConjugateGradient::ConjugateGradient(ObjectiveFunction* objective) :
  VertexMover(objective),
  PatchSetUser(true),
  pMemento(NULL),
  conjGradDebug(0)
{ }  

ConjugateGradient::ConjugateGradient(ObjectiveFunction* objective,
                                     MsqError &err) :
  VertexMover(objective),
  PatchSetUser(true),
  pMemento(NULL),
  conjGradDebug(0)
{
    //Michael:: default to global?
  set_debugging_level(0);
    //set the default inner termination criterion
  TerminationCriterion* default_crit=get_inner_termination_criterion();
  if(default_crit==NULL){
    MSQ_SETERR(err)("QualityImprover did not create a default inner "
                    "termination criterion.", MsqError::INVALID_STATE);
    return;
  }
  else{
    default_crit->add_iteration_limit( 5 );
    MSQ_ERRRTN(err);
  }
  
}  


ConjugateGradient::~ConjugateGradient()
{
  // Checks that cleanup() has been called. 
  assert(pMemento==NULL);
}
  
void ConjugateGradient::initialize(PatchData &pd, MsqError &err)
{
  if (get_parallel_size())
    MSQ_DBGOUT(2) << "\nP[" << get_parallel_rank() << "] " << "o   Performing Conjugate Gradient optimization.\n";
  else
    MSQ_DBGOUT(2) << "\no   Performing Conjugate Gradient optimization.\n";
  pMemento=pd.create_vertices_memento(err);
}


void ConjugateGradient::initialize_mesh_iteration(PatchData &/*pd*/,
                                                  MsqError &/*err*/)
{
  
}

/*!Performs Conjugate gradient minimization on the PatchData, pd.*/
void ConjugateGradient::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err){
  // pd.reorder();

  MSQ_FUNCTION_TIMER( "ConjugateGradient::optimize_vertex_positions" );

  Timer c_timer;
  size_t num_vert=pd.num_free_vertices();
  if(num_vert<1){
     MSQ_DBGOUT(1) << "\nEmpty free vertex list in ConjugateGradient\n";
     return;
  }
/*  
    //zero out arrays
  int zero_loop=0;
  while(zero_loop<arraySize){
    fGrad[zero_loop].set(0,0,0);
    pGrad[zero_loop].set(0,0,0);
    fNewGrad[zero_loop].set(0,0,0);
    ++zero_loop;
  }
*/
  
  // get OF evaluator
  OFEvaluator& objFunc = get_objective_function_evaluator();
  
  size_t ind;
    //Michael cull list:  possibly set soft_fixed flags here
  
   //MsqFreeVertexIndexIterator free_iter(pd, err);  MSQ_ERRRTN(err);
  
      
   double f=0;
   //Michael, this isn't equivalent to CUBIT because we only want to check
   //the objective function value of the 'bad' elements
   //if invalid initial patch set an error.
   bool temp_bool = objFunc.update(pd, f, fGrad, err);
   assert(fGrad.size() == num_vert);
   if(MSQ_CHKERR(err))
      return;
  if( ! temp_bool){
    MSQ_SETERR(err)("Conjugate Gradient not able to get valid gradient "
                    "and function values on intial patch.", 
                    MsqError::INVALID_MESH);
    return;
  }
  double grad_norm=MSQ_MAX_CAP;
  
  if(conjGradDebug>0){
    MSQ_PRINT(2)("\nCG's DEGUB LEVEL = %i \n",conjGradDebug);
    grad_norm=Linf(arrptr(fGrad),fGrad.size());
    MSQ_PRINT(2)("\nCG's FIRST VALUE = %f,grad_norm = %f",f,grad_norm);
    MSQ_PRINT(2)("\n   TIME %f",c_timer.since_birth());
    grad_norm=MSQ_MAX_CAP;
  }

    //Initializing pGrad (search direction).  
  pGrad.resize(fGrad.size());
  for (ind = 0; ind < num_vert; ++ind)
    pGrad[ind]=(-fGrad[ind]);

  int j=0; // total nb of step size changes ... not used much
  int i=0; // iteration counter
  unsigned m=0; // 
  double alp=MSQ_MAX_CAP; // alp: scale factor of search direction
    //we know inner_criterion is false because it was checked in
    //loop_over_mesh before being sent here.
  TerminationCriterion* term_crit=get_inner_termination_criterion();
  
    //while ((i<maxIteration && alp>stepBound && grad_norm>normGradientBound)
    //     && !inner_criterion){
  while(!term_crit->terminate()){
    ++i;
      //std::cout<<"\Michael delete i = "<<i;
    int k=0;
    alp=get_step(pd,f,k,err);
    j+=k;
    if(conjGradDebug>2){
      MSQ_PRINT(2)("\n  Alp initial, alp = %20.18f",alp);
    }

     // if alp == 0, revert to steepest descent search direction
    if(alp==0){
      for (m = 0; m < num_vert; ++m) {
        pGrad[m]=(-fGrad[m]);
      }
      alp=get_step(pd,f,k,err);
      j+=k;
      if(conjGradDebug>1){
        MSQ_PRINT(2)("\n CG's search direction reset.");
        if(conjGradDebug>2)
          MSQ_PRINT(2)("\n  Alp was zero, alp = %20.18f",alp);
      }
      
    }
    if(alp!=0){
      pd.move_free_vertices_constrained( arrptr(pGrad), num_vert, alp, err );
      MSQ_ERRRTN(err);
      
      if (! objFunc.update(pd, f, fNewGrad, err)){
        MSQ_SETERR(err)("Error inside Conjugate Gradient, vertices moved "
                        "making function value invalid.", 
                        MsqError::INVALID_MESH);
        return;
      }
      assert(fNewGrad.size() == (unsigned)num_vert);
      
      if(conjGradDebug>0){
        grad_norm=Linf(arrptr(fNewGrad),num_vert);
        MSQ_PRINT(2)("\nCG's VALUE = %f,  iter. = %i,  grad_norm = %f,  alp = %f",f,i,grad_norm,alp);
        MSQ_PRINT(2)("\n   TIME %f",c_timer.since_birth());
      }
      double s11=0;
      double s12=0;
      double s22=0;
      //free_iter.reset();
      //while (free_iter.next()) {
      //  m=free_iter.value();
      for (m = 0; m < num_vert; ++m) {
        s11+=fGrad[m]%fGrad[m];
        s12+=fGrad[m]%fNewGrad[m];
        s22+=fNewGrad[m]%fNewGrad[m];
      }
      
        // Steepest Descent (takes 2-3 times as long as P-R)
        //double bet=0;
      
        // Fletcher-Reeves (takes twice as long as P-R)
        //double bet = s22/s11;

        // Polack-Ribiere  
      double bet;  
      if (!divide( s22-s12, s11, bet ))
        return; // gradient is zero    
      //free_iter.reset();
      //while (free_iter.next()) {
      //  m=free_iter.value();
      for (m = 0; m < num_vert; ++m) {
        pGrad[m]=(-fNewGrad[m]+(bet*pGrad[m]));
        fGrad[m]=fNewGrad[m];
      }
      if(conjGradDebug>2){
        MSQ_PRINT(2)(" \nSEARCH DIRECTION INFINITY NORM = %e",
                   Linf(arrptr(fNewGrad),num_vert));
      }
      
    }//end if on alp == 0

    term_crit->accumulate_patch( pd, err ); MSQ_ERRRTN(err);
    term_crit->accumulate_inner( pd, f, arrptr(fGrad), err );  MSQ_ERRRTN(err);
  }//end while
  if(conjGradDebug>0){
    MSQ_PRINT(2)("\nConjugate Gradient complete i=%i ",i);
    MSQ_PRINT(2)("\n-  FINAL value = %f, alp=%4.2e grad_norm=%4.2e",f,alp,grad_norm);
    MSQ_PRINT(2)("\n   FINAL TIME %f",c_timer.since_birth());
  }
}


void ConjugateGradient::terminate_mesh_iteration(PatchData &/*pd*/,
                                                 MsqError &/*err*/)
{
    //  cout << "- Executing ConjugateGradient::iteration_complete()\n";
}


void ConjugateGradient::cleanup()
{
    //  cout << "- Executing ConjugateGradient::iteration_end()\n";
  fGrad.clear();
  pGrad.clear();
  fNewGrad.clear();
    //pMemento->~PatchDataVerticesMemento();
  delete pMemento;
  pMemento = NULL;
}

//!Computes a distance to move vertices given an initial position and search direction (stored in data member pGrad).
/*!Returns alp, the double which scales the search direction vector
  which when added to the old nodal positions yields the new nodal
  positions.*/
/*!\todo Michael NOTE:  ConjugateGradient::get_step's int &j is only
  to remain consisitent with CUBIT for an initial test.  It can be
  removed.*/

double ConjugateGradient::get_step(PatchData &pd,double f0,int &j,
                                   MsqError &err)
{
  // get OF evaluator
  OFEvaluator& objFunc = get_objective_function_evaluator();

  size_t num_vertices=pd.num_free_vertices();
    //initial guess for alp
  double alp=1.0;
  int jmax=100;
  double rho=0.5;
    //feasible=false implies the mesh is not in the feasible region
  bool feasible=false;
  int found=0;
    //f and fnew hold the objective function value
  double f=0;
  double fnew=0;
    //Counter to avoid infinitly scaling alp
  j=0;
  //save memento
  pd.recreate_vertices_memento(pMemento, err);
    //if we must check feasiblility
    //while step takes mesh into infeasible region and ...
  while (j<jmax && !feasible && alp>MSQ_MIN) {
    ++j;
    pd.set_free_vertices_constrained(pMemento,arrptr(pGrad),num_vertices,alp,err);
    feasible=objFunc.evaluate(pd,f,err);
    if(err.error_code() == err.BARRIER_VIOLATED)
      err.clear();  // barrier violation does not represent an actual error here
    MSQ_ERRZERO(err);
      //if not feasible, try a smaller alp (take smaller step)
    if(!feasible){
      alp*=rho;
    }
  }//end while ...
  
    //if above while ended due to j>=jmax, no valid step was found.
  if(j>=jmax){
    MSQ_PRINT(2)("\nFeasible Point Not Found");
    return 0.0;
  }
    //Message::print_info("\nOriginal f %f, first new f = %f, alp = %f",f0,f,alp);
    //if new f is larger than original, our step was too large
  if(f>=f0){
    j=0;
    while (j<jmax && found == 0){
      ++j;
      alp *= rho;
      pd.set_free_vertices_constrained(pMemento,arrptr(pGrad),num_vertices,alp,err);
        //Get new obj value
        //if patch is now invalid, then the feasible region is  convex or
        //we have an error.  For now, we assume an error.
      if(! objFunc.evaluate(pd,f,err) ){
        MSQ_SETERR(err)("Non-convex feasiblility region found.",MsqError::INVALID_MESH);
      }
      pd.set_to_vertices_memento(pMemento,err);MSQ_ERRZERO(err);
        //if our step has now improved the objective function value
      if(f<f0){
        found=1;
      }
    }//   end while j less than jmax
      //Message::print_info("\nj = %d found = %d f = %20.18f f0 = %20.18f\n",j,found,f,f0);
      //if above ended because of j>=jmax, take no step
    if(found==0){
        //Message::print_info("alp = %10.8f, but returning zero\n",alp);
      alp=0.0; 
      return alp;
    }

    j=0;
      //while shrinking the step improves the objFunc value further,
      //scale alp down.  Return alp, when scaling once more would
      //no longer improve the objFunc value.  
    while(j<jmax){
      ++j;
      alp*=rho;
      //step alp in search direction from original positions
      pd.set_free_vertices_constrained(pMemento,arrptr(pGrad),num_vertices,alp,err);MSQ_ERRZERO(err);

        //get new objective function value
      if (! objFunc.evaluate(pd,fnew,err))
        MSQ_SETERR(err)("Non-convex feasiblility region found while "
                        "computing new f.",MsqError::INVALID_MESH);
      if(fnew<f){
        f=fnew;
      }
      else{
	//Reset the vertices to original position
	pd.set_to_vertices_memento(pMemento,err);MSQ_ERRZERO(err);
	alp/=rho;
	return alp;
      }
    }
    //Reset the vertices to original position and return alp
    pd.set_to_vertices_memento(pMemento,err);MSQ_ERRZERO(err);
    return alp;
  }
    //else our new f was already smaller than our original
  else{
    j=0;
      //check to see how large of step we can take
    while (j<jmax && found == 0) {
      ++j;
        //scale alp up (rho must be less than 1)
      alp /= rho;
      //step alp in search direction from original positions
      pd.set_free_vertices_constrained(pMemento,arrptr(pGrad),num_vertices,alp,err);MSQ_ERRZERO(err);

      feasible = objFunc.evaluate(pd,fnew, err);
      if(err.error_code() == err.BARRIER_VIOLATED)
        err.clear();  // evaluate() error does not represent an actual problem here
      MSQ_ERRZERO(err);
      if ( ! feasible ){
        alp *= rho;
          //Reset the vertices to original position and return alp
        pd.set_to_vertices_memento(pMemento,err);MSQ_ERRZERO(err);
        return alp;
      }
      if (fnew<f) { 
          f = fnew; 
       }
       else {
         found=1;
         alp *= rho;
       }
    }

    //Reset the vertices to original position and return alp
    pd.set_to_vertices_memento(pMemento,err);MSQ_ERRZERO(err);
    return alp;
  }
}

/*!Quadratic one-dimensional line search.*/
/*
double ConjugateGradient::get_step(PatchData &pd,double f0,int &j,
                                   MsqError &err){
  const double CGOLD = 0.3819660;
  const double ZEPS = 1.0e-10;
  int n=pd.num_free_vertices();
  MsqVertex* vertices=pd.get_vertex_array(err);
  double a,b,d,etemp,fb,fu,fv,fw,fx,p,q,r,tol,tol1,tol2,u,v,w,x,xm;
  double e=0.0;
  d=0.0;
  tol=.001;
  int iter, maxiter;
  maxiter=100;
  a=0;
  b=.125;
  int m=0;
  fb=f0-1.0;
  iter=0;
  //find b such that a b 'should' bracket the min
  while (fb<=f0 && iter<maxiter){
    ++iter;
    b*=2.0;
    for(m=0;m<n;++m){
      mCoord[m]=mCoord[m] + (b*pGrad[m]);
      vertices[m]=(mCoord[m]);
    }
    fb=objFunc->evaluate(pd,err);
  }
  iter=0;
  x=w=v=(b/2.0);
  for(m=0;m<n;++m){
    mCoord[m]=mCoord[m] + (x*pGrad[m]);
    vertices[m]=(mCoord[m]);
  }
  fw=fv=fx=objFunc->evaluate(pd,err);
  for(iter=0;iter<maxiter;++iter){
      //Message::print_info("a=%f,b=%f,x=%f,iter=%i\n",a,b,x,iter);
    xm=(a+b)*.5;
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if(fabs(x-xm)<= (tol2-0.5*(b-a))){
      return x;
    }
    if(fabs(e)>tol1){
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if(q>0.0)
        p=-p;
      q=fabs(q);
      etemp=e;
      e=d;
      if(fabs(p)>=fabs(0.5*q*etemp)||(p<=q*(a-x))||(p>=q*(b-x))){
        d=CGOLD*(e=(x>=xm?a-x:b-x));
      }
      else{
        d=p/q;
        u=x+d;
        if(u-a<tol2||b-u<tol2)
        {
          if(tol1<0.0)
            d=x-xm;
          else
            d=xm-x;
        }
      }
    }
    
    else{
      d=CGOLD*(e=(x>=xm?a-x:b-x));
    }
    if(tol<0.0)
      u=(fabs(d)>=tol1?x+d:x-d);
    else
      u=(fabs(d)>=tol1?x+d:x+d);
    for(m=0;m<n;++m){
      mCoord[m]=mCoord[m] + (u*pGrad[m]);
      vertices[m]=(mCoord[m]);
    }
    fu=objFunc->evaluate(pd,err);
    if(fu<fx){
      if(u>=x)
        a=x;
      else
        b=x;
      v=w;
      w=x;
      x=u;
      fv=fw;
      fw=fx;
      fx=fu;
    }
    else{
      if(u<x)
        a=u;
      else
        b=u;
      if(fu<=fw||w==x){
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      }
      else if (fu<=fv||v==x||v==w){
        v=u;
        fv=fu;
      }
    }
  }
  for(m=0;m<n;++m){   
    vertices[m]=(mCoord[m]);
  }
    //PRINT_WARNING("TOO MANY ITERATIONS IN QUADRATIC LINE SEARCH");
  return x;
}
*/

} //namespace Mesquite
    



