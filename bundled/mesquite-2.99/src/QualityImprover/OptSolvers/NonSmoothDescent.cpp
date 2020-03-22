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
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*!  \File NonSmoothDescent.cpp \brief
  
  Implements the NonSmoothDescent class member functions.
  
  \author Lori Freitag
  \date 2002-07-20 */

#include <stdlib.h>
#include <stdio.h>
#include "NonSmoothDescent.hpp"
#include "MsqTimer.hpp"

#undef NUMERICAL_GRADIENT

namespace MESQUITE_NS {

const double EPSILON = 1e-16;
const int MSQ_MAX_OPT_ITER = 20;

enum Rotate {
  COUNTERCLOCKWISE = 1,
  CLOCKWISE = 0
};


NonSmoothDescent::NonSmoothDescent(QualityMetric* qm)
  : currentQM(qm)
{

  MSQ_DBGOUT(1) << "- Executed NonSmoothDescent::NonSmoothDescent()\n";
}  

std::string NonSmoothDescent::get_name() const
  { return "NonSmoothDescent"; }

PatchSet* NonSmoothDescent::get_patch_set()
  { return &patchSet; }
  
void NonSmoothDescent::initialize(PatchData &/*pd*/, MsqError &err)
{
  minStepSize = 1e-6;
  MSQ_DBGOUT(1) << "- Executed NonSmoothDescent::initialize()\n";
}

void NonSmoothDescent::initialize_mesh_iteration(PatchData &/*pd*/,
                                                         MsqError &/*err*/)
{
}
void NonSmoothDescent::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  MSQ_FUNCTION_TIMER( "NonSmoothDescent" );

  //  cout << "- Executing NonSmoothDescent::optimize_node_positions()\n";
  /* perform the min max smoothing algorithm */
  MSQ_PRINT(2)("\nInitializing the patch iteration\n");

  MSQ_PRINT(3)("Number of Vertices: %d\n",(int)pd.num_nodes());
  MSQ_PRINT(3)("Number of Elements: %d\n",(int)pd.num_elements());
    //Michael: Note: is this a reliable way to get the dimension?
  // NOTE: Mesquite only does 3-dimensional (including surface) meshes.
  // mDimension = pd.get_mesh()->get_geometric_dimension(err); MSQ_ERRRTN(err);
  // MSQ_PRINT(3)("Spatial Dimension: %d\n",mDimension);
  MSQ_PRINT(3)("Spatial Dimension: 3\n");

  MSQ_PRINT(3)("Num Free = %d\n",(int)pd.num_free_vertices());

  MsqFreeVertexIndexIterator free_iter(pd, err); MSQ_ERRRTN(err);
  free_iter.reset();
  free_iter.next(); 
  freeVertexIndex = free_iter.value();
  MSQ_PRINT(3)("Free Vertex Index = %lu\n",(unsigned long)freeVertexIndex);

  // TODO - need to switch to validity via metric evaluations should
  // be associated with the compute_function somehow
  /* check for an invalid mesh; if it's invalid return and ask the user 
     to use untangle */
  if (!this->validity_check(pd,err)) {
      MSQ_PRINT(1)("ERROR: Invalid mesh\n");
      MSQ_SETERR(err)("Invalid Mesh: Use untangle to create a valid "
                      "triangulation", MsqError::INVALID_MESH);
      return;
  }

  /* initialize the optimization data up to numFunctionValues */
  this->init_opt(pd,err); MSQ_ERRRTN(err);
  this->init_max_step_length(pd,err);  MSQ_ERRRTN(err);
  MSQ_PRINT(3)("Done initializing optimization\n");

  /* compute the initial function values */
  //TODO this should return a bool with the validity
  this->compute_function(&pd, originalFunction, err);  MSQ_ERRRTN(err);
 
  // find the initial active set
  this->find_active_set(originalFunction, mActive);

  this->minmax_opt(pd,err);  MSQ_ERRRTN(err);
}


void NonSmoothDescent::terminate_mesh_iteration( PatchData &/*pd*/,
                                                 MsqError &/*err*/)
{
}
  
void NonSmoothDescent::cleanup()
{
  MSQ_DBGOUT(1) << "- Executing NonSmoothDescent::cleanup()\n";
  MSQ_DBGOUT(1) << "- Done with NonSmoothDescent::cleanup()\n";
}





void NonSmoothDescent::find_plane_points( Direction dir1, 
                                          Direction dir2,
                                          const std::vector<Vector3D>& vec,
                                          Vector3D& pt1,
                                          Vector3D& pt2,
                                          Vector3D& /*pt3*/,
                                          Status& status,
                                          MsqError& )
{
    int i;
    int num_min, num_max;
    Rotate rotate=CLOCKWISE;
    int num_rotated=0;
    double pt_1, pt_2;
    double min, inv_slope;
    double min_inv_slope=0.;
    double max; 
    double max_inv_slope=0;
    double inv_origin_slope=0;
    const int num_vec = vec.size();
    const int auto_ind_size = 50;
    int auto_ind[auto_ind_size];
    std::vector<int> heap_ind;
    int* ind;
    if (num_vec <= auto_ind_size)
      ind = auto_ind;
    else {
      heap_ind.resize(num_vec);
      ind = &heap_ind[0];
    }

    status = MSQ_CHECK_BOTTOM_UP;
    /* find the minimum points in dir1 starting at -1 */
    num_min = 0; ind[0]=-1; ind[1]=-1; ind[2]=-1; min=1.0;
    for (i=0;i<num_vec;i++) {
      if (vec[i][dir1]<min) {
	min = vec[i][dir1]; ind[0] = i; num_min = 1;
      } else if (fabs(vec[i][dir1] - min) < EPSILON) {
	ind[num_min++] = i;
      }
    }
    if (min >= 0) status = MSQ_NO_EQUIL;
 
    if (status != MSQ_NO_EQUIL) {
      switch(num_min) {
      case 1: /* rotate to find the next point */
        pt1 = vec[ind[0]];
	pt_1 = pt1[dir1]; pt_2 = pt1[dir2];
	if (pt1[dir2] <= 0){rotate=COUNTERCLOCKWISE; max_inv_slope=-HUGE_VAL;}
	if (pt1[dir2] > 0){rotate=CLOCKWISE; min_inv_slope=HUGE_VAL;}
	switch(rotate) {
	case COUNTERCLOCKWISE:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if ((inv_slope>max_inv_slope) &&  
		  (fabs(inv_slope - max_inv_slope) > EPSILON)) {
		ind[1] = i; max_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - max_inv_slope) < EPSILON) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	  break;
	case CLOCKWISE:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if ((inv_slope<min_inv_slope) && 
		  (fabs(inv_slope - max_inv_slope) > EPSILON)){
		ind[1] = i; min_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - min_inv_slope) < EPSILON) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	}
	switch(num_rotated) {
	case 0:
	  MSQ_PRINT(3)("No points in the rotation ... odd\n");
	    status = MSQ_HULL_TEST_ERROR;
	  break;
	case 1:
	  MSQ_PRINT(3)("Found a line in the convex hull\n");
          pt2 = vec[ind[1]];
	  status = MSQ_TWO_PT_PLANE;
	  break;
	default:
	  MSQ_PRINT(3)("Found 2 or more points in the rotation\n");
	    if (fabs(pt_1) > EPSILON) inv_origin_slope = pt_2/pt_1;
	  switch(rotate) {
	  case COUNTERCLOCKWISE:
	    if (inv_origin_slope >= max_inv_slope) status=MSQ_NO_EQUIL;
	    else status=MSQ_CHECK_TOP_DOWN;
	    break;
	  case CLOCKWISE:
	    if (inv_origin_slope <= min_inv_slope) status=MSQ_NO_EQUIL;
	    else status=MSQ_CHECK_TOP_DOWN;
	  }
	}
	break;
      case 2: /* use these two points to define the plane */
	MSQ_PRINT(3)("Found two minimum points to define the plane\n");
        pt1 = vec[ind[0]];
        pt2 = vec[ind[1]];
	status = MSQ_TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	MSQ_PRINT(3)("Found 3 or more points in min plane %f\n",min);
	  if (vec[ind[0]][dir1] >= 0) status = MSQ_NO_EQUIL;
	  else status = MSQ_CHECK_TOP_DOWN;
    }
    }

    /***************************/
    /*  failed to find any information, checking top/down this coord*/
    /***************************/

    if (status == MSQ_CHECK_TOP_DOWN) {
    /* find the maximum points in dir1 starting at 1 */
    num_max = 0; ind[0]=-1; ind[1]=-1; ind[2]=-1; max=-1.0;
    for (i=0;i<num_vec;i++) {
      if (vec[i][dir1] > max) {
	max = vec[i][dir1]; ind[0] = i; num_max = 1;
      } else if (fabs(vec[i][dir1] - max) < EPSILON) {
	ind[num_max++] = i;
      }
    }
    if (max <= 0) status = MSQ_NO_EQUIL;
 
    if (status != MSQ_NO_EQUIL) {
      switch(num_max) {
      case 1: /* rotate to find the next point */
        pt1 = vec[ind[0]];
	pt_1 = pt1[dir1];  pt_2 = pt1[dir2];
	if (pt1[dir2] < 0){rotate=CLOCKWISE; min_inv_slope=HUGE_VAL;}
	if (pt1[dir2] >= 0){rotate=COUNTERCLOCKWISE; max_inv_slope=-HUGE_VAL;}
	switch(rotate) {
	case COUNTERCLOCKWISE:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if (inv_slope>max_inv_slope) {
		ind[1] = i; max_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - max_inv_slope) < EPSILON) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	  break;
	case CLOCKWISE:
	  for (i=0;i<num_vec;i++) {
	    if (i!=ind[0]) {
	      inv_slope = (vec[i][dir2] - pt_2)/(vec[i][dir1]-pt_1);
	      if (inv_slope<min_inv_slope) {
		ind[1] = i; min_inv_slope=inv_slope; num_rotated = 1;
	      } else if (fabs(inv_slope - min_inv_slope) < EPSILON) {
		ind[2] = i; num_rotated++;
	      }
	    }
	  }
	}
	switch(num_rotated) {
	case 0:
	  MSQ_PRINT(3)("No points in the rotation ... odd\n");
	  status = MSQ_HULL_TEST_ERROR;
	  break;
	case 1:
	  MSQ_PRINT(3)("Found a line in the convex hull\n");
          pt2 = vec[ind[1]];
	  status = MSQ_TWO_PT_PLANE;
	  break;
	default:
	  MSQ_PRINT(3)("Found 2 or more points in the rotation\n");
	    /* check to see if rotation got past origin */
	  inv_origin_slope = pt_2/pt_1;
	  switch(rotate) {
	  case COUNTERCLOCKWISE:
	    if (inv_origin_slope >= max_inv_slope) status=MSQ_NO_EQUIL;
	    else if (dir1 == 2) status=MSQ_CHECK_Y_COORD_DIRECTION;
	    else if (dir1 == 1) status=MSQ_CHECK_X_COORD_DIRECTION;
	    else status=MSQ_EQUIL;
	    break;
	  case CLOCKWISE:
	    if (inv_origin_slope <= min_inv_slope) status=MSQ_NO_EQUIL;
	    else if (dir1 == 2) status=MSQ_CHECK_Y_COORD_DIRECTION;
	    else if (dir1 == 1) status=MSQ_CHECK_X_COORD_DIRECTION;
	    else status=MSQ_EQUIL;
	  }
	}
	break;
      case 2: /* use these two points to define the plane */
        pt1 = vec[ind[0]];
        pt2 = vec[ind[1]];
	status = MSQ_TWO_PT_PLANE;
	break;
      default: /* check to see if all > 0 */
	MSQ_PRINT(3)("Found 3 in max plane %f\n",max);
	if (vec[ind[0]][dir1] <= 0) status = MSQ_NO_EQUIL;
	else if (dir1==2) status=MSQ_CHECK_Y_COORD_DIRECTION;
	else if (dir1==1) status=MSQ_CHECK_X_COORD_DIRECTION;
	else status = MSQ_EQUIL;
      }
    }
  }

}

void NonSmoothDescent::search_direction( PatchData &/*pd*/,
                                         Vector3D& mSearch,
                                         MsqError &err)
{
   bool       viable;
   double     a, b, c, denom;
   std::vector<Vector3D> dir;
   double     R0, R1;
   SymmetricMatrix P;
   double     x[2];
   double     search_mag;

   const int num_active = mActive.active_ind.size();;

   //TODO This might be o.k. actually - i don't see any dependence
   // on the element geometry here... try it and see if it works.
   // if not, try taking all of the gradients in the active set
   // and let the search direction be the average of those.
//   MSQ_FUNCTION_TIMER( "Search Direction" );

   MSQ_PRINT(2)("\nIn Search Direction\n");
   this->print_active_set(mActive, mFunction);
   
   if (num_active==0) {
       MSQ_SETERR(err)("No active values in search",MsqError::INVALID_STATE);
       return;
    }

    switch(num_active) {
    case 1: 
        mSearch = mGradient[mActive.active_ind[0]];
        mSteepest = mActive.active_ind[0];
        break;
    case 2:
        /* if there are two active points, move in the direction of the
	   intersection of the planes.  This is the steepest descent
           direction found by analytically solving the QP */
        
        /* set up the active gradient directions */
        this->get_active_directions(mGradient,dir,err); MSQ_ERRRTN(err);

        /* form the grammian */
        this->form_grammian(dir,err); MSQ_ERRRTN(err);
        this->form_PD_grammian(err); MSQ_ERRRTN(err);

        denom = (mG(0,0) + mG(1,1) - 2*mG(0,1));
        viable = true;
        if (fabs(denom) > EPSILON) {
	  /* gradients are LI, move along their intersection */
           b = (mG(0,0) - mG(0,1))/denom;  
           a = 1 - b;
           if ((b < 0) || (b > 1)) viable=false;  /* 0 < b < 1 */
           if (viable) {
             mSearch = a*dir[0] + b*dir[1];
           } else {
             /* the gradients are dependent, move along one face */
             mSearch = dir[0];
           }
        } else {
	   /* the gradients are dependent, move along one face */
           mSearch = dir[0];
        }
        mSteepest = mActive.active_ind[0];

        break;
    default:
        /* as in case 2: solve the QP problem to find the steepest
           descent direction.  This can be done analytically - as
           is done in Gill, Murray and Wright 
             for 3 active points in 3 directions - test PD of G
             otherwise we know it's SP SD so search edges and faces */

        /* get the active gradient directions */
        this->get_active_directions(mGradient,dir,err);  MSQ_ERRRTN(err);

        /* form the entries of the grammian matrix */
        this->form_grammian(dir,err); MSQ_ERRRTN(err);
        this->form_PD_grammian(err);  MSQ_ERRRTN(err);

        if (num_active == 3) {
          if (mG.condition3x3() < 1e14) { // if not singular
            /* form the entries of P=Z^T G Z where Z = [-1...-1; I ] */
            this->form_reduced_matrix(P); 
              /* form  the RHS and solve the system for the coeffs */
            R0 = mG(0,0) - mG(1,0);  R1 = mG(0,0) - mG(2,0);
            bool ok = this->solve2x2(P(0,0),P(0,1),P(1,0),P(1,1),R0,R1,x,err); MSQ_ERRRTN(err);
            if (ok) {
              a = 1 - x[0] - x[1];  b = x[0];  c = x[1];
              mSearch = a*dir[0] + b*dir[1] + c*dir[2];
              mSteepest = mActive.active_ind[0];
            } 
            else { 
              this->search_edges_faces(&dir[0], mSearch, err); MSQ_ERRRTN(err);
            }
          } 
          else {
            this->search_edges_faces(&dir[0], mSearch, err);  MSQ_ERRRTN(err);
          }
        } 
        else {
          this->search_edges_faces(&dir[0], mSearch, err); MSQ_ERRRTN(err);
        }
    }

    /* if the search direction is essentially zero, equilibrium pt */
    search_mag = mSearch % mSearch;
    MSQ_PRINT(3)("  Search Magnitude %g \n",search_mag);

    if (fabs(search_mag)<1E-13) 
      optStatus = MSQ_ZERO_SEARCH;
    else 
      mSearch /= std::sqrt(search_mag);
    MSQ_PRINT(3)("  Search Direction %g %g  Steepest %lu\n",mSearch[0],mSearch[1],(unsigned long)mSteepest);
}

void NonSmoothDescent::minmax_opt(PatchData &pd, MsqError &err)
{
      bool equilibriumPt = false;
      Vector3D mSearch(0,0,0);

//      int valid;
      MSQ_FUNCTION_TIMER( "Minmax Opt" );
      MSQ_PRINT(2)("In minmax_opt\n");

      mFunction = originalFunction;
      originalValue = mActive.true_active_value;

      int iterCount = 0;
      int optIterCount = 0;

      MSQ_PRINT(3)("Done copying original function to function\n");

      this->find_active_set(mFunction, mActive);
      prevActiveValues.clear();
      prevActiveValues.push_back( mActive.true_active_value );

     /* check for equilibrium point */
     /* compute the gradient */
     this->compute_gradient(&pd, mGradient, err); MSQ_ERRRTN(err);
     
     if (mActive.active_ind.size() >= 2) {
	MSQ_PRINT(3)("Testing for an equilibrium point \n");
	equilibriumPt = this->check_equilibrium( optStatus, err); MSQ_ERRRTN(err);

	if (MSQ_DBG(2) && equilibriumPt ) 
	  MSQ_PRINT(2)("Optimization Exiting: An equilibrium point \n");
     }

    /* terminate if we have found an equilibrium point or if the step is
       too small to be worthwhile continuing */
    while ((optStatus != MSQ_EQUILIBRIUM) && 
	   (optStatus != MSQ_STEP_TOO_SMALL) &&
	   (optStatus != MSQ_IMP_TOO_SMALL) &&
	   (optStatus != MSQ_FLAT_NO_IMP) &&
           (optStatus != MSQ_ZERO_SEARCH) &&
	   (optStatus != MSQ_MAX_ITER_EXCEEDED)) {

	/* increase the iteration count by one */
        /* smooth_param->iter_count += 1; */
        iterCount += 1;
        optIterCount += 1;
        if (iterCount > MSQ_MAX_OPT_ITER) optStatus = MSQ_MAX_ITER_EXCEEDED;

	MSQ_PRINT(3)("\nITERATION %d \n",iterCount);
	    
	/* compute the gradient */
	this->compute_gradient(&pd, mGradient, err); MSQ_ERRRTN(err);
        
	MSQ_PRINT(3)("Computing the search direction \n");
	this->search_direction(pd, mSearch, err); MSQ_ERRRTN(err);

	/* if there are viable directions to search */
	if ((optStatus != MSQ_ZERO_SEARCH) &&
            (optStatus != MSQ_MAX_ITER_EXCEEDED)) {

	    MSQ_PRINT(3)("Computing the projections of the gradients \n");
	    this->get_gradient_projections(mSearch, err); MSQ_ERRRTN(err);

	    MSQ_PRINT(3)("Computing the initial step size \n");
	    this->compute_alpha(err); MSQ_ERRRTN(err);

	    MSQ_PRINT(3)("Testing whether to accept this step \n");
	    this->step_acceptance(pd, iterCount, mSearch, err); MSQ_ERRRTN(err);
            //MSQ_PRINT(3)("The new free vertex position is %f %f %f\n",
            //  mCoords[freeVertexIndex][0],mCoords[freeVertexIndex][1],mCoords[freeVertexIndex][2]);

	    if (MSQ_DBG(3)) {
     		/* Print the active set */
	     	this->print_active_set(mActive, mFunction);
	    }

	    /* check for equilibrium point */
	    if (mActive.active_ind.size() >= 2) {
		MSQ_PRINT(3)("Testing for an equilibrium point \n");
                equilibriumPt = this->check_equilibrium(optStatus, err); MSQ_ERRRTN(err);

		if (MSQ_DBG(2) && equilibriumPt) 
		    MSQ_PRINT(2)("Optimization Exiting: An equilibrium point \n");
	    }

	    /* record the values */
            prevActiveValues.push_back( mActive.true_active_value );

	} else {
	    /* decrease the iteration count by one */
	    /* smooth_param->iter_count -= 1; */
	    iterCount -= 1;
	    if (MSQ_DBG(2)) {
		MSQ_PRINT(2)("Optimization Exiting: No viable directions; equilibrium point \n");
		/* Print the old active set */
		this->print_active_set(mActive,mFunction);
	    }
	}
      }

      MSQ_PRINT(2)("Checking the validity of the mesh\n");
      if (!this->validity_check(pd,err)) MSQ_PRINT(2)("The final mesh is not valid\n");
      MSQ_ERRRTN(err);

      MSQ_PRINT(2)("Number of optimization iterations %d\n", iterCount);
 
      switch(optStatus) {
        default:
          MSQ_PRINT(2)("Optimization Termination OptStatus: Invalid early termination\n"); break;
	case MSQ_EQUILIBRIUM:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Equilibrium\n"); break;
	case MSQ_STEP_TOO_SMALL:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Step Too Small\n"); break;
	case MSQ_IMP_TOO_SMALL:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Improvement Too Small\n"); break;
	case MSQ_FLAT_NO_IMP:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Flat No Improvement\n"); break;
	case MSQ_ZERO_SEARCH:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Zero Search\n"); break;
	case MSQ_MAX_ITER_EXCEEDED:
	  MSQ_PRINT(2)("Optimization Termination OptStatus: Max Iter Exceeded\n"); break;
      }
}


void NonSmoothDescent::step_acceptance( PatchData &pd, 
                                        int iterCount, 
                                        const Vector3D& mSearch, 
                                        MsqError &err )
{
  const double minAcceptableImprovement = 1e-6;

//  int        ierr;
//  int        num_values;
  bool       step_done;
  bool       valid = true, accept_alpha;
  double     estimated_improvement;
  double     current_improvement = HUGE_VAL;
  double     previous_improvement = HUGE_VAL;
  double     current_percent_diff = HUGE_VAL;
  Vector3D   original_point;

//  MSQ_FUNCTION_TIMER( "Step Acceptance" );
//  num_values = qmHandles.size();

  optStatus = MSQ_NO_STATUS;

  if (mAlpha < minStepSize) {
      optStatus = MSQ_IMP_TOO_SMALL;
      step_done = true;
      MSQ_PRINT(3)("Alpha starts too small, no improvement\n");
  }

  const MsqVertex* coords = pd.get_vertex_array(err);

  /* save the original function and active set */
  original_point = coords[freeVertexIndex];
  originalFunction = mFunction;
  originalActive = mActive;

  step_done = false;
  for (int num_steps = 0; !step_done && num_steps < 100; ++num_steps) {
    accept_alpha = false;

    while (!accept_alpha && mAlpha>minStepSize) {
 
      /* make the step */
      pd.move_vertex( -mAlpha*mSearch, freeVertexIndex, err );
        //pd.set_coords_array_element(coords[freeVertexIndex],0,err);

      MSQ_PRINT(2)("search direction %f %f \n",mSearch[0],mSearch[1]); 
      MSQ_PRINT(2)("new vertex position %f %f \n",coords[freeVertexIndex][0],coords[freeVertexIndex][1]); 

      /* assume alpha is acceptable */
      accept_alpha=true;

      /* never take a step that makes a valid mesh invalid or worsens the quality */
      // TODO Validity check revision -- do the compute function up here
      // and then the rest based on validity
      valid = validity_check(pd,err); MSQ_ERRRTN(err);
      if (valid) {
        valid = improvement_check(err); MSQ_ERRRTN(err);
      }
      if (!valid) {
          accept_alpha=false;
          pd.move_vertex( mAlpha * mSearch, freeVertexIndex, err );
            //pd.set_coords_array_element(coords[freeVertexIndex],0,err);
          mAlpha = mAlpha/2;
          MSQ_PRINT(2)("Step not accepted, the new alpha %f\n",mAlpha); 

          if (mAlpha < minStepSize) {
 	        optStatus = MSQ_STEP_TOO_SMALL;
                step_done = true;
                MSQ_PRINT(2)("Step too small\n");
 	        /* get back the original point, mFunction, and active set */
                pd.set_vertex_coordinates( original_point, freeVertexIndex, err );
                mFunction = originalFunction;
                mActive = originalActive;
	  }
       }
    } 
         
    if (valid  && (mAlpha > minStepSize)) {
      /* compute the new function and active set */
      this->compute_function(&pd, mFunction, err); MSQ_ERRRTN(err);
      this->find_active_set(mFunction, mActive);
	
      /* estimate the minimum improvement by taking this step */
      this->get_min_estimate(&estimated_improvement, err); MSQ_ERRRTN(err);
      MSQ_PRINT(2)("The estimated improvement for this step: %f\n",
		   estimated_improvement); 
	
      /* calculate the actual increase */
      current_improvement = mActive.true_active_value - prevActiveValues[iterCount-1];

      MSQ_PRINT(3)("Actual improvement %f\n",current_improvement);

      /* calculate the percent difference from estimated increase */
      current_percent_diff = fabs(current_improvement-estimated_improvement)/
	fabs(estimated_improvement);

      /* determine whether to accept a step */
      if ((fabs(previous_improvement) > fabs(current_improvement)) && 
	  (previous_improvement < 0)) {
	/* accept the previous step - it was better */
	     MSQ_PRINT(2)("Accepting the previous step\n");
 
	/* subtract alpha in again (previous step) */
        pd.move_vertex( -mAlpha * mSearch, freeVertexIndex, err );
            //pd.set_coords_array_element(coords[freeVertexIndex],0,err);

	/* does this make an invalid mesh valid? */
   //TODO Validity check revisison
        valid = validity_check(pd,err); MSQ_ERRRTN(err);
        if (valid) valid=improvement_check(err); MSQ_ERRRTN(err);

	/* copy test function and active set */
        mFunction = testFunction;
        mActive = testActive;
 
	optStatus = MSQ_STEP_ACCEPTED;  step_done = true;
            
	/* check to see that we're still making good improvements */
	if (fabs(previous_improvement) < minAcceptableImprovement) {
	  optStatus = MSQ_IMP_TOO_SMALL; step_done = true;
	  MSQ_PRINT(2)("Optimization Exiting: Improvement too small\n");
	}

      } else if (((fabs(current_improvement) > fabs(estimated_improvement)) ||
		  (current_percent_diff < .1)) && (current_improvement<0)) {
	/* accept this step, exceeded estimated increase or was close */
	optStatus = MSQ_STEP_ACCEPTED;  step_done = true;

	/* check to see that we're still making good improvements */
	if (fabs(current_improvement) < minAcceptableImprovement) {
	  MSQ_PRINT(2)("Optimization Exiting: Improvement too small\n");
	  optStatus = MSQ_IMP_TOO_SMALL; step_done = true;
	}

      } else if ((current_improvement > 0) && (previous_improvement > 0) &&
		 (fabs(current_improvement) < minAcceptableImprovement) &&
		 (fabs(previous_improvement) < minAcceptableImprovement)) {

	/* we are making no progress, quit */
	optStatus = MSQ_FLAT_NO_IMP; step_done = true;
	MSQ_PRINT(2)("Opimization Exiting: Flat no improvement\n");
           
	/* get back the original point, function, and active set */
        pd.set_vertex_coordinates( original_point, freeVertexIndex, err ); MSQ_ERRRTN(err);
        mFunction = originalFunction;
        mActive = originalActive;
      }
      else
      {
	/* halve alpha and try again */
	/* add out the old step */
        pd.move_vertex( mAlpha * mSearch, freeVertexIndex, err );
            //pd.set_coords_array_element(coords[freeVertexIndex],0,err);

	/* halve step size */
	mAlpha = mAlpha/2; 
	MSQ_PRINT(3)("Step not accepted, the new alpha %f\n",mAlpha);

	if (mAlpha < minStepSize)
          {
	  /* get back the original point, function, and active set */
	  MSQ_PRINT(2)("Optimization Exiting: Step too small\n");
          pd.set_vertex_coordinates( original_point, freeVertexIndex, err ); MSQ_ERRRTN(err);
          mFunction = originalFunction;
          mActive = originalActive;
	  optStatus = MSQ_STEP_TOO_SMALL;  step_done = true;
	}
        else
        {
          testFunction = mFunction;
          testActive = mActive;
	  previous_improvement = current_improvement;
	}
      }
    }
  }
  if (current_improvement>0 && optStatus==MSQ_STEP_ACCEPTED) {
    MSQ_PRINT(2)("Accepted a negative step %f \n",current_improvement);
  }

}


bool NonSmoothDescent::compute_function(PatchData *patch_data, std::vector<double>& func, MsqError &err)
{
  // NEED 1.0/FUNCTION WHICH IS ONLY TRUE OF CONDITION NUMBER

  func.resize(qmHandles.size());
  bool valid_bool=true;
  for (size_t i = 0; i < qmHandles.size(); i++) {
    valid_bool = valid_bool &&
         currentQM->evaluate(*patch_data, qmHandles[i], func[i], err); MSQ_ERRZERO(err);
  }

  return valid_bool;
}


bool NonSmoothDescent::compute_gradient(PatchData *patch_data, 
                                        std::vector<Vector3D>& gradient_out,
                                        MsqError &err)
{
  MSQ_DBGOUT(2) << "Computing Gradient\n";

  bool valid = true;

#ifdef NUMERICAL_GRADIENT

  const double delta = 10e-6;
  std::vector<double> func(qmHandles.size()), fdelta(qmHandles.size());

  valid = this->compute_function(patch_data, func, err); MSQ_ERRZERO(err);

  /* gradient in the x, y, z direction */

  for (int j=0;j<3;j++) {

    // perturb the coordinates of the free vertex in the j direction by delta
    Vector3D delta_3( 0, 0, 0 );
    Vector3D orig_pos = patch_data->vertex_by_index(freeVertexIndex);
    delta_3[j] = delta;
    patch_data->move_vertex( delta_3, freeVertexIndex, err ); 
 
    //compute the function at the perturbed point location
    valid = valid && this->compute_function(patch_data, fdelta, err);  MSQ_ERRZERO(err);

    //compute the numerical gradient
    for (size_t i=0;i<qmHandles.size();i++) {
       mGradient[i][j] = (fdelta[i] - func[i])/delta;
       // MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  Gradient value[%d][%d]=%g\n",i,j,mGradient[i][j]);});
    }

    // put the coordinates back where they belong
    patch_data->set_vertex_coordinates( orig_pos, freeVertexIndex, err );
  }

#else
  double value;
  gradient_out.resize( qmHandles.size() );
  for (size_t i = 0; i < qmHandles.size(); i++) {
    valid = valid &&
         currentQM->evaluate_with_gradient(*patch_data, qmHandles[i], value, tmpIdx, tmpGrad, err); MSQ_ERRZERO(err);
    assert(tmpIdx.size() == 1 && tmpIdx[0] == freeVertexIndex);
    gradient_out[i] = tmpGrad[0];
  }
#endif

  return valid;
}

void NonSmoothDescent::find_active_set( const std::vector<double>& function,
                                        ActiveSet& active_set )
{ 
  // local parameter initialization
  const double activeEpsilon = .3e-4;
  //  activeEpsilon = .3e-8;

    double      function_val;
    double      active_value0;
    double      temp;

//    FUNCTION_TIMER_START("Find Active Set");
    MSQ_DBGOUT(2) << "\nFinding the active set\n";

    /* the first function value is our initial active value */
    active_set.set( 0 );
    active_set.true_active_value = function[0]; 
    //    MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  function value[0]: %g\n",function[0]);});

    /* first sort out the active set... 
       all vals within active_epsilon of largest val */

    for (size_t i = 1; i < qmHandles.size(); i++) {
	function_val = function[i];
        if (active_set.true_active_value < function_val)
          active_set.true_active_value = function_val;
	active_value0 = function[active_set.active_ind[0]];
	temp = fabs(function_val - active_value0);
	//        MSQ_DEBUG_ACTION(3,{fprintf(stdout,"  function value[%d]: %g\n",i,function[i]);});
	if ( function_val > active_value0 ) {  // seek max_i function[i]
	    if ( temp >= activeEpsilon) {
                active_set.set( i );   // new max
	    } 
            else {
                active_set.add( i, fabs(function_val - active_value0) < EPSILON );
	    }
	} else {
	    if (temp < activeEpsilon) {
                active_set.add( i, fabs(function_val - active_value0) < EPSILON );
	    }
	}
    }

}


bool NonSmoothDescent::validity_check(PatchData& pd, MsqError &err)
        
{
//  FUNCTION_TIMER_START("Validity Check");

  // ONLY FOR SIMPLICIAL MESHES - THERE SHOULD BE A VALIDITY CHECKER ASSOCIATED
  // WITH MSQ ELEMENTS
  
  /* check that the simplicial mesh is still valid, based on right handedness. 
       Returns a 1 or a 0 */

  // TODO as a first step we can switch this over to the function
  // evaluation and use the rest of the code as is
  bool valid = true;
  double dEps = 1.e-13;

  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  const MsqVertex* coords = pd.get_vertex_array(err);

  for (size_t i=0;i<pd.num_elements();i++)
  {
    const size_t* conn = pd.element_by_index(i).get_vertex_index_array();
    coords[conn[0]].get_coordinates(x1, y1, z1);
    coords[conn[1]].get_coordinates(x2, y2, z2);
    coords[conn[2]].get_coordinates(x3, y3, z3);
    coords[conn[3]].get_coordinates(x4, y4, z4);

    double dDX2 = x2 - x1;
    double dDX3 = x3 - x1;
    double dDX4 = x4 - x1;

    double dDY2 = y2 - y1;
    double dDY3 = y3 - y1;
    double dDY4 = y4 - y1;

    double dDZ2 = z2 - z1;
    double dDZ3 = z3 - z1;
    double dDZ4 = z4 - z1;

      /* dDet is proportional to the cell volume */
    double dDet = dDX2*dDY3*dDZ4 + dDX3*dDY4*dDZ2 + dDX4*dDY2*dDZ3
      - dDZ2*dDY3*dDX4 - dDZ3*dDY4*dDX2 - dDZ4*dDY2*dDX3 ;

      /* Compute a length scale based on edge lengths. */
    double dScale = ( sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) +
                           (z1-z2)*(z1-z2)) +
                      sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) +
                           (z1-z3)*(z1-z3)) +
                      sqrt((x1-x4)*(x1-x4) + (y1-y4)*(y1-y4) +
                           (z1-z4)*(z1-z4)) +
                      sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3) +
                           (z2-z3)*(z2-z3)) +
                      sqrt((x2-x4)*(x2-x4) + (y2-y4)*(y2-y4) +
                           (z2-z4)*(z2-z4)) +
                      sqrt((x3-x4)*(x3-x4) + (y3-y4)*(y3-y4) +
                           (z3-z4)*(z3-z4)) ) / 6.;

      /* Use the length scale to get a better idea if the tet is flat or
         just really small. */
    if (fabs(dScale) < EPSILON)
    {
      return false;
    }
    else
    {
      dDet /= (dScale*dScale*dScale);
    }

    if (dDet > dEps)
    {
      valid = true;
    }
    else if (dDet < -dEps)
    {
      return false;
    }
    else
    {
      return false;
    }
  }  // end for i=1,numElements
  
  //  MSQ_DEBUG_ACTION(2,{fprintf(stdout,"Mesh Validity is: %d \n",valid);});
  
//  FUNCTION_TIMER_END();
  return(valid);
}


void NonSmoothDescent::get_active_directions( const std::vector<Vector3D>& mGradient, 
                                              std::vector<Vector3D>& dir,
                                              MsqError &/*err*/)
{
    const size_t num_active = mActive.active_ind.size();
    dir.resize(num_active);
    for (size_t i = 0; i < num_active; i++) {
      dir[i] = mGradient[mActive.active_ind[i]];
    }
}


bool NonSmoothDescent::check_vector_dots( const std::vector<Vector3D>& vec,
                                          const Vector3D& normal,
                                          MsqError &/*err*/)
{
    double test_dot = vec[0] % normal;
    unsigned ind = 1;
    while ( (fabs(test_dot) < EPSILON) && (ind < vec.size()) ) {
      test_dot = vec[ind] % normal;
      ind++;
    }
      
    for (unsigned i=ind;i<vec.size();i++) {
       double dot = vec[i] % normal;
       if ( ((dot>0 && test_dot<0) || (dot<0 && test_dot>0)) &&
            (fabs(dot)>EPSILON)) {
          return true;

       }
    }
    return false;
}



bool NonSmoothDescent::convex_hull_test(const std::vector<Vector3D>& vec, MsqError &err)
{
//    int ierr;
    bool equil = false;
    Direction dir_done;
    Status status = MSQ_CHECK_Z_COORD_DIRECTION;
    Vector3D pt1, pt2, pt3, normal;

    /* tries to determine equilibrium for the 3D case */

    if (vec.size() <= 2) status = MSQ_NO_EQUIL;

    while ((status != MSQ_EQUIL) && (status != MSQ_NO_EQUIL) && 
           (status != MSQ_HULL_TEST_ERROR)) {
       if (status == MSQ_CHECK_Z_COORD_DIRECTION) {
          this->find_plane_points(MSQ_ZDIR, MSQ_YDIR, 
                          vec, pt1, pt2, pt3, status, err);
          dir_done = MSQ_ZDIR;
       }
       if (status == MSQ_CHECK_Y_COORD_DIRECTION) {
          this->find_plane_points(MSQ_YDIR, MSQ_XDIR, 
                          vec, pt1, pt2, pt3, status, err);
          dir_done = MSQ_YDIR;
       }
       if (status == MSQ_CHECK_X_COORD_DIRECTION) {
          this->find_plane_points(MSQ_XDIR, MSQ_ZDIR, 
                          vec, pt1, pt2, pt3, status, err);
          dir_done = MSQ_XDIR;
       }
       if (status == MSQ_TWO_PT_PLANE) {
          pt3 = Vector3D(0,0,0);
       }
       if ((status == MSQ_TWO_PT_PLANE) || (status == MSQ_THREE_PT_PLANE)){
           this->find_plane_normal(pt1,pt2,pt3,normal,err); 
           equil = this->check_vector_dots(vec,normal,err); 
           if (equil) {
             switch(dir_done){
             case MSQ_ZDIR:
               equil = false; status = MSQ_CHECK_Y_COORD_DIRECTION;
               break;
             case MSQ_YDIR:
               equil = false; status = MSQ_CHECK_X_COORD_DIRECTION;
               break;
             case MSQ_XDIR:
               equil = true; status = MSQ_EQUIL;
             }
           } else if (equil == 0) {
               status = MSQ_NO_EQUIL;
           } else {
              MSQ_SETERR(err)("equil flag not set to 0 or 1",MsqError::INVALID_STATE);
           }
       }
    }
    switch (status){
    case MSQ_NO_EQUIL:
      MSQ_PRINT(3)("Not an equilibrium point\n");
      equil = false; 
      break;
    case MSQ_EQUIL:
      MSQ_PRINT(3)("An equilibrium point\n");
      equil = true;
      break;
    default:
      MSQ_PRINT(3)("Failed to determine equil or not; status = %d\n",status);
    }
//    FUNCTION_TIMER_END();
    return (equil);
}

void NonSmoothDescent::form_grammian(const std::vector<Vector3D>& vec, MsqError &err)
{
   /* form the grammian with the dot products of the gradients */
   const size_t num_active = mActive.active_ind.size();
   mG.resize( num_active );
   for (size_t i = 0; i < num_active; i++) 
      for (size_t j = i; j < num_active; j++) 
         mG(i,j) = vec[i] % vec[j];
}

bool NonSmoothDescent::check_equilibrium(OptStatus& status, MsqError &err)
{
    std::vector<Vector3D> dir;
 
    //TODO - this subroutine is no longer clear to me... is it still
    // appropriate for quads and hexes?  I think it might be in 2D, but
    // 3D is less clear.  Is there a more general algorithm to use?
    // ask Todd/check in numerical optimization

    bool equil = false;
    const size_t num_active = mActive.active_ind.size();;

    if (num_active == qmHandles.size())
    {
         equil = true; 
         status = MSQ_EQUILIBRIUM;
         MSQ_PRINT(3)("All the function values are in the active set\n"); 
    }

    /* set up the active mGradient directions */
    this->get_active_directions(mGradient,dir,err); MSQ_ERRZERO(err);

    /* normalize the active directions */
    for (size_t j=0;j<num_active;j++)
      dir[j] /= dir[j].length();
 
    equil = this->convex_hull_test(dir,err);
    if (equil) 
      status = MSQ_EQUILIBRIUM;
    
    return equil;
}


static double condition3x3(const double A[9]) 
{
//   int ierr;
   double a11, a12, a13;
   double a21, a22, a23;
   double a31, a32, a33;
//   double s1, s2, s4, s3, t0;
   double s1, s2, s3;
   double denom;
//   double one = 1.0;
   double temp;
   bool zero_denom = true;

   a11 = A[0]; a12=A[1]; a13=A[2];
   a21 = A[3]; a22=A[4]; a23=A[5];
   a31 = A[6]; a32=A[7]; a33=A[8];


   denom = -a11*a22*a33+a11*a23*a32+a21*a12*a33-a21*a13*a32-
            a31*a12*a23+a31*a13*a22;

   if ( (fabs(a11) > EPSILON) && 
        (fabs(denom/a11) > EPSILON)) {
         zero_denom = false;
   }
   if ( (fabs(a22) > EPSILON) && 
        (fabs(denom/a22) > EPSILON)) {
         zero_denom = false;
   }       
   if ( (fabs(a33) > EPSILON) && 
        (fabs(denom/a33) > EPSILON)) {
         zero_denom = false;
   }

   if (zero_denom) {
     return HUGE_VAL;
   } 
   else {
     s1 = sqrt(a11*a11 + a12*a12 + a13*a13 + 
               a21*a21 + a22*a22 + a23*a23 + 
               a31*a31 + a32*a32 + a33*a33);


     temp = (-a22*a33+a23*a32)/denom;
     s3 = temp*temp;
     temp =(a12*a33-a13*a32)/denom;
     s3 += temp*temp;
     temp = (a12*a23-a13*a22)/denom;
     s3 += temp*temp;
     temp = (a21*a33-a23*a31)/denom;
     s3 += temp*temp;
     temp = (a11*a33-a13*a31)/denom;
     s3 += temp*temp;
     temp = (a11*a23-a13*a21)/denom;
     s3 += temp*temp;
     temp = (a21*a32-a22*a31)/denom;
     s3 += temp*temp;
     temp = (-a11*a32+a12*a31)/denom;
     s3 += temp*temp;
     temp = (-a11*a22+a12*a21)/denom;
     s3 += temp*temp;


     s2 = sqrt(s3);
     return s1*s2;
   }
}

double NonSmoothDescent::SymmetricMatrix::condition3x3() const
{
  double values[9] = { 
    operator()(0,0), operator()(0,1), operator()(0,2),
    operator()(1,0), operator()(1,1), operator()(1,2),
    operator()(2,0), operator()(2,1), operator()(2,2)
  };
  return Mesquite::condition3x3( values );
}

void NonSmoothDescent::singular_test(int n, const Matrix3D&  A, bool& singular, MsqError &err) 
{
//    int test;
//    double determinant;
    double cond;

    if ((n>3) || (n<1)) {
      MSQ_SETERR(err)("Singular test works only for n=1 to n=3",MsqError::INVALID_ARG);
      return;
    }

    singular = true;
    switch(n) {
    case 1:
        if (A[0][0] > 0) singular = false;
        break;
    case 2:
        if (fabs(A[0][0]*A[1][1] - A[0][1]*A[1][0]) > EPSILON)           
            singular = false;
        break;
    case 3:
       /* calculate the condition number */
        cond = condition3x3(A[0]); 
        if (cond < 1E14) singular=false;
        break;
    }
}


void NonSmoothDescent::form_PD_grammian(MsqError &err)
{
    int  i,j,k;
    int  g_ind_1;
    bool  singular = false;

    const int num_active = mActive.active_ind.size();
        
    /* this assumes the grammian has been formed */
    for (i=0;i<num_active;i++) {
      for (j=i;j<num_active;j++) {
        if (mG(i,j)==-1) {
          MSQ_SETERR(err)("Grammian not computed properly",MsqError::INVALID_STATE);
          return;
        }
      }
    }

    /* use the first gradient in the active set */
    g_ind_1 = 0;
    mPDG[0][0] = mG(0,0);
    pdgInd[0] = mActive.active_ind[0];

    /* test the rest and add them as appropriate */
    k = 1; i = 1;
    while( (k<3) && (i < num_active) ) {
        mPDG[0][k] = mPDG[k][0] = mG(0,i);
        mPDG[k][k] = mG(i,i);
        if ( k == 2) { /* add the dot product of g1 and g2 */
           mPDG[1][k] = mPDG[k][1] = mG(g_ind_1,i);
        }
        this->singular_test(k+1,mPDG,singular,err);
        if (!singular) {
           pdgInd[k] = mActive.active_ind[i];
           if (k==1) g_ind_1 = i;
           k++;
        }
        i++;
    }
}


void NonSmoothDescent::search_edges_faces( const Vector3D* dir, Vector3D& mSearch, MsqError &err)
{
    bool viable;
    double a,b,denom;
    Vector3D g_bar;
    Vector3D temp_search(0,0,0); /* initialize the search direction to 0,0 */
    double projection, min_projection;

    const size_t num_active = mActive.active_ind.size();;

    /* Check for viable faces */
    min_projection = HUGE_VAL;
    for (size_t i=0; i<num_active; i++) {
        /* FACE I */
        viable = true;

        /* test the viability */
        for (size_t j=0;j<num_active;j++) {       /* lagrange multipliers>0 */
             if (mG(j,i) < 0) 
                  viable = false;
        }
       
        /* find the minimum of viable directions */
        if ((viable) && (mG(i,i) < min_projection)) {
            min_projection = mG(i,i);
            temp_search = dir[i];
            mSteepest = mActive.active_ind[i];
        }
    
       /* INTERSECTION IJ */
       for (size_t j=i+1; j<num_active; j++) {
          viable = true;

          /* find the coefficients of the intersection 
             and test the viability */
          denom = 2*mG(i,j) - mG(i,i) - mG(j,j);
          a = b = 0;
          if (fabs(denom) > EPSILON) {
             b = (mG(i,j) - mG(i,i))/denom;
             a = 1 - b;
             if ((b < 0) || (b > 1)) /* 0 < b < 1 */
                 viable = false;  
	     for (size_t k=0;k<num_active;k++) {       /* lagrange multipliers>0 */
                 if ((a*mG(k,i) + b*mG(k,j)) <= 0) 
                     viable = false;
             }
          } else {
             viable = false;                        /* Linearly dependent */
          }

          /* find the minimum of viable directions */
          if (viable) {
             g_bar = a * dir[i] + b * dir[j];
             projection = g_bar % g_bar;
             if (projection < min_projection) {
	        min_projection = projection;
                temp_search = g_bar;
                mSteepest = mActive.active_ind[i];
             }
          }
       }
    }
    if (optStatus != MSQ_EQUILIBRIUM) {
        mSearch = temp_search;
    }
}         


 bool NonSmoothDescent::solve2x2( double a11, double a12,
                                  double a21, double a22, 
                                  double b1, double b2,
                                  double x[2], MsqError &/*err*/)
{
    double factor;

    /* if the system is not singular, solve it */
    if (fabs(a11*a22 - a21*a12) > EPSILON) {
	if (fabs(a11) > EPSILON) {
	    factor = (a21/a11);
	    x[1] = (b2 - factor*b1)/(a22 - factor*a12);
	    x[0] = (b1 - a12*x[1])/a11;
	} else if (fabs(a21) > EPSILON) {
	    factor = (a11/a21);
	    x[1] = (b1 - factor*b2)/(a12 - factor*a22);
	    x[0] = (b2 - a22*x[1])/a21;
	}
        return true;
    } else {
	return false;
    }
}


void NonSmoothDescent::form_reduced_matrix(SymmetricMatrix& P)
{
    const size_t P_size = mActive.active_ind.size() - 1;
    P.resize( P_size );
    for (size_t i = 0; i < P_size; i++) {
        P(i,i) = mG(0,0) - 2*mG(0,i+1) + mG(i+1,i+1);
        for (size_t j = i+1; j < P_size; j++) {
            P(i,j) = mG(0,0) - mG(0,j+1) - mG(i+1,0) + mG(i+1,j+1);
        }
    }
}


void  NonSmoothDescent::get_min_estimate( double *final_est, MsqError &/*err*/)
{
    double est_imp;

    *final_est = -HUGE_VAL;
    for (size_t i=0;i<mActive.active_ind.size();i++) {
	est_imp = -mAlpha*mGS[mActive.active_ind[i]];
        if (est_imp>*final_est) *final_est = est_imp;
    }
    if (*final_est == 0) {
	*final_est = -HUGE_VAL;
	for (size_t i=0;i<qmHandles.size();i++) {
	    est_imp = -mAlpha*mGS[i];
	    if ((est_imp>*final_est) && (fabs(est_imp) > EPSILON)) {
		*final_est = est_imp;
	    }
	}
    }
}


void NonSmoothDescent::get_gradient_projections(const Vector3D& mSearch, MsqError &/*err*/)
{
    for (size_t i=0;i<qmHandles.size();i++) 
	mGS[i] = mGradient[i] % mSearch;

    MSQ_PRINT(3)("steepest in get_gradient_projections %lu\n",(unsigned long)mSteepest);
}


void NonSmoothDescent::compute_alpha(MsqError &/*err*/)
{
    double    steepest_function;
    double    steepest_grad;
    double    alpha_i;
    double    min_positive_value=HUGE_VAL;

//    FUNCTION_TIMER_START("Compute Alpha");

    MSQ_PRINT(2)("In compute alpha\n");

    mAlpha = HUGE_VAL;

    steepest_function = mFunction[mSteepest];
    steepest_grad = mGS[mSteepest];
    for (size_t i=0;i<qmHandles.size();i++)
    {
        /* if it's not active */
      if (i!=mSteepest)
      {
	  alpha_i = steepest_function - mFunction[i];
	   
	  if (fabs(mGS[mSteepest] - mGS[i])>1E-13) {
	     /* compute line intersection */
	     alpha_i = alpha_i/(steepest_grad - mGS[i]);
	  } else {
	     /* the lines don't intersect - it's not under consideration*/
	     alpha_i = 0;
	  }
	  if ((alpha_i > minStepSize ) && (fabs(alpha_i) < fabs(mAlpha))) {
	    mAlpha = fabs(alpha_i); 
            MSQ_PRINT(3)("Setting alpha %lu %g\n",(unsigned long)i,alpha_i);
	  }
          if ((alpha_i > 0) && (alpha_i < min_positive_value)) {
            min_positive_value = alpha_i;
          }
       }
    }

    if ((mAlpha == HUGE_VAL) && (min_positive_value != HUGE_VAL)) {
      mAlpha = min_positive_value;
    }

    /* if it never gets set, set it to the default */
    if (mAlpha == HUGE_VAL) {
      mAlpha = maxAlpha;
      MSQ_PRINT(3)("Setting alpha to the maximum step length\n");
    }

    MSQ_PRINT(3)("  The initial step size: %f\n",mAlpha);

//    FUNCTION_TIMER_END();
}


void NonSmoothDescent::print_active_set( const ActiveSet& active_set, const std::vector<double>& func )
{
    if (active_set.active_ind.empty()) MSQ_DBGOUT(3)<< "No active values\n";
    /* print the active set */
    for (size_t i = 0; i < active_set.active_ind.size(); i++) {
     MSQ_PRINT(3)("Active value %lu:   %f \n", (unsigned long)i+1,func[active_set.active_ind[i]]);
    }
}


void NonSmoothDescent::init_opt(PatchData& pd, MsqError &err)
{
    qmHandles.clear();
    currentQM->get_evaluations( pd, qmHandles, true, err ); MSQ_ERRRTN(err);

    MSQ_PRINT(2)("\nInitializing Optimization \n");

    /* for the purposes of initialization will be set to zero after */
    optStatus = MSQ_NO_STATUS;
    mSteepest = 0;
    mAlpha = 0;
    maxAlpha = 0;

    MSQ_PRINT(3)("  Initialized Constants \n");
    pdgInd[0] = pdgInd[1] = pdgInd[2] = -1;
    mPDG = Matrix3D(0.0);

    MSQ_PRINT(3)("  Initialized search and PDG \n");
    mFunction.clear(); 
    mFunction.resize( qmHandles.size(), 0.0 );
    testFunction.clear(); 
    testFunction.resize( qmHandles.size(), 0.0 );
    originalFunction.clear(); 
    originalFunction.resize( qmHandles.size(), 0.0 );
    mGradient.clear();
    mGradient.resize( qmHandles.size(), Vector3D(0.0) );
    mGS.resize( qmHandles.size() );

    MSQ_PRINT(3)("  Initialized function/gradient \n");
    mG.resize( qmHandles.size() );
    mG.fill(-1);
    MSQ_PRINT(3)("  Initialized G\n");
 
    prevActiveValues.clear();
    prevActiveValues.reserve(32);
    MSQ_PRINT(3)("  Initialized prevActiveValues\n");
}


void NonSmoothDescent::init_max_step_length( PatchData& pd, MsqError &err )
{
  size_t i, j;
  double max_diff = 0;
  double diff=0;

  MSQ_PRINT(2)("In init_max_step_length\n");

  /* check that the input data is correct */
  if (pd.num_elements()==0) {
    MSQ_SETERR(err)("Num incident vtx = 0\n",MsqError::INVALID_MESH);
    return;
  }

  /* find the maximum distance between two incident vertex locations */
  const MsqVertex* coords = pd.get_vertex_array(err);
  for (i=0;i<pd.num_nodes()-1;i++) {
    for (j=i;j<pd.num_nodes();j++) {
      diff = (coords[i]-coords[j]).length_squared();
      if (max_diff < diff) max_diff=diff;
    } 
  }
  max_diff = sqrt(max_diff);
  if (max_diff==0) {
     MSQ_SETERR(err)("Maximum distance between incident vertices = 0\n",
                     MsqError::INVALID_MESH);
    return;
  }
  maxAlpha = max_diff/100;

  MSQ_PRINT(3)("  Maximum step is %g\n",maxAlpha);
}

} // namespace Mesquite

