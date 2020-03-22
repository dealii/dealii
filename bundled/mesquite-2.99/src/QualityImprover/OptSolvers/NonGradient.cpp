/****************************************************************** 
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
   
  ******************************************************************/
/*!
  \file   NonGradient.cpp
  \brief  
  The NonGradient class is also a concrete vertex mover
  which performs derivative free minimization
  based on the Amoeba Method, as implemented in Numerical Recipes in C.

  The other optimization methods either accept or reject a point immediately.
  Amoeba may not accept a point for several iterations/evaluations.
  Amoeba does not check for BOUNDED_VERTEX_MOVERMENT. 
  rtol=2.0*fabs( height[ihi]-height[ilo] )/
         ( fabs(height[ihi])+fabs(height[ilo])+threshold );

  We are aware of the problem of choosing the initial simplex, but have not
  sovled it.  How has this been solved?
  There is a problem of letting Mesquite handle convergence.  The idea here is to 
  call accumulate_? only if the best point changes. But does it matter?  Where 
  are the errors determined?  What does it mean to cull?
  What does Pat... compare and contrast the differnt vertex movement criteria
  with the Amoeba criterion.
*/

#include "NonGradient.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"
#include "MsqError.hpp"
#include <cmath>
#include <iostream>
#include "../../ObjectiveFunction/ObjectiveFunction.hpp"
#include <float.h>

namespace MESQUITE_NS {

std::string NonGradient::get_name() const
  { return "NonGradient"; }
  
PatchSet* NonGradient::get_patch_set()
  { return PatchSetUser::get_patch_set(); }

NonGradient::NonGradient(ObjectiveFunction* of) 
  : VertexMover(of),
    PatchSetUser(true),
    projectGradient(false),
    mDimension(0),
    mThreshold(0.0),
    mTolerance(0.0),
    mMaxNumEval(0),
    mNonGradDebug(0),
    mUseExactPenaltyFunction(true),
    mScaleDiameter(0.1)
{
  set_debugging_level(2);
  //set the default inner termination criterion
  TerminationCriterion* default_crit=get_inner_termination_criterion();
  if(default_crit==NULL){
    return;
  }
  else{
    default_crit->add_iteration_limit( 5 );
  }
}  

NonGradient::NonGradient(ObjectiveFunction* of, MsqError &err) 
  : VertexMover(of),
    PatchSetUser(true),
    projectGradient(false),
    mDimension(0),
    mThreshold(0.0),
    mTolerance(0.0),
    mMaxNumEval(0),
    mNonGradDebug(0),
    mUseExactPenaltyFunction(true),
    mScaleDiameter(0.1)
{
  set_debugging_level(2);
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

bool
NonGradient::testRowSum( int numRow, int numCol, double* matrix, double* oldRowSum)
{
  bool same = true;
  std::vector<double> rowSum(numRow,0.);
  double maxVal = 0.;
  for (int col=0;col<numCol;col++)    
  {
    for (int row=0;row<numRow;row++)
    {
      rowSum[row] += matrix[row+col*numRow];
      if( fabs( matrix[row+col*numRow] ) > maxVal )
      {
         maxVal = fabs( matrix[row+col*numRow] );
      }
    }
  }
  double machEps = 1.e-14 * static_cast<double>(numRow); // better to use system parameters
  double upperBound = machEps * maxVal + 1.e-15;
  for (int row=0;row<numRow;row++)
  {
    if( fabs( rowSum[row] -  oldRowSum[row]) >  upperBound ) 
    { 
         same = false;
         if( mNonGradDebug >= 2 )
         {
           std::cout << " Warning: NonGradient Row Sum " << row << " Test: value " << rowSum[row] << " Discrepancy " <<  rowSum[row] -  oldRowSum[row] << " maxVal " << maxVal << " numRow " << numRow << " numCol " << numCol <<std::endl;
         }
         MSQ_PRINT(2)("NonGradient Row Sum [%d] Test failure: value %22.15e  difference %22.15e \n", row, rowSum[row], rowSum[row] -  oldRowSum[row]); 
    }
  }
  return(same);
}

void 
NonGradient::getRowSum( int numRow, int numCol, std::vector<double>& matrix, std::vector<double>& rowSum)
{
  for (int row=0;row<numRow;row++)
  {
      rowSum[row] = 0.;
  }
  for (int col=0;col<numCol;col++)    
  {
    for (int row=0;row<numRow;row++)
    {
      rowSum[row] += matrix[row+col*numRow];
    }
  }
}

double 
NonGradient::evaluate( double *point,  PatchData &pd, MsqError &err )
{
  if( pd.num_free_vertices() > 1 )
  {
    MSQ_SETERR(err)("Only one free vertex per patch implemented", MsqError::NOT_IMPLEMENTED);
  }
  const size_t vertexIndex = 0; 
  Vector3D originalVec = pd.vertex_by_index(vertexIndex);
  Vector3D pointVec;
  for( int index = 0; index<3; index++)
  {
    pointVec[index] = point[index];
  }
  pd.set_vertex_coordinates( pointVec, vertexIndex, err ); 
  pd.snap_vertex_to_domain( vertexIndex, err );  

  OFEvaluator& obj_func = get_objective_function_evaluator();
  TerminationCriterion* term_crit=get_inner_termination_criterion();

  double value;
  bool feasible = obj_func.evaluate( pd, value, err ); // MSQ_ERRZERO(err);
  term_crit->accumulate_inner( pd, value, 0, err );  //MSQ_CHKERR(err);
  if (err.error_code() == err.BARRIER_VIOLATED)
    err.clear();   // barrier violation not an error here
  MSQ_ERRZERO(err);
  pd.set_vertex_coordinates( originalVec, vertexIndex, err ); 
  if( !feasible && mUseExactPenaltyFunction )  
  { // "value" undefined btw
    double ensureFiniteRtol= .25;
    value = DBL_MAX * ensureFiniteRtol;
    //std::cout << " Infeasible patch: " << value << std::endl; printPatch( pd, err );
  }
  return value;
}

// Extrapolate by a factor of fac through the face of the simplex
// opposite the high point to a new point.  If the new point is 
// better, swap it with the high point.
double
NonGradient::amotry( std::vector<double>& simplex, 
                 std::vector<double>& height,
                 double psum[], int ihi, double fac, PatchData &pd, MsqError &err)
{
  int numRow = getDimension();
  int numCol = numRow + 1;
  std::vector<double> ptry(numRow); // does this make sense?
  double fac1=(1.0-fac)/static_cast<double>(numRow);
  double fac2=fac1-fac;
  for (int row=0;row<numRow;row++)
  {
    ptry[row]=psum[row]*fac1-simplex[row+ihi*numRow]*fac2;
  }
  if( mNonGradDebug >= 3 ) 
  {      
    std::cout << "Try ";
  }      
  MSQ_PRINT(3)("Try");
  double ytry = evaluate(&ptry[0], pd, err); // value at trial point
  if( mNonGradDebug >= 3 ) 
  {      
    std::cout << ytry << std::endl;
  }      
  MSQ_PRINT(3)("yTry");
  if (ytry < height[ihi]) // better than highest (worst)
  {
    height[ihi]=ytry;     // swap ihi and ytry
    for (int row=0;row<numRow;row++)
    {
      psum[row] += (ptry[row]-simplex[row+ihi*numRow]);
      simplex[row+ihi*numRow]=ptry[row];
    }
  }
  return ytry;
}

void NonGradient::printPatch(const PatchData &pd, MsqError &err)
{
  if( mNonGradDebug==0 )
  {
    return;
  }
  const size_t numNode = pd.num_nodes();   //27,  27 what?
  MSQ_PRINT(3)("Number of Vertices: %d\n",(int)pd.num_nodes());
  const size_t numVert = pd.num_free_vertices(); // 1
  MSQ_PRINT(3)("Num Free = %d\n",(int)pd.num_free_vertices());
  const size_t numSlaveVert = pd.num_slave_vertices(); //0
  const size_t numCoin = pd.num_corners(); // 64
  const MsqVertex* coord = pd.get_vertex_array(err);
  MSQ_PRINT(3)("Number of Vertices: %d\n",(int)pd.num_nodes());

  std::cout << "Patch " << numNode << "  " << numVert << "  " << numSlaveVert << "  " << numCoin << std::endl;
  MSQ_PRINT(3)("");
  std::cout << "Coordinate ";
  std::cout << "         " << std::endl;
  for( size_t index = 0; index < numVert; index++ )
  {
    std::cout << coord[index][0] << "  " << coord[index][1] << "  " << coord[index][2] << std::endl;
  }
  //const size_t numElt = pd.num_elements();
  if( mNonGradDebug >= 3 ) 
  {      
         std::cout << "Number of Elements: " << pd.num_elements() << std::endl;
  }      
  MSQ_PRINT(3)("Number of Elements: %d\n",(int)pd.num_elements());
}

void
NonGradient::printSimplex( std::vector<double>& simplex , std::vector<double>& height )
{
  int numRow = getDimension();
  int numCol = numRow + 1;
  for (int col=0;col<numCol;col++)    
  {
    //std::cout << "simplex[ " << col << "]= " ;
    for (int row=0;row<numRow;row++)
    {
      std::cout << simplex[row+col*numRow] << " "; 
    }
    //std::cout << "           "  << height[col] << std::endl;
  }
  for (int col=0;col<numCol;col++)    
  {
    std::cout << height[col] << " ";
  }
  std::cout << std::endl;
}

int NonGradient::getPatchDimension(const PatchData &pd, MsqError &err)
{
  const size_t numElt = pd.num_elements();
  unsigned edimMax = 0;  // in case numElt == 0
  for( size_t elementId = 0; elementId < numElt; elementId++)
  {
    const MsqMeshEntity& element = pd.element_by_index( elementId );
    EntityTopology type = element.get_element_type();
    unsigned edim = TopologyInfo::dimension( type );
    if( elementId == 0 ) 
    {
      edimMax = edim;
    } 
    else
    { 
      if(edimMax != edim)
      {
        MSQ_SETERR(err)("A patch has elements of different dimensions", MsqError::INVALID_MESH);
        std::cout << "A patch has elements of different dimensions" << edimMax << "  "  << edim << std::endl;
        if(edimMax < edim)
        {
          edimMax = edim;
        }
      }
    }
  }
  return( edimMax );
}
  
void NonGradient::initialize(PatchData &/*pd*/, MsqError &/*err*/)
{
}

void NonGradient::initialize_mesh_iteration(PatchData &pd, MsqError &err)
{
  int elementDimension = getPatchDimension( pd, err );  // to do: react to error
  int dimension = elementDimension * pd.num_free_vertices();
  //printPatch( pd, err );
  setDimension(dimension);
  int maxNumEval = 100*dimension;  // 1. Make this a user parameter
  setMaxNumEval(maxNumEval);
  double threshold = 1.e-10; // avoid division by zero
  setThreshold(threshold);
  double minEdgeLen = 0.0;
  double maxEdgeLen = 0.0;
//  double ftol = 0.;
  if( dimension > 0 )
  {
    pd.get_minmax_edge_length( minEdgeLen, maxEdgeLen );
    //ftol = minEdgeLen * 1.e-4; // Turn off Amoeba convergence criterion
    if( mNonGradDebug >= 1 ) 
    {      
         std::cout << "minimum edge length " << minEdgeLen << " maximum edge length " << maxEdgeLen << std::endl;
    }      
    MSQ_PRINT(3)("minimum edge length %e    maximum edge length %e\n", minEdgeLen,  maxEdgeLen);
  }
//  setTolerance(ftol);
  int numRow = dimension;
  int numCol = numRow+1;  
  if( numRow*numCol <= simplex.max_size() )
  { 
    simplex.assign(numRow*numCol, 0.);  // guard against previous simplex value
    double scale = minEdgeLen * mScaleDiameter;; 
    const MsqVertex* coord = pd.get_vertex_array(err);
    if( pd.num_free_vertices() > 1 )
    {
      MSQ_SETERR(err)("Only one free vertex per patch implemented", MsqError::NOT_IMPLEMENTED);
    }
    size_t index = 0;
    for( int col = 0; col < numCol; col++ )
    {
      for (int row=0;row<numRow;row++)
      {
        simplex[ row + col*numRow ] = coord[index][row];
        if( row == col-1 )
        {
          simplex[ row + col*numRow ] += scale/ static_cast<double>(numCol);
        }
      }
    }
  }
  else
  {
    MSQ_SETERR(err)("Use patch with fewer free vertices", MsqError::OUT_OF_MEMORY);
    if( mNonGradDebug >= 1 ) 
    {      
      std::cout << "ERROR: Too many free vertices in patch" << std::endl;
    }      
    //MSQ_PRINT(1)("ERROR: Too many free vertices in patch\n");
  }
}

void NonGradient::optimize_vertex_positions(PatchData &pd, 
                                            MsqError &err)
{
  MSQ_FUNCTION_TIMER( "NonGradient::optimize_vertex_positions" );
  int numRow = getDimension();
  int numCol = numRow+1;  
  std::vector<double> height(numCol); 

  for(int col = 0; col < numCol; col++)
  {
    height[col] =  evaluate(&simplex[col*numRow], pd, err);//  eval patch stuff
  }
  if(mNonGradDebug > 0)
  {
    printSimplex( simplex, height );
  }

  // standardization
  TerminationCriterion* term_crit=get_inner_termination_criterion();
  int maxNumEval = getMaxNumEval();
  double threshold = getThreshold();
//  double ftol = getTolerance();
  int ilo = 0;  //height[ilo]<=...
  int inhi = 0; //...<=height[inhi]<=
  int ihi = 0;  //<=height[ihi] 
//  double rtol = 2.*ftol;
  double ysave;
  double ytry;
  bool afterEvaluation = false;
  std::vector<double> rowSum(numRow);
  getRowSum( numRow, numCol, simplex, rowSum);
  while( !( term_crit->terminate() ) )
  {

    if(mNonGradDebug > 0)
    {
      printSimplex( simplex, height );
    }
    //std::cout << "rtol " << rtol << " ftol " << ftol << " MesquiteIter " << term_crit->get_iteration_count() << " Done " << term_crit->terminate() << std::endl;

    if( afterEvaluation )   
    {
      // reflect highPt through opposite face
      // height[0] may vanish
/*
      if( !testRowSum( numRow, numCol, &simplex[0], &rowSum[0]) )
      {
        // Before uncommenting here and ... 
        //MSQ_SETERR(err)("Internal check sum test A failed", MsqError::INTERNAL_ERROR);
        //MSQ_ERRRTN(err);
      }
*/
      ytry=amotry(simplex,height,&rowSum[0],ihi,-1.0, pd, err);
/*
      if( !testRowSum( numRow, numCol, &simplex[0], &rowSum[0]) )
      {
         // ... here, determine a maxVal majorizing the previous as well as the current simplex.
         //MSQ_SETERR(err)("Internal check sum test B failed", MsqError::INTERNAL_ERROR);
         //MSQ_ERRRTN(err);
      }   
*/
  
/*
      if( height[0] == 0.)
      {
         MSQ_SETERR(err)("(B) Zero objective function value", MsqError::INTERNAL_ERROR);
         exit(-1);
      }
*/
      if (ytry <= height[ilo])   
      {
        ytry=amotry(simplex,height,&rowSum[0],ihi,-2.0,pd,err);
        if( mNonGradDebug >= 3 ) 
        {      
         std::cout << "Reflect and Expand from highPt " << ytry << std::endl;
        }      
        //MSQ_PRINT(3)("Reflect and Expand from highPt : %e\n",ytry);
      }
      else 
      {
        if (ytry >= height[inhi]) 
        {
          ysave=height[ihi]; // Contract along highPt
          ytry=amotry(simplex,height,&rowSum[0],ihi,0.5,pd,err);
          if (ytry >= ysave)
          { // contract all directions toward lowPt
            for (int col=0;col<numCol;col++)
            {
              if (col != ilo)
              {
                for (int row=0;row<numRow;row++)
                {
                  rowSum[row]=0.5*(simplex[row+col*numRow]+simplex[row+ilo*numRow]);
                  simplex[row+col*numRow]=rowSum[row];
                }
                height[col] = evaluate(&rowSum[0], pd, err); 
                if( mNonGradDebug >= 3 ) 
                {      
                  std::cout << "Contract all directions toward lowPt value( " << col << " ) = " << height[col] << " ilo = " << ilo << std::endl;
                }      
                //MSQ_PRINT(3)("Contract all directions toward lowPt value( %d ) = %e    ilo = %d\n", col, height[col], ilo);
              }
            }
          }
        }   
      } // ytri > h(ilo) 
    } // after evaluation
    ilo=1; // conditional operator or inline if 
    ihi = height[0] > height[1] ? (inhi=1,0) : (inhi=0,1);
    for (int col=0;col<numCol;col++)
    {
      if (height[col] <= height[ilo])
      {
        ilo=col;  // ilo := argmin height
      }
      if (height[col] > height[ihi])
      {
        inhi=ihi;
        ihi=col;
      } 
      else  // height[ihi] >= height[col]
        if (col != ihi && height[col] > height[inhi] ) inhi=col;
    }
//    rtol=2.0*fabs( height[ihi]-height[ilo] )/
//         ( fabs(height[ihi])+fabs(height[ilo])+threshold );
    afterEvaluation = true;
  } //  while not converged 

  // Always set to current best mesh.
  { 
    if( ilo != 0 )
    {
      double yTemp = height[0];
      height[0] = height[ilo]; // height dimension numCol
      height[ilo] = yTemp;
      for (int row=1;row<numRow;row++)
      { 
          yTemp = simplex[row];
          simplex[row] = simplex[row+ilo*numRow];
          simplex[row+ilo*numRow] = yTemp;
      }
    }
  }
  if( pd.num_free_vertices() > 1 )
  {
    MSQ_SETERR(err)("Only one free vertex per patch implemented", MsqError::NOT_IMPLEMENTED);
  }

  Vector3D newPoint( &simplex[0] ); 
  size_t vertexIndex = 0; // fix c.f. freeVertexIndex
  pd.set_vertex_coordinates( newPoint, vertexIndex, err ); 
  pd.snap_vertex_to_domain( vertexIndex, err );
  if( term_crit->terminate() )
  {
    if( mNonGradDebug >= 1 ) 
    {      
         std::cout << "Optimization Termination OptStatus: Max Iter Exceeded" << std::endl;
    }      
    //MSQ_PRINT(1)("Optimization Termination OptStatus: Max Iter Exceeded\n"); 
  }
}

void NonGradient::terminate_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
  if( mNonGradDebug >= 2 )
  {
    std::cout << "- Executing NonGradient::iteration_complete()" << std::endl;
  }
  //MSQ_PRINT(2)("\n - Executing NonGradient::iteration_complete() \n");
}
  
void NonGradient::cleanup()
{
  if( mNonGradDebug >= 2 )
  {
    std::cout << " - Executing NonGradient::iteration_end()" << std::endl;
  }
  //MSQ_PRINT(2)("\n - Executing NonGradient::iteration_end() \n");
}
  
} // namespace Mesquite
