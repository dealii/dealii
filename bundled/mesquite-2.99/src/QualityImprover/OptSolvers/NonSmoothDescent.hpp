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
  \file   NonSmoothDescent.hpp
  \brief  

  The NonSmoothDescent Class implements the steepest descent algorythm in
  order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QaulityMetric object.

  \author Thomas Leurent
  \date   2002-06-13
*/

#ifndef Mesquite_NonSmoothDescent_hpp 
#define Mesquite_NonSmoothDescent_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "ObjectiveFunction.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqDebug.hpp"
#include "QualityMetric.hpp"
#include "VertexPatches.hpp"
#include <vector>
#include "Vector3D.hpp"

namespace MESQUITE_NS
{

  class NonSmoothDescent : public VertexMover 
  {
  public:
    MESQUITE_EXPORT
    NonSmoothDescent(QualityMetric* qm);

    MESQUITE_EXPORT virtual
    ~NonSmoothDescent() { }
    
    MESQUITE_EXPORT virtual
    std::string get_name() const;
    
    MESQUITE_EXPORT virtual
    PatchSet* get_patch_set();
    
  protected:


    struct ActiveSet
    {
      int num_equal;
      double true_active_value;
      std::vector<int> active_ind;
      
      void set( int index )
      {
        active_ind.clear();
        active_ind.push_back(index);
        num_equal = 0;
      }
      void add( int index, bool equal )
      {
        active_ind.push_back(index);
        num_equal += equal;
      }
    };

    class SymmetricMatrix
    {
      double* storage;
      size_t size;
      size_t index(size_t r, size_t c) const
        { return (r <= c) ? size*r - r*(r+1)/2+c : size*c - c*(c+1)/2+r; }
      
      public:
      
      SymmetricMatrix() : storage(0), size(0) {}
      
      ~SymmetricMatrix() { free(storage); }
      
      void resize( size_t new_size ) {
        size = new_size;
        storage = (double*)realloc( storage, sizeof(double) * size * (size+1) / 2 );
      }
      
      void fill( double value ) {
        std::fill( storage, storage + (size * (size+1) / 2), value );
      }
      
      double& operator()(size_t r, size_t c)
        { return storage[index(r,c)]; }
        
      double operator()(size_t r, size_t c) const
        { return storage[index(r,c)]; }
        
      double condition3x3() const;
    };

    MESQUITE_EXPORT virtual
    void initialize(PatchData &pd, MsqError &err);
    
    MESQUITE_EXPORT virtual
    void optimize_vertex_positions(PatchData &pd, MsqError &err);
    
    MESQUITE_EXPORT virtual
    void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    
    MESQUITE_EXPORT virtual
    void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    
    MESQUITE_EXPORT virtual
    void cleanup();
    
  private:

    enum OptStatus {
       MSQ_NO_STATUS = 0,
       MSQ_STEP_ACCEPTED = 100,  
       MSQ_IMP_TOO_SMALL,
       MSQ_FLAT_NO_IMP,
       MSQ_STEP_TOO_SMALL,
       MSQ_EQUILIBRIUM,
       MSQ_ZERO_SEARCH,
       MSQ_MAX_ITER_EXCEEDED
    };

    enum Status {
      MSQ_NO_EQUIL = 101,
      MSQ_CHECK_TOP_DOWN,
      MSQ_CHECK_BOTTOM_UP,
      MSQ_TWO_PT_PLANE,
      MSQ_THREE_PT_PLANE,
      MSQ_CHECK_Y_COORD_DIRECTION,
      MSQ_CHECK_X_COORD_DIRECTION,
      MSQ_CHECK_Z_COORD_DIRECTION,
      MSQ_EQUIL,
      MSQ_HULL_TEST_ERROR
    };

    enum Direction {
      MSQ_XDIR = 0,
      MSQ_YDIR = 1,
      MSQ_ZDIR = 2 };


      /* local copy of patch data */
      //    PatchData patch_data;
    size_t freeVertexIndex;
    
      /* smoothing parameters */
    double minStepSize;
    
      /* optimization data */
    VertexPatches patchSet;
    QualityMetric* currentQM;
    double originalValue;
    std::vector<double> mFunction, testFunction, originalFunction;
    std::vector<double> mGS;
    SymmetricMatrix mG;
    Matrix3D mPDG;
    std::vector<double> prevActiveValues;
    std::vector<Vector3D> mGradient, tmpGrad;
    std::vector<size_t> tmpIdx, qmHandles;
    OptStatus optStatus;
    size_t mSteepest;
    double mAlpha;
    double maxAlpha;
    int pdgInd[3];
    ActiveSet mActive, testActive, originalActive;
    
      /* functions */
    void init_opt(PatchData& pd, MsqError &err);
    void init_max_step_length(PatchData& pd, MsqError &err);
    
      /* optimize */
    void minmax_opt(PatchData &pd, MsqError &err);
    void step_acceptance(PatchData &pd, int iter_count, const Vector3D& search_dir, MsqError &err); 
    void get_min_estimate(double *final_est, MsqError &err);
    void get_gradient_projections(const Vector3D& mSearch, MsqError &err);
    void compute_alpha(MsqError &err);
    
      /* function/gradient/active set computations */
    bool compute_function(PatchData *pd, std::vector<double>& function, MsqError &err);
    bool compute_gradient(PatchData *pd, std::vector<Vector3D>& grad_out, MsqError &err);
    void find_active_set( const std::vector<double>& function, ActiveSet& active_set);
    void print_active_set(const ActiveSet& active_set, const std::vector<double>& func);
    
      /* checking validity/improvement */
    bool improvement_check(MsqError &err);
    bool validity_check(PatchData& pd, MsqError &err);
    
      /* checking equilibrium routines */
    bool check_equilibrium(OptStatus& opt_status, MsqError &err);
    bool convex_hull_test(const std::vector<Vector3D>& vec, MsqError &err);
    bool check_vector_dots(const std::vector<Vector3D>& vec, const Vector3D& normal, MsqError &err);
    void find_plane_normal( const Vector3D& pt1,
                            const Vector3D& pt2,
                            const Vector3D& pt3,
                            Vector3D& cross,
                            MsqError &/*err*/ );
    void find_plane_points( Direction dir1, 
                            Direction dir2,
                            const std::vector<Vector3D>& vec,
                            Vector3D& pt1,
                            Vector3D& pt2,
                            Vector3D& pt3,
                            Status& status,
                            MsqError& );
    
      /* from the matrix file */
    void form_grammian(const std::vector<Vector3D>& vec, MsqError &err);
    void form_PD_grammian(MsqError &err);
    void singular_test(int n, const Matrix3D& A, bool& singular, MsqError &err);
    bool solve2x2(double a11, double a12, double a21, double a22, 
                  double b1, double b2, double x[2],MsqError &err);
    void form_reduced_matrix(SymmetricMatrix& P);
    
      /* search direction */
    void search_direction(PatchData &pd, Vector3D& search_dir_out, MsqError &err);
    void search_edges_faces(const Vector3D* dir, Vector3D& search, MsqError &err);
    void get_active_directions( const std::vector<Vector3D>& gradient,
                                std::vector<Vector3D>& dir, 
                                MsqError &err );
    NonSmoothDescent(const NonSmoothDescent &pd); //disable copying
    NonSmoothDescent& operator=(const NonSmoothDescent &pd);  //disable assignment

  };
  

inline void NonSmoothDescent::find_plane_normal( const Vector3D& pt1,
                                                 const Vector3D& pt2,
                                                 const Vector3D& pt3,
                                                 Vector3D& cross,
                                                 MsqError &/*err*/ )
{
  Vector3D vecA = pt2 - pt1;
  Vector3D vecB = pt3 - pt1;
  cross = vecA * vecB;
  cross /= cross.length();
}

inline bool NonSmoothDescent::improvement_check(MsqError &/*err*/)
{
  /* check to see that the mesh didn't get worse */
  if (originalValue < mActive.true_active_value) {
     MSQ_PRINT(2)("The local mesh got worse; initial value %f; final value %f\n",
	       originalValue,  mActive.true_active_value );
       return false;
   }

  return true;
}


} // namespace

#endif  // Mesquite_NonSmoothDescent_hpp 
