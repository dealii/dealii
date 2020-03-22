/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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

    (2006) kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_MAPPING_FUNCTION_HPP
#define MSQ_MAPPING_FUNCTION_HPP

/** \file MappingFunction.hpp
 *  \brief Header containg defintion of MappingFunction
 *  \author Jason Kraftcheck
 */


#include "Mesquite.hpp"
#include <vector>
#include "MsqMatrix.hpp"
#include "TopologyInfo.hpp"
#include "NodeSet.hpp"

namespace MESQUITE_NS {

class MsqError;
class PatchData;

/**\brief An interface for a mapping function of the form 
 * \f$\vec{x}(\vec{\xi})=\sum_{i=1}^n N_i(\vec{\xi})\vec{x_i}\f$,
 * where \f$\vec{x_i}\f$ is a point
 * in \f$\mathbf{R}^3\f$ (i.e. \f$x_i,y_i,z_i\f$),
 * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \end{array}\right\}\f$ 
 * for surface elements and 
 * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \zeta_i\\ \end{array}\right\}\f$ 
 * for volume elements. 
 *
 * This is an interface for describing a mapping function for a
 * single element topology.  A mapping function is assumed to be
 * of the following form: 
 * \f$\vec{x}(\vec{\xi})=\sum_{i=1}^n N_i(\vec{\xi})\vec{x_i}\f$
 * where \f$n\f$ is the number of nodes in the element,
 * \f$\vec{x_i}\f$ is a point
 * in \f$\mathbf{R}^3\f$ (i.e. \f$x_i,y_i,z_i\f$), and
 * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \end{array}\right\}\f$ 
 * for surface elements and 
 * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \zeta_i\\ \end{array}\right\}\f$ 
 * for volume elements.  For example,
 * for a linear quadrilateral element, the mapping function will be
 * of the form: 
 * \f$\vec{x}(\xi,\eta)=N_1(\xi,\eta)\vec{x_1}
 *                     +N_2(\xi,\eta)\vec{x_2}
 *                     +N_3(\xi,\eta)\vec{x_3}
 *                     +N_4(\xi,\eta)\vec{x_4}\f$
 *
 * A single implementation of this interface may support multiple 
 * types of elements of the same topology.  Element types within
 * a topology may vary by the presences or lack there of of mid-edge,
 * mid-face, and mid-element nodes.
 */
class MESQUITE_EXPORT MappingFunction 
{
public:

  virtual
  ~MappingFunction() {}

  /**\brief Get Mesquite::EntityTopology handled by this mapping function */
  virtual 
  EntityTopology element_topology() const = 0;
  
  /**\brief Get number of nodes in the element type 
   *
   * Get the number of nodes in the element type that the mapping
   * function implements.  It is assumed that the result of this
   * function, in combination with the element topology, is sufficient
   * to determine the element type.
   */
  virtual
  int num_nodes() const = 0;
  
  /**\brief Get sample points at which to evaluate mapping function
   *
   * Get the points within the element at which TMP quality metrics
   * that are a function of the mapping function Jacobian should
   * be evaluated.  The default (which may be overridden by individual
   * mapping functions) is to evaluate at all nodes.
   */
  virtual
  NodeSet sample_points( NodeSet higher_order_nodes ) const;

  /**\brief Mapping Function Coefficients
   *
   * This function returns the list of scalar values (\f$N_i\f$'s) resulting
   * from the evaluation of the mapping function coefficient terms
   * \f$N_1(\vec{\xi}), N_2(\vec{\xi}), \ldots, N_n(\vec{\xi})\f$
   * for a given \f$\vec{\xi}\f$.
   *\param location Where within the element at which to evaluate the coefficients.
   *\param nodeset  List of which nodes are present in the element.  
   *\param coefficients_out The coefficients (\f$N_i(\vec{\xi})\f$) for each 
   *                vertex in the element.
   *\param indices_out  The index ($i$ in $N_i$) for each term in 'coeffs_out'.
   *                The assumption is that mapping function implementations
   *                will not return zero coefficients.  This is not required,
   *                but for element types with large numbers of nodes it may
   *                have a significant impact on performance.
   */
  virtual 
  void coefficients( Sample location,
                     NodeSet nodeset,
                     double* coeff_out,
                     size_t* indices_out,
                     size_t& num_coeff_out,
                     MsqError& err ) const = 0;
                     
    /*\brief Convert connectivity list indices for different element types.
     *
     * Given two elements of the same topology but different types
     * (number of nodes) and a list of indices into the connectivity
     * list for one element type, convert the list to be indices
     * into a second element type such that the node in the same logical
     * position (e.g. middle of edge 1) is indicated.
     */
  static inline
  void convert_connectivity_indices( EntityTopology topology,
                                     int num_nodes_in_input_elem_type,
                                     int num_nodes_in_output_elem_type,
                                     size_t* index_list,
                                     unsigned num_indices,
                                     MsqError& err );
                                     
    /*\brief Convert connectivity list indices for different element types.
     *
     * Given an element type with the same topology as that of this
     * mapping function but with a different number of nodes, convert
     * indices into the connectivity list of this element type to
     * those of the specified element type such that indices indicate
     * nodes at the corresponding logical locations (e.g. middle of edge 1).
     */
  inline
  void convert_connectivity_indices( int num_nodes_in_output_element_type,
                                     size_t* index_list, 
                                     unsigned num_indices,
                                     MsqError& err ) const
    { convert_connectivity_indices( element_topology(), num_nodes(),
                                    num_nodes_in_output_element_type,
                                    index_list, num_indices, err ); }

private:
  static
  void convert_connectivity_indices_impl( EntityTopology topology,
                                     int num_nodes_in_input_elem_type,
                                     int num_nodes_in_output_elem_type,
                                     size_t* index_list,
                                     unsigned num_indices,
                                     MsqError& err );
};

/**\brief MappingFunction for topologically 2D (surface) elements. */
class MESQUITE_EXPORT MappingFunction2D : public MappingFunction
{
public:

  virtual
  ~MappingFunction2D() {}

  /**\brief Mapping Function Derivatives
   *
   * This function returns the partial derivatives of the mapping
   * function coefficient terms
   * \f$\nabla N_1(\vec{\xi}), \nabla N_2(\vec{\xi}), \ldots, \nabla N_n(\vec{\xi})\f$
   * evaluated for a given \f$\vec{\xi}\f$, where \f$\vec{x_i}\f$ is a point
   * in \f$\mathbf{R}^3\f$ (i.e. \f$x_i,y_i,z_i\f$).  
   * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \end{array}\right\}\f$ 
   * for surface elements and 
   * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \zeta_i\\ \end{array}\right\}\f$ 
   * for volume elements.
   *
   * The list of returned partial derivatives may be considered list of elements 
   * of a matrix \f$\mathbf{D}\f$ in row major order.  For surface elements,
   * \f$\mathbf{D}\f$ is a \f$n\times 2\f$ matrix and for volume elements it
   * is a \f$n \times 3\f$ matrix.  Each row of 
   * \f$\mathbf{D}\f$ corresponds to one of the
   * coefficient functions \f$N_i(\vec{\xi})\f$ and each column corresponds
   * to one of the components of \f$\vec{\xi}\f$ 
   * that the corresponding coefficient function is differentiated with
   * respect to. 
   *
   * \f$ \mathbf{D} = \left[ \begin{array}{ccc}
   *     \frac{\delta N_1}{\delta \xi} & \frac{\delta N_1}{\delta \eta} & \ldots \\
   *     \frac{\delta N_2}{\delta \xi} & \frac{\delta N_2}{\delta \eta} & \ldots \\
   *     \vdots & \vdots & \ddots \end{array} \right]\f$
   *
   * The Jacobian matrix (\f$\mathbf{J}\f$) of the mapping function can be calculated
   * as follows. Define a matrix \f$\mathbf{X}\f$ such that each column contains
   * the coordinates of the element nodes.
   *
   * \f$ \mathbf{X} = \left[ \begin{array}{ccc}
   *                   x_1 & x_2 & \ldots \\
   *                   y_1 & y_2 & \ldots \\
   *                   z_1 & z_2 & \ldots 
   *                  \end{array}\right]\f$
   *
   * The Jacobian matrix is then:
   *
   * \f$\mathbf{J} = \mathbf{X} \times \mathbf{D}\f$
   *
   * \f$\mathbf{X}\f$ is always \f$3\times n\f$, so \f$\mathbf{J}\f$ is
   * either \f$3\times 2\f$ (surface elements) or \f$3\times 3\f$ (volume
   * elements) depending on the dimensions of \f$\mathbf{D}\f$.
   *
   * If the Jacobian matrix of the mapping function is considered as a 
   * function of the element vertex coordinates \f$\mathbf{J}(\vec{x_1},\vec{x_2},\ldots)\f$ 
   * with \f$\vec{\xi}\f$ constant, then the gradient of that Jacobian matrix 
   * function (with respect
   * to the vertex coordinates) can be obtained from the same output list of
   * partial deravitves.
   *
   * \f$\frac{\delta \mathbf{J}}{\delta x_i} = 
   *         \left[ \begin{array}{ccc}
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots \\
   *         0 & 0 & \ldots \\ 
   *         0 & 0 & \ldots 
   *         \end{array} \right]\f$
   * \f$\frac{\delta \mathbf{J}}{\delta y_i} = 
   *         \left[ \begin{array}{ccc}
   *         0 & 0 & \ldots \\ 
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots \\
   *         0 & 0 & \ldots 
   *         \end{array} \right]\f$
   * \f$\frac{\delta \mathbf{J}}{\delta z_i} = 
   *         \left[ \begin{array}{ccc}
   *         0 & 0 & \ldots \\ 
   *         0 & 0 & \ldots \\
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots 
   *         \end{array} \right]\f$
   * 
   *
   *\param location Where within the element at which to evaluate the derivatives.
   *\param nodeset  List of which nodes are present in the element.  
   *\param vertices_out The list of vertices for which the corresponding
   *                coefficient in the mapping function is non-zero.  The
   *                vertices are specified by their index in the canonical
   *                ordering for an element with all mid-nodes present (i.e.
   *                first all the corner nodes, then the mid-edge nodes, ...).
   *\param d_coeff_d_xi_out The mapping function is composed of a series of 
   *                coefficient functions \f$N_i(\vec{\xi})\f$, one correspoding
   *                to the position \f$\vec{x_i}\f$ of each node in the
   *                element such that the mapping function is of the form:
   *                \f$\vec{x}(\vec{\xi})=\sum_{i=1}^n N_i(\vec{\xi})\vec{x_i}\f$.
   *                For each vertex indicated in vertex_indices_out, 
   *                this list contains the partial derivatives of the cooresponding
   *                coefficient function \f$N_i\f$ with respect to each 
   *                component of \f$\vec{\xi}\f$ in the same order as the
   *                corresponding nodes in vertex_indices_out. 
   *\param num_vtx  Output: The number of vertex indices and derivitive
   *                tuples returned in vertices_out and d_coeff_d_xi_out,
   *                respectively.
   */
  virtual 
  void derivatives( Sample location,
                    NodeSet nodeset,
                    size_t* vertex_indices_out,
                    MsqVector<2>* d_coeff_d_xi_out,
                    size_t& num_vtx,
                    MsqError& err ) const = 0;

  /**\brief Mapping function derivatives and Jacobian
   *
   * This function returns the partial derivatives of the mapping
   * function coefficient terms and the Jacobian calculated from
   * those terms and the cooresponding vertex coordinates.
   *
   * This function returns the same logical data as 'derivatives',
   * except that it also calculates the Jacobian from the actual
   * vertex coordinates.  Also, unlike the 'derivatives' function
   * which returns the vertex indices as positions in the element
   * connectivity list, this function is expected to 
   * a) return the actual indices of the vertices in the PatchData
   *    vertex list and
   * b) remove from the list of indices and derivatives and values
   *    corresponding to fixed vertices.
   *
   * The default implementation of this function will calculate the
   * Jacobian and modify the vertex and derivative lists returned
   * from "derivatives".  The default implementation serves as a
   * utility function for other classes using this one.  The function
   * is virtual to allow mapping function implementations to provide
   * an optimized version that avoids extra calculations for zero terms
   * in the derivative list.
   *
   *\param pd  The PatchData instance containing the vertex coordinates
   *           and element connectcivity.
   *\param element_number  The index of the mesh element in the PatchData.
   *\param nodeset         List of which nodes are present in the element.  
   *\param location Where within the element at which to evaluate the Jacobian.
   *\param vertex_patch_indices_out  For each free vertex in the element
   *                       the influences the mapping function value at
   *                       the specified logical location, the index of
   *                       that vertex in the PatchData.
   *\param d_coeff_d_xi_out For each vertex in 'vertex_patch_indices_out',
   *                       the partial derivatives of the corresponding
   *                       coefficient of the mapping function.
   *\param num_vtx_out     The number of values passed back in 
   *                       'vertex_patch_indices_out' and 'd_coeff_d_xi_out'.
   *\param jacobian_out    The Jacobian of the mapping function at the
   *                       specified logical location.
   */             
  virtual 
  void jacobian( const PatchData& pd,
                 size_t element_number,
                 NodeSet nodeset,
                 Sample location,
                 size_t* vertex_patch_indices_out,
                 MsqVector<2>* d_coeff_d_xi_out,
                 size_t& num_vtx_out,
                 MsqMatrix<3,2>& jacobian_out,
                 MsqError& err ) const;

  /**\brief Get ideal Jacobian matrix
   *
   * Returns the Jacobian matrix of an ideal element.  The orientation
   * of element or corresponding matrix is arbitrary.  The "ideal" element
   * should be scaled such the Jacobian (determinant of the Jacobian
   * matrix) is 1.0.
   *
   *\param location Where within the element at which to evaluate the Jacobian.
   *                Typically doesn't matter except for degenerate elements
   *                (e.g. pyramid as degenerate hex.)
   *\param jacobian_out    The Jacobian of the mapping function at the
   *                       specified logical location.
   */
  virtual
  void ideal( Sample location, 
              MsqMatrix<3,2>& jacobian_out,
              MsqError& err ) const;
};

/**\brief MappingFunction for topologically 3D (volume) elements. */
class MESQUITE_EXPORT MappingFunction3D : public MappingFunction
{
public:

  virtual
  ~MappingFunction3D() {}

  /**\brief Mapping Function Derivatives
   *
   * This group of methods return the partial derivatives of the mapping
   * function coefficient terms
   * \f$\nabla N_1(\vec{\xi}), \nabla N_2(\vec{\xi}), \ldots, \nabla N_n(\vec{\xi})\f$
   * evaluated for a given \f$\vec{\xi}\f$, where \f$\vec{x_i}\f$ is a point
   * in \f$\mathbf{R}^3\f$ (i.e. \f$x_i,y_i,z_i\f$).  
   * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \end{array}\right\}\f$ 
   * for surface elements and 
   * \f$\vec{\xi_i} = \left\{\begin{array}{c}\xi_i\\ \eta_i\\ \zeta_i\\ \end{array}\right\}\f$ 
   * for volume elements.
   *
   * The list of returned partial derivatives may be considered list of elements 
   * of a matrix \f$\mathbf{D}\f$ in row major order.  For surface elements,
   * \f$\mathbf{D}\f$ is a \f$n\times 2\f$ matrix and for volume elements it
   * is a \f$n \times 3\f$ matrix.  Each row of 
   * \f$\mathbf{D}\f$ corresponds to one of the
   * coefficient functions \f$N_i(\vec{\xi})\f$ and each column corresponds
   * to one of the components of \f$\vec{\xi}\f$ 
   * that the corresponding coefficient function is differentiated with
   * respect to. 
   *
   * \f$ \mathbf{D} = \left[ \begin{array}{ccc}
   *     \frac{\delta N_1}{\delta \xi} & \frac{\delta N_1}{\delta \eta} & \ldots \\
   *     \frac{\delta N_2}{\delta \xi} & \frac{\delta N_2}{\delta \eta} & \ldots \\
   *     \vdots & \vdots & \ddots \end{array} \right]\f$
   *
   * The Jacobian matrix (\f$\mathbf{J}\f$) of the mapping function can be calculated
   * as follows. Define a matrix \f$\mathbf{X}\f$ such that each column contains
   * the coordinates of the element nodes.
   *
   * \f$ \mathbf{X} = \left[ \begin{array}{ccc}
   *                   x_1 & x_2 & \ldots \\
   *                   y_1 & y_2 & \ldots \\
   *                   z_1 & z_2 & \ldots 
   *                  \end{array}\right]\f$
   *
   * The Jacobian matrix is then:
   *
   * \f$\mathbf{J} = \mathbf{X} \times \mathbf{D}\f$
   *
   * \f$\mathbf{X}\f$ is always \f$3\times n\f$, so \f$\mathbf{J}\f$ is
   * either \f$3\times 2\f$ (surface elements) or \f$3\times 3\f$ (volume
   * elements) depending on the dimensions of \f$\mathbf{D}\f$.
   *
   * If the Jacobian matrix of the mapping function is considered as a 
   * function of the element vertex coordinates \f$\mathbf{J}(\vec{x_1},\vec{x_2},\ldots)\f$ 
   * with \f$\vec{\xi}\f$ constant, then the gradient of that Jacobian matrix 
   * function (with respect
   * to the vertex coordinates) can be obtained from the same output list of
   * partial deravitves.
   *
   * \f$\frac{\delta \mathbf{J}}{\delta x_i} = 
   *         \left[ \begin{array}{ccc}
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots \\
   *         0 & 0 & \ldots \\ 
   *         0 & 0 & \ldots 
   *         \end{array} \right]\f$
   * \f$\frac{\delta \mathbf{J}}{\delta y_i} = 
   *         \left[ \begin{array}{ccc}
   *         0 & 0 & \ldots \\ 
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots \\
   *         0 & 0 & \ldots 
   *         \end{array} \right]\f$
   * \f$\frac{\delta \mathbf{J}}{\delta z_i} = 
   *         \left[ \begin{array}{ccc}
   *         0 & 0 & \ldots \\ 
   *         0 & 0 & \ldots \\
   *         \frac{\delta N_i}{\delta \xi} & \frac{\delta N_i}{\delta \eta} & \ldots 
   *         \end{array} \right]\f$
   * 
   *
   *\param location Where within the element at which to evaluate the derivatives.
   *\param nodeset  List of which nodes are present in the element.  
   *\param vertices_out The list of vertices for which the corresponding
   *                coefficient in the mapping function is non-zero.  The
   *                vertices are specified by their index in the canonical
   *                ordering for an element with all mid-nodes present (i.e.
   *                first all the corner nodes, then the mid-edge nodes, ...).
   *\param d_coeff_d_xi_out The mapping function is composed of a series of 
   *                coefficient functions \f$N_i(\vec{\xi})\f$, one correspoding
   *                to the position \f$\vec{x_i}\f$ of each node in the
   *                element such that the mapping function is of the form:
   *                \f$\vec{x}(\vec{\xi})=\sum_{i=1}^n N_i(\vec{\xi})\vec{x_i}\f$.
   *                For each vertex indicated in vertex_indices_out, 
   *                this list contains the partial derivatives of the cooresponding
   *                coefficient function \f$N_i\f$ with respect to each 
   *                component of \f$\vec{\xi}\f$ in the same order as the
   *                corresponding nodes in vertex_indices_out. 
   *\param num_vtx  Output: The number of vertex indices and derivitive
   *                tuples returned in vertices_out and d_coeff_d_xi_out,
   *                respectively.
   */
  virtual 
  void derivatives( Sample location,
                    NodeSet nodeset,
                    size_t* vertex_indices_out,
                    MsqVector<3>* d_coeff_d_xi_out,
                    size_t& num_vtx,
                    MsqError& err ) const = 0;
                    
                    
  /**\brief Mapping function derivatives and Jacobian
   *
   * This function returns the partial derivatives of the mapping
   * function coefficient terms and the Jacobian calculated from
   * those terms and the cooresponding vertex coordinates.
   *
   * This function returns the same logical data as 'derivatives',
   * except that it also calculates the Jacobian from the actual
   * vertex coordinates.  Also, unlike the 'derivatives' function
   * which returns the vertex indices as positions in the element
   * connectivity list, this function is expected to 
   * a) return the actual indices of the vertices in the PatchData
   *    vertex list and
   * b) remove from the list of indices and derivatives and values
   *    corresponding to fixed vertices.
   *
   * The default implementation of this function will calculate the
   * Jacobian and modify the vertex and derivative lists returned
   * from "derivatives".  The default implementation serves as a
   * utility function for other classes using this one.  The function
   * is virtual to allow mapping function implementations to provide
   * an optimized version that avoids extra calculations for zero terms
   * in the derivative list.
   *\param pd  The PatchData instance containing the vertex coordinates
   *           and element connectcivity.
   *\param element_number  The index of the mesh element in the PatchData.
   *\param nodeset         List of which nodes are present in the element.  
   *\param location Where within the element at which to evaluate the Jacobian.
   *\param vertex_patch_indices_out  For each free vertex in the element
   *                       the influences the mapping function value at
   *                       the specified logical location, the index of
   *                       that vertex in the PatchData.
   *\param d_coeff_d_xi_out For each vertex in 'vertex_patch_indices_out',
   *                       the partial derivatives of the corresponding
   *                       coefficient of the mapping function.
   *\param num_vtx_out     The number of values passed back in 
   *                       'vertex_patch_indices_out' and 'd_coeff_d_xi_out'.
   *\param jacobian_out    The Jacobian of the mapping function at the
   *                       specified logical location.
   */             
  virtual 
  void jacobian( const PatchData& pd,
                 size_t element_number,
                 NodeSet nodeset,
                 Sample location,
                 size_t* vertex_patch_indices_out,
                 MsqVector<3>* d_coeff_d_xi_out,
                 size_t& num_vtx_out,
                 MsqMatrix<3,3>& jacobian_out,
                 MsqError& err ) const;

  /**\brief Get ideal Jacobian matrix
   *
   * Returns the Jacobian matrix of an ideal element.  The orientation
   * of element or corresponding matrix is arbitrary.  The "ideal" element
   * should be scaled such the Jacobian (determinant of the Jacobian
   * matrix) is 1.0.
   *
   *\param location Where within the element at which to evaluate the Jacobian.
   *                Typically doesn't matter except for degenerate elements
   *                (e.g. pyramid as degenerate hex.)
   *\param jacobian_out    The Jacobian of the mapping function at the
   *                       specified logical location.
   */
  virtual
  void ideal( Sample location, 
              MsqMatrix<3,3>& jacobian_out,
              MsqError& err ) const;
};


inline void
MappingFunction::convert_connectivity_indices( EntityTopology topo,
                                               int input_type,
                                               int output_type,
                                               size_t* index_list,
                                               unsigned num_indices,
                                               MsqError& err )
{
    // If the types are the same or either type has only corner
    // vertices, then no conversion is necessary.
  const int num_corners = TopologyInfo::corners(topo);
  if (input_type != output_type && input_type != num_corners && output_type != num_corners)
    convert_connectivity_indices_impl( topo, input_type, output_type, index_list, num_indices, err );
} 
  
                                                  

} // namespace Mesquite

#endif
