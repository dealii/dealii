/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2010 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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

    (2010) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file CompareQM.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "CompareQM.hpp"
#include "MsqError.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

namespace MESQUITE_NS {

CompareQM::CompareQM( QualityMetric* primary, QualityMetric* other, 
                      const char* primary_name, const char* other_name )
  : primaryMetric( primary ), 
    otherMetric( other ),
    abortOnMismatch( false ),
    toleranceFactor( 1e-6 )
{
  if (primary_name)
    primaryName = primary_name;
  if (other_name)
    otherName = other_name;
}

CompareQM::~CompareQM()
{}

void CompareQM::abort_on_mismatch( double tolerance ) 
{
  abortOnMismatch = true;
  toleranceFactor = tolerance;
}

void CompareQM::do_not_abort()
{ 
  abortOnMismatch = false;
}

bool CompareQM::will_abort_on_mismatch() const
{
  return abortOnMismatch;
}

QualityMetric::MetricType CompareQM::get_metric_type() const 
{
  MetricType type1, type2;
  type1 = primaryMetric->get_metric_type();
  type2 = otherMetric->get_metric_type();
  if (type1 != type2) {
    std::cerr << "Incompatible metric types in CompareQM" << std::endl
              << __FILE__ << ':' << __LINE__ << std::endl;
    if (abortOnMismatch)
      abort();
    // otherwise just return some type because this function can't
    // flag an error.  The mismatch should cause get_evaluations
    // to fail anyway.
  }
  return type1;  
}

std::string CompareQM::get_name() const 
{
  return primaryMetric->get_name();
}

int CompareQM::get_negate_flag() const 
{
  return primaryMetric->get_negate_flag();
}

void CompareQM::get_evaluations( PatchData& pd,
                                 std::vector<size_t>& handles,
                                 bool free_only,
                                 MsqError& err )
{
  primaryMetric->get_evaluations( pd, handles, free_only, err ); MSQ_ERRRTN(err);
  
  std::vector<size_t> handles2;
  otherMetric->get_evaluations( pd, handles2, free_only, err ); MSQ_ERRRTN(err);
  
  std::vector<size_t> handles1( handles );
  std::sort( handles1.begin(), handles1.end() );
  std::sort( handles2.begin(), handles2.end() );
  if (handles1 != handles2) {
    MSQ_SETERR(err)("Incompatible metrics cannot be compared", MsqError::INVALID_STATE);
  }
}

bool CompareQM::check_valid( size_t handle, bool rval1, bool rval2 ) 
{
  if (rval1 != rval2) {
    std::cerr << "Metrics returned conflicting validity at location " << std::hex << handle << std::dec
              << std::endl << __FILE__ << ":" << __LINE__ << std::endl;
    if (abortOnMismatch) 
      abort();
    else
      return false;
  }
  
  return rval1;
}

double CompareQM::epsilon( double a, double b )
{ 
  return toleranceFactor * std::max( 1.0, std::max( a, b ) );
}

void CompareQM::check_value( size_t handle, double value, double value2 )
{
  valPrimary.add_value( value );
  valOther.add_value( value2 );
  valDiff.add_value( fabs(value - value2) );
  
  if (abortOnMismatch) {
    if (fabs(value - value2) > epsilon(value, value2)) {
      std::cerr << "Metric values to not match at location " << std::hex << handle << std::dec << std::endl
                << "Primary: " << value << std::endl
                << "Other  : " << value2 << std::endl
                << __FILE__ << ":" << __LINE__ << std::endl;
      abort();
    }
  }
}

bool CompareQM::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  double value2;
  bool rval1, rval2;
  rval1 = primaryMetric->evaluate( pd, handle, value, err ); MSQ_ERRZERO(err);
  rval2 = otherMetric->evaluate( pd, handle, value2, err ); MSQ_ERRZERO(err);
  if (!check_valid( handle, rval1, rval2 ))
    return false;
  
  check_value( handle, value, value2 );
  return true;
}

void CompareQM::index_mismatch( size_t handle,
                                const std::vector<size_t>& idx1,
                                const std::vector<size_t>& idx2,
                                MsqError& err )
{
  std::vector<size_t>::const_iterator i;

  std::cerr << "Metrics cannot be compared at location " << std::hex << handle << std::dec << 
                  " because they are incompatible." << std::endl
                  << "Primary metric vertices: ";
  if (idx1.empty())
    std::cerr << "(empty)";
  else {
    i = idx1.begin();
    std::cerr << *i;
    for (++i; i != idx1.end(); ++i) 
      std::cerr << ',' << *i;
  }
  std::cerr << std::endl << "Other metric vertices: " ;
  if (idx2.empty())
    std::cerr << "(empty)";
  else {
    i = idx2.begin();
    std::cerr << *i;
    for (++i; i != idx2.end(); ++i) 
      std::cerr << ',' << *i;
  }
  std::cerr << std::endl;

  if (abortOnMismatch)
    abort();
  MSQ_SETERR(err)("Cannot compare incompatible metrics", MsqError::INVALID_STATE);
}
                               

void CompareQM::check_indices( size_t handle,
                               const std::vector<size_t>& idx1,
                               const std::vector<size_t>& idx2,
                               std::vector<size_t>& map_out,
                               MsqError& err )
{
  if (idx1.size() != idx2.size()) {
    index_mismatch( handle, idx1, idx2, err );
    MSQ_ERRRTN(err);
  }
  
  std::vector<size_t>::const_iterator i, j;
  map_out.clear();
  for (i = idx1.begin(); i != idx1.end(); ++i) {
    j = std::find( idx2.begin(), idx2.end(), *i );
    if (j == idx2.end()) {
      index_mismatch( handle, idx1, idx2, err );
      MSQ_ERRRTN(err);
    }
    map_out.push_back( j - idx2.begin() );
  }
}

bool CompareQM::evaluate_with_indices( PatchData& pd,
                                       size_t handle,
                                       double& value,
                                       std::vector<size_t>& indices,
                                       MsqError& err )
{
  bool valid1, valid2;
  double value2;
  std::vector<size_t> indices2, junk;
  valid1 = primaryMetric->evaluate_with_indices( pd, handle, value, indices, err ); MSQ_ERRZERO(err);
  valid2 = otherMetric->evaluate_with_indices( pd, handle, value2, indices2, err ); MSQ_ERRZERO(err);
  if (!check_valid( handle, valid1, valid2 ))
    return false;
  
  check_value( handle, value, value2 );
  check_indices( handle, indices, indices2, junk, err ); MSQ_ERRZERO(err);
  return true;
}

void CompareQM::check_grad( size_t handle,
                            const std::vector<size_t>& indices,
                            const std::vector<size_t>& index_map,
                            const std::vector<Vector3D>& grad1,
                            const std::vector<Vector3D>& grad2 )
{
  assert(index_map.size() == indices.size());
  assert(index_map.size() == grad1.size());
  assert(index_map.size() == grad2.size());
  for (size_t i = 0; i < index_map.size(); ++i) {
    gradPrimary.add( grad1[i] );
    gradOther.add( grad2[i] );
    gradDiff.add_diff( grad1[i], grad2[index_map[i]] );
    
    if (abortOnMismatch) {
      if ((grad1[i] - grad2[index_map[i]]).length() > epsilon(grad1[i].length(), grad2[index_map[i]].length()) ) {
        std::cerr << "Gradients differ for metric evaluation at " << std::hex << handle << std::dec << std::endl
                  << "Primary metric derivs with respect to vertex " << indices[i] << ": " << grad1[i] << std::endl
                  << "Other metric derivs with presect to vertex " << indices[i] << ": " << grad2[index_map[i]] << std::endl
                  << __FILE__ << ":" << __LINE__ << std::endl;
        abort();
      }
    }
  }
}

bool CompareQM::evaluate_with_gradient( PatchData& pd,
                                        size_t handle,
                                        double& value,
                                        std::vector<size_t>& indices,
                                        std::vector<Vector3D>& grad,
                                        MsqError& err )
{
  bool valid1, valid2;
  double value2;
  std::vector<size_t> indices2, map;
  std::vector<Vector3D> grad2;
  valid1 = primaryMetric->evaluate_with_gradient( pd, handle, value, indices, grad, err ); MSQ_ERRZERO(err);
  valid2 = otherMetric->evaluate_with_gradient( pd, handle, value2, indices2, grad2, err ); MSQ_ERRZERO(err);
  if (!check_valid( handle, valid1, valid2 ))
    return false;
  
  check_value( handle, value, value2 );
  check_indices( handle, indices, indices2, map, err ); MSQ_ERRZERO(err);
  check_grad( handle, indices, map, grad, grad2 );
  return true;
}

void CompareQM::check_hess_diag( size_t handle,
                                 const std::vector<size_t>& indices,
                                 const std::vector<size_t>& index_map,
                                 const std::vector<SymMatrix3D>& hess1,
                                 const std::vector<SymMatrix3D>& hess2 )
{
  assert(index_map.size() == indices.size());
  assert(index_map.size() == hess1.size());
  assert(index_map.size() == hess2.size());
  for (size_t i = 0; i < index_map.size(); ++i) {
    hessPrimary.add_diag( hess1[i] );
    hessOther.add_diag( hess2[i] );
    hessDiff.add_diag_diff( hess1[i], hess2[index_map[i]] );
    
    if (abortOnMismatch) {
      double eps = epsilon( Frobenius(hess1[i]), Frobenius(hess2[index_map[i]]) );
      if (Frobenius(hess1[i] - hess2[index_map[i]]) > eps) {
        std::cerr << "Hessian diagonal blocks differ for metric evaluation at " << std::hex << handle << std::dec << std::endl
                  << "For second derivatives with repsect to vertex " << indices[i] << "twice" << std::endl
                  << "Primary metric derivs: " << hess1[i] << std::endl
                  << "Other metric derivs: " << hess2[index_map[i]] << std::endl
                  << __FILE__ << ":" << __LINE__ << std::endl;
        abort();
      }
    }
  }
}

bool CompareQM::evaluate_with_Hessian_diagonal( PatchData& pd,
                                        size_t handle,
                                        double& value,
                                        std::vector<size_t>& indices,
                                        std::vector<Vector3D>& grad,
                                        std::vector<SymMatrix3D>& hess,
                                        MsqError& err )
{
  bool valid1, valid2;
  double value2;
  std::vector<size_t> indices2, map;
  std::vector<Vector3D> grad2;
  std::vector<SymMatrix3D> hess2;
  valid1 = primaryMetric->evaluate_with_Hessian_diagonal( pd, handle, value, indices, grad, hess, err ); MSQ_ERRZERO(err);
  valid2 = otherMetric->evaluate_with_Hessian_diagonal( pd, handle, value2, indices2, grad2, hess2, err ); MSQ_ERRZERO(err);
  if (!check_valid( handle, valid1, valid2 ))
    return false;
  
  check_value( handle, value, value2 );
  check_indices( handle, indices, indices2, map, err ); MSQ_ERRZERO(err);
  check_grad( handle, indices, map, grad, grad2 );
  check_hess_diag( handle, indices, map, hess, hess2 );
  return true;
}

void CompareQM::check_hess( size_t handle,
                            const std::vector<size_t>& indices,
                            const std::vector<size_t>& index_map,
                            const std::vector<Matrix3D>& hess1,
                            const std::vector<Matrix3D>& hess2 )
{
  const size_t n = index_map.size();
  const size_t N = (n + 1) * n / 2;
  assert(n == indices.size());
  assert(N == hess1.size());
  assert(N == hess2.size());
  
  for (size_t r = 0; r < n; ++r) {
    const size_t r2 = index_map[r];
    for (size_t c = r; c < n; ++c) {
      const size_t c2 = index_map[c];
      size_t idx1 = n*r - r*(r+1)/2 + c;
      Matrix3D h2;
      if (r2 <= c2) {
        size_t idx2 = n*r2 - r2*(r2+1)/2 + c2;
        h2 = hess2[idx2];
      }
      else {
        size_t idx2 = n*c2 - c2*(c2+1)/2 + r2;
        h2 = transpose(hess2[idx2]);
      }
      
      if (r == c) {
        hessPrimary.add_diag( hess1[idx1] );
        hessOther.add_diag( h2 );
        hessDiff.add_diag_diff( hess1[idx1], h2 );
      }
      else {
        hessPrimary.add_nondiag( hess1[idx1] );
        hessOther.add_nondiag( h2 );
        hessDiff.add_nondiag_diff( hess1[idx1], h2 );
      }
      
      if (abortOnMismatch) {
        double eps = epsilon( sqrt(Frobenius_2(hess1[idx1])), sqrt(Frobenius_2(h2)) );
        if (sqrt(Frobenius_2(hess1[idx1] - h2)) > eps) {
          std::cerr << "Hessian blocks differ for metric evaluation at " << std::hex << handle << std::dec << std::endl
                    << "For second derivatives with repsect to vertices " << indices[r] << " and " << indices[c] << std::endl
                    << "Primary metric derivs: " << hess1[idx1] << std::endl
                    << "Other metric derivs: " << h2 << std::endl
                    << __FILE__ << ":" << __LINE__ << std::endl;
          abort();
        }
      }
    }
  }
}

bool CompareQM::evaluate_with_Hessian( PatchData& pd,
                                       size_t handle,
                                       double& value,
                                       std::vector<size_t>& indices,
                                       std::vector<Vector3D>& grad,
                                       std::vector<Matrix3D>& hess,
                                       MsqError& err )
{
  bool valid1, valid2;
  double value2;
  std::vector<size_t> indices2, map;
  std::vector<Vector3D> grad2;
  std::vector<Matrix3D> hess2;
  valid1 = primaryMetric->evaluate_with_Hessian( pd, handle, value, indices, grad, hess, err ); MSQ_ERRZERO(err);
  valid2 = otherMetric->evaluate_with_Hessian( pd, handle, value2, indices2, grad2, hess2, err ); MSQ_ERRZERO(err);
  if (!check_valid( handle, valid1, valid2 ))
    return false;
  
  check_value( handle, value, value2 );
  check_indices( handle, indices, indices2, map, err ); MSQ_ERRZERO(err);
  check_grad( handle, indices, map, grad, grad2 );
  check_hess( handle, indices, map, hess, hess2 );
  return true;
}

void CompareQM::GradStat::add( Vector3D grad ) 
{
  x.add_value(grad[0]);
  y.add_value(grad[1]);
  z.add_value(grad[2]);
}

void CompareQM::GradStat::add_diff( Vector3D grad1, Vector3D grad2 ) 
{
  x.add_value(fabs(grad1[0]-grad2[0]));
  y.add_value(fabs(grad1[1]-grad2[1]));
  z.add_value(fabs(grad1[2]-grad2[2]));
}

void CompareQM::HessStat::add_diag( Matrix3D hess ) 
{
  xx.add_value(  hess[0][0] );
  xy.add_value( (hess[0][1]+hess[1][0])/2 );
  xz.add_value( (hess[0][2]+hess[2][0])/2 );
  yy.add_value(  hess[1][1] );
  yz.add_value( (hess[1][2]+hess[2][1])/2 );
  zz.add_value(  hess[2][2] );
}

void CompareQM::HessStat::add_diag( SymMatrix3D hess ) 
{
  xx.add_value( hess[0] );
  xy.add_value( hess[1] );
  xz.add_value( hess[2] );
  yy.add_value( hess[3] );
  yz.add_value( hess[4] );
  zz.add_value( hess[5] );
}

void CompareQM::HessStat::add_diag_diff( Matrix3D hess1, Matrix3D hess2 ) 
{
  Matrix3D d = hess1 - hess2;
  Matrix3D hess( fabs(d[0][0]), fabs(d[0][1]), fabs(d[0][2]),
                 fabs(d[1][0]), fabs(d[1][1]), fabs(d[1][2]),
                 fabs(d[2][0]), fabs(d[2][1]), fabs(d[2][2]) );
  add_diag( hess );
}

void CompareQM::HessStat::add_diag_diff( SymMatrix3D hess1, SymMatrix3D hess2 ) 
{
  SymMatrix3D d = hess1 - hess2;
  SymMatrix3D hess( fabs(d[0]), fabs(d[1]), fabs(d[2]),
                                fabs(d[3]), fabs(d[4]),
                                            fabs(d[5]) );
  add_diag( hess );
}

void CompareQM::HessStat::add_nondiag( Matrix3D hess ) 
{
  xx.add_value( hess[0][0] );
  xy.add_value( hess[0][1] );
  xy.add_value( hess[1][0] );
  xz.add_value( hess[0][2] );
  xz.add_value( hess[2][0] );
  yy.add_value( hess[1][1] );
  yz.add_value( hess[1][2] );
  yz.add_value( hess[2][1] );
  zz.add_value( hess[2][2] );
}

void CompareQM::HessStat::add_nondiag_diff( Matrix3D hess1, Matrix3D hess2 ) 
{
  Matrix3D d = hess1 - hess2;
  Matrix3D hess( fabs(d[0][0]), fabs(d[0][1]), fabs(d[0][2]),
                 fabs(d[1][0]), fabs(d[1][1]), fabs(d[1][2]),
                 fabs(d[2][0]), fabs(d[2][1]), fabs(d[2][2]) );
  add_nondiag( hess );
}

static void print( const char* title,
                   const char* name1,
                   const char* name2,
                   const SimpleStats& s1,
                   const SimpleStats& s2,
                   const SimpleStats& sd )
{
  const char named[] = "difference";
  int len = std::max( std::max( strlen(named), strlen(title) ), std::max( strlen(name1), strlen(name2) ) );
  std::vector<char> dashes(len, '-');
  dashes.push_back('\0');
  printf( "%*s %12s %12s %12s %12s %12s\n", len, title, 
     "minimum", "average", "rms", "maximum", "std.dev." );
  printf( "%s ------------ ------------ ------------ ------------ ------------\n", &dashes[0] );
  printf( "%*s % 12g % 12g % 12g % 12g % 12g\n", len, name1, 
                                     s1.minimum(),
                                     s1.average(),
                                     s1.rms(),
                                     s1.maximum(),
                                     s1.standard_deviation() );
  printf( "%*s % 12g % 12g % 12g % 12g % 12g\n", len, name2, 
                                     s2.minimum(),
                                     s2.average(),
                                     s2.rms(),
                                     s2.maximum(),
                                     s2.standard_deviation() );
  printf( "%*s % 12g % 12g % 12g % 12g % 12g\n", len, named, 
                                     sd.minimum(),
                                     sd.average(),
                                     sd.rms(),
                                     sd.maximum(),
                                     sd.standard_deviation() );
  printf("\n");
}

void CompareQM::print_stats() const
{
  std::string name1 = primaryName.empty() ? primaryMetric->get_name() : primaryName;
  std::string name2 =   otherName.empty() ?   otherMetric->get_name() :   otherName;
  print( "Values", name1.c_str(), name2.c_str(), valPrimary, valOther, valDiff );
  print( "Gradient X", name1.c_str(), name2.c_str(), gradPrimary.x, gradOther.x, gradDiff.x );
  print( "Gradient Y", name1.c_str(), name2.c_str(), gradPrimary.y, gradOther.y, gradDiff.y );
  print( "Gradient Z", name1.c_str(), name2.c_str(), gradPrimary.z, gradOther.z, gradDiff.z );
  print( "Hessian XX", name1.c_str(), name2.c_str(), hessPrimary.xx, hessOther.xx, hessDiff.xx );
  print( "Hessian XY", name1.c_str(), name2.c_str(), hessPrimary.xy, hessOther.xy, hessDiff.xy );
  print( "Hessian XZ", name1.c_str(), name2.c_str(), hessPrimary.xz, hessOther.xz, hessDiff.xz );
  print( "Hessian YY", name1.c_str(), name2.c_str(), hessPrimary.yy, hessOther.yy, hessDiff.yy );
  print( "Hessian YZ", name1.c_str(), name2.c_str(), hessPrimary.yz, hessOther.yz, hessDiff.yz );
  print( "Hessian ZZ", name1.c_str(), name2.c_str(), hessPrimary.zz, hessOther.zz, hessDiff.zz );
}

} // namespace MESQUITE_NS
