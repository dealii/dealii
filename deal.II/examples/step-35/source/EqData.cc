/*
    Implementation of the classes that represent the exact solution
    and external force FOR THE PROBLEM ON THE DISK


    There's not much to say here. All these classes just represent
    mathematical formulae

    by Abner Salgado.
*/



#include "../include/EqData.h"



// Force function methods
template<int dim> Force<dim>::Force( const double initial_time ): Multi_Component_Function<dim>( initial_time ){
}



template<int dim> void  Force<dim>::value_list( const std::vector<Point<dim> > &points, std::vector<double> &values,
                                                  const unsigned int ) const{
  const unsigned int n_points = points.size();
  Assert( values.size() == n_points, ExcDimensionMismatch( values.size(), n_points ) );
  for (unsigned int i=0; i<n_points; ++i)
    values[i] = Force<dim>::value( points[i] );
}



template<int dim> inline double Force<dim>::value(const Point<dim> &p, const unsigned int) const{
  double t = FunctionTime::get_time(),
         cosx = std::cos( p(0) ),
         sinx = std::sin( p(0) ),
         cosy = std::cos( p(1) ),
         siny = std::sin( p(1) ),
         cost = std::cos(t),
         sint = std::sin(t),
         return_value = 0.;
  switch( Multi_Component_Function<dim>::component ){
    case 0:
    // y*sin(t)-x*cos(t)^2+cos(x)*sin(y)*sin(t)
      return_value = p(1)*sint - p(0)*cost*cost + cosx*siny*sint;
      break;
    case 1:
    // -x*sin(t)-y*cos(t)^2+sin(x)*cos(y)*sin(t)
      return_value = -p(0)*sint - p(1)*cost*cost + sinx*cosy*sint ;

      break;
    default:
      Assert( false, ExcNotImplemented() );
  };
  return return_value;
}



// Velocity function methods
template<int dim> Velocity<dim>::Velocity(const double initial_time): Multi_Component_Function<dim>( initial_time ){
}



template<int dim> void  Velocity<dim>::value_list( const std::vector<Point<dim> > &points, std::vector<double> &values,
                                                  const unsigned int ) const{
  const unsigned int n_points = points.size();
  Assert( values.size() == n_points, ExcDimensionMismatch( values.size(), n_points ) );
  for (unsigned int i=0; i<n_points; ++i)
    values[i] = Velocity<dim>::value( points[i] );
}



template<int dim> inline double Velocity<dim>::value(const Point<dim> &p, const unsigned int) const{
  double return_value = std::cos( Function<dim>::get_time() );
  switch( Multi_Component_Function<dim>::component ){
    case 0:
     // -y*cos(t)
       return_value *= -p(1);
      break;
    case 1:
      // x*cos(t)
      return_value *= p(0);
      break;
    default:
      Assert( false, ExcNotImplemented() );
  };
  return return_value;
}



template<int dim> inline Tensor<1,dim> Velocity<dim>::gradient(const Point<dim> &p, const unsigned int) const{
  Tensor<1,dim> return_value;
  switch( Multi_Component_Function<dim>::component ){
    // [0, -cos(t)]
    case 0:
     return_value[0] = 0.;
     return_value[1] = -std::cos( Function<dim>::get_time() );
      break;
    case 1:
    // [cos(t), 0]
     return_value[0] = std::cos( Function<dim>::get_time() );
     return_value[1] = 0.;
      break;
    default:
      Assert( false, ExcNotImplemented() );
  };
  return return_value;
}



template<int dim> void Velocity<dim>::gradient_list( const std::vector<Point<dim> > &points, std::vector< Tensor<1,dim> > &gradients,
                                                        const unsigned int ) const{
  const unsigned int n_points = points.size();
  Assert( gradients.size() == n_points, ExcDimensionMismatch( gradients.size(), n_points ) );
  for (unsigned int i=0; i<n_points; ++i)
    gradients[i] = Velocity<dim>::gradient( points[i] );
}



// Pressure function methods
template<int dim> Pressure<dim>::Pressure(const double initial_time): Function<dim>(1,initial_time){}



template<int dim> inline double Pressure<dim>::value(const Point<dim> &p, const unsigned int) const{
  // sin(x-y+t)
  return std::sin( p(0) )*std::sin( p(1) )*std::sin( Function<dim>::get_time() );
}



template<int dim> inline Tensor<1,dim> Pressure<dim>::gradient(const Point<dim> &p, const unsigned int) const{
  // [cos(x)*sin(y)*sin(t), sin(x)*cos(y)*sin(t)]
  return Point<dim>( std::cos( p(0) )*std::sin( p(1) )*std::sin( Function<dim>::get_time() ),
                       std::sin( p(0) )*std::cos( p(1) )*std::sin( Function<dim>::get_time() ) );
}



template<int dim> void Pressure<dim>::value_list( const std::vector<Point<dim> > &points, std::vector<double> &values,
                                                    const unsigned int ) const{
  const unsigned int n_points = points.size();
  Assert( values.size() == n_points, ExcDimensionMismatch( values.size(), n_points ) );
  for (unsigned int i=0; i<n_points; ++i)
    values[i] = Pressure<dim>::value( points[i] );
}



template<int dim> inline void Pressure<dim>::gradient_list( const std::vector<Point<dim> > &points, std::vector< Tensor<1,dim> > &gradients,
                                                                const unsigned int ) const{
  const unsigned int n_points = points.size();
  Assert( gradients.size() == n_points, ExcDimensionMismatch( gradients.size(), n_points ) );
  for (unsigned int i=0; i<n_points; ++i)
    gradients[i] = Pressure<dim>::gradient( points[i] );
}



// explicit template instantiation
template class Force<deal_II_dimension>;
template class Velocity<deal_II_dimension>;
template class Pressure<deal_II_dimension>;
