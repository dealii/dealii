#
# Check for various compiler bugs:
#


#
# Versions of GCC before 3.0 had a problem with the explicit
# instantiation of member templates when the member was in fact
# an operator. In that case, they needed the "template" keyword,
# which is actually not allowed at this place. Test case is
#
# /* ----------------------------------------------- */
# struct X
# {
#     template <typename T2>
#     X operator = (T2 &) { return X(); };
# };
#
# template X X::operator=<float> (float &);
# /* ---------------------------------------------------------- */
#
# The compiler only groks this if the "operator=" is prepended
# by "template". We detect this, and either set the
# DEAL_II_MEMBER_OP_TEMPLATE_INST to "template" or nothing, so
# that it gets expanded to the right string needed in this place.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  struct X
  {
      template <typename T2>
      X operator = (T2 &) { return X(); }
  };

  template X X::operator=<float> (float &);
  int main(){return 0;}
  "
  DEAL_II_MEMBER_OP_TEMPLATE_INST_OK)

IF(DEAL_II_MEMBER_OP_TEMPLATE_INST_OK)
  SET(DEAL_II_MEMBER_OP_TEMPLATE_INST "")
ELSE()
  SET(DEAL_II_MEMBER_OP_TEMPLATE_INST "template")
ENDIF()



#
# Some compiler versions, notably ICC, have trouble with the
# following code in which we explicitly call a destructor.
# This has to be worked around with a typedef. The problem is
# that the workaround fails with some other compilers, so that
# we can not unconditionally use the workaround...
#
CHECK_CXX_COMPILER_BUG(
  "
  namespace dealii
  {
    namespace FEValuesViews
    {
      template <int dim, int spacedim> struct Scalar {};
    }

    template <int dim, int spacedim>
    struct X
    {
        FEValuesViews::Scalar<dim,spacedim> scalars[dim*spacedim];

        void f()
          {
            scalars[0].dealii::FEValuesViews::Scalar<dim,spacedim>::~Scalar ();
          }
    };

    template struct X<2,2>;
  }
  int main() { return 0; }
  "
  DEAL_II_EXPLICIT_DESTRUCTOR_BUG)
