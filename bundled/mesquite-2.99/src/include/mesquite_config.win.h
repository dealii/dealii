

#if _MSC_VER <= 1200
#  define MSQ_USE_OLD_C_HEADERS
#  define MSQ_FUNCTION ""
#  include <float.h>
#else
#  define MSQ_FUNCTION __FUNCDNAME__
#  include <cfloat>
#endif

inline int finite( double d ) { return _finite(d); }

