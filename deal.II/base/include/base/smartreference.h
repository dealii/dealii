/*----------------------------   smartreference.h     ---------------------------*/
/*      $Id$                 */
#ifndef __smartreference_H
#define __smartreference_H
/*----------------------------   smartreference.h     ---------------------------*/

#include <base/subscriptor.h>

/**
 * Smart references avoid destruction of a referenced object.  This
 * class has not been fully developed, since the compiler could not
 * resolve the dot operator in a convenient manner. The use of
 * #SmartPointer# is recommended, instead.
 */
template<typename T>
class SmartReference
{
    T& t;

  public:
				     /**
				      * Constructor taking a normal reference.
				      */
    SmartReference(const T& tt)
		    : t(tt) 
      {
	t.subscribe();
      }
    
				     /**
				      * Destructor, removing the subscription.
				      */
    ~SmartReference()
      {
	t.unsubscribe();
      }
  
				     /**
				      * Conversion to normal reference
				      */
    operator T& () const
      {
	return t;
      }
};

/*----------------------------   smartreference.h     ---------------------------*/
/* end of #ifndef __smartreference_H */
#endif
/*----------------------------   smartreference.h     ---------------------------*/
