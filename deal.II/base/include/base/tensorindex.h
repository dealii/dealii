/*----------------------------   tensorindex.h     ---------------------------*/
/*      $Id$ */
#ifndef __tensorindex_H
#define __tensorindex_H
/*----------------------------   tensorindex.h     ---------------------------*/

#include <base/exceptions.h>


/**
 * Rank-independent access to elements of #Tensor#.  A little class
 * template remembering #rank# integers.
 *
 * This template is implemented only for ranks between 1 and 4. If you
 * really need higher ranks, please modify the source code accordingly
 * or contact the developer.
 *
 * @author Guido Kanschat, 1999
 */
template<int rank>
class TensorIndex
{
  private:
				     /**
				      * Field of indices.
				      */
    unsigned int index[rank];
  public:
				     /**
				      * Constructor taking #rank# indices. 
				      */
    TensorIndex(...);
    
				     /**
				      * Access operator returning index
				      * in #n#th component
				      */
    unsigned int operator () (const unsigned int n) const;
    
				     /**
				      * Exception.
				      */
    DeclException1(ExcRank, int,
		   << "Index " << arg1 << " higher than maximum " << rank-1);  
};



template<>
class TensorIndex<4>
{
  private:
				     /**
				      * Field of indices.
				      */
    unsigned int index[4];
  public:
				     /**
				      * Constructor taking #rank# indices.
				      */
    TensorIndex (const unsigned int i0,
		 const unsigned int i1,
		 const unsigned int i2,
		 const unsigned int i3)
      {
	index[0] = i0;
	index[1] = i1;
	index[2] = i2;
	index[3] = i3;
      }
  
  
				     /**
				      * Access operator returning index
				      * in #n#th component
				      */
    unsigned int operator () (const unsigned int n) const
      {
	Assert(n<4, ExcRank(n));
	return index[n];

      }

				     /**
				      * Exception
				      */
    DeclException1(ExcRank, unsigned int,
		   << "Index " << arg1 << " higher than maximum 3");  
};



template<>
class TensorIndex<3>
{
  private:
				     /**
				      * Field of indices.
				      */
    unsigned int index[3];
  public:
				     /**
				      * Constructor taking #rank# indices.
				      */
    TensorIndex(const unsigned int i0,
		const unsigned int i1,
		const unsigned int i2)
      {
	index[0] = i0;
	index[1] = i1;
	index[2] = i2;
      }
  
  
				     /**
				      * Access operator returning index
				      * in #n#th component
				      *
				      */
    unsigned int operator () (const unsigned int n) const
      {
	Assert(n<3, ExcRank(n));
	return index[n];

      }
  
				     /**
				      * Exception
				      */
    DeclException1(ExcRank, unsigned int,
		   << "Index " << arg1 << " higher than maximum 2");  
};



template<>
class TensorIndex<2>
{
  private:
				     /**
				      * Field of indices.
				      */
    unsigned int index[4];
  public:
				     /**
				      * Constructor taking #rank# indices.
				      */
    TensorIndex(const unsigned int i0,
		const unsigned int i1)
      {
	index[0] = i0;
	index[1] = i1;
      }
  
  
				     /**
				      * Access operator returning index
				      * in #n#th component
				      */
    unsigned int operator () (const unsigned int n) const
      {
	Assert(n<2, ExcRank(n));
	return index[n];

      }

    				     /**
				      * Exception
				      */
    DeclException1(ExcRank, unsigned int,
		   << "Index " << arg1 << " higher than maximum 1");  
};



template<>
class TensorIndex<1>
{
  private:
				     /**
				      * Field of indices.
				      */
    unsigned int index[1];
  public:
				     /**
				      * Constructor taking #rank# indices.
				      */
    TensorIndex(const unsigned int i0)
      {
	index[0] = i0;
      }
  
  
				     /**
				      * Access operator returning index
				      * in #n#th component
				      */
    unsigned int operator () (const unsigned int n) const
      {
	Assert(n<1, ExcRank(n));
	return index[n];

      }

				     /**
				      * Exception
				      */
    DeclException1(ExcRank, unsigned int,
		   << "Index " << arg1 << " higher than maximum 0");  
};



template<int rank>
inline unsigned int
TensorIndex<rank>::
operator() (unsigned int n) const
{
  Assert(n<rank, ExcRank(n));
  return index[n];
}



/*----------------------------   tensorindex.h     ---------------------------*/
/* end of #ifndef __tensorindex_H */
#endif
/*----------------------------   tensorindex.h     ---------------------------*/
