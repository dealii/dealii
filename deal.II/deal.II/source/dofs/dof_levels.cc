

#include <base/exceptions.h>
#include <base/memory_consumption.h>
#include <dofs/dof_levels.h>


namespace internal
{
  namespace DoFHandler
  {

#if deal_II_dimension == 1

    unsigned int
    DoFLevel<1>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (lines));
    }

#endif

#if deal_II_dimension == 2
    
    unsigned int
    DoFLevel<2>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (quads));
    }

#endif

#if deal_II_dimension == 3

    unsigned int
    DoFLevel<3>::memory_consumption () const
    {
      return (DoFLevel<0>::memory_consumption () +
              MemoryConsumption::memory_consumption (hexes));
    }

#endif
    
  }
  
}
