/*----------------------------   data_out.h     ---------------------------*/
/*      $Id$                 */
#ifndef __data_out_H
#define __data_out_H
/*----------------------------   data_out.h     ---------------------------*/


#include <basic/data_out_base.h>


class DataOut1 : protected DataOutBase
{
  public:
  protected:
    virtual void make_patch_list () const;
};



/*----------------------------   data_out.h     ---------------------------*/
/* end of #ifndef __data_out_H */
#endif
/*----------------------------   data_out.h     ---------------------------*/
