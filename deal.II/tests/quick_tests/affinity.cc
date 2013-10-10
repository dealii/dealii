#include <deal.II/grid/tria.h>
#include <stdio.h>
#include <sched.h>

int main ()
{
  // we need this, otherwise gcc will not link against deal.II
  dealii::Triangulation<2> test;

  cpu_set_t my_set;
  CPU_ZERO(&my_set);

  unsigned int len = sizeof(my_set);
  int   ret = sched_getaffinity(0, len, &my_set);

  if (ret!=0)
    {
      printf("sched_getaffinity() failed, return value: %d\n", ret);
      return -1;
    }

  unsigned int bits_set = 0;//not supported on old kernels: CPU_COUNT(&my_set);
  for (int i=0;i<CPU_SETSIZE;++i)
    bits_set += CPU_ISSET(i,&my_set);

  if (bits_set==1)
    {
      printf("Warning: sched_getaffinity() returns that we can only use one CPU.\n");
      return 1;
    }
  printf("ncpus=%d, mask=%08X\n", bits_set, *(unsigned int*)(&my_set));
}
