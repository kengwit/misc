#include "auto_handle.h"

int 
main (void)
{
  int i = 5;
  auto_handle<int *> ih = &i;
  return 0;
}
