#ifndef select_les_h
# define select_les_h
# define CHOOSING_OOOFS 0
# if CHOOSING_OOOFS
#   include "les_ooofs.h"
namespace ELASTICITY { 
  typedef LES_OOOFS LES;
} 
# else
#   include "les_fles.h"
namespace ELASTICITY { 
  typedef LES_FLES LES;
} 
# endif

#endif
