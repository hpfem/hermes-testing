#include "hermes2d.h"

#ifdef _WINDOWS
#include "windows.h"
#include "psapi.h"
#else
  #include "stdlib.h"
  #include "stdio.h"
  #include "string.h"
#endif

using namespace Hermes;
using namespace Hermes::Hermes2D;


namespace Hermes
{
  namespace Testing
  {
    long get_current_virtual_memory();

#ifdef _WINDOWS
    long get_peak_virtual_memory();
#endif

    bool check_expected_memory(long expected_memory);
  }
}