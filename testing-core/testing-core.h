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
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

typedef ::std::complex<double> complex;

#define SHOW_OUTPUT

namespace Hermes
{
  namespace Testing
  {
    extern Hermes::Mixins::Loggable logger;

    long get_current_virtual_memory();

#ifdef _WINDOWS
    long get_peak_virtual_memory();
#endif

    bool check_expected_memory(long expected_memory);

    bool test_value(double obtained_value, double expected_value, const char* identifier, double absolute_precision = 1e-4);

    bool compare_files(const char* filename_1, const char* filename_2);
  }
}

using namespace Hermes::Testing;
