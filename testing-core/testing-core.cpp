#include "testing-core.h"

namespace Hermes
{
  namespace Testing
  {
#ifdef _WINDOWS
    long get_current_virtual_memory()
    {
      PROCESS_MEMORY_COUNTERS pmc;
      GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
      return pmc.WorkingSetSize;
    }

    long get_peak_virtual_memory()
    {
      PROCESS_MEMORY_COUNTERS pmc;
      GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
      return pmc.PeakWorkingSetSize;
    }
#else
    static int parseLine(char* line)
    {
      int i = strlen(line);
      while (*line < '0' || *line > '9') line++;
      line[i-3] = '\0';
      i = atoi(line);
      return i;
    }

    long get_current_virtual_memory()
    { //Note: this value is in KB!
      FILE* file = fopen("/proc/self/status", "r");
      int result = -1;
      char line[128];


      while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
          result = parseLine(line);
          break;
        }
      }
      fclose(file);
      return result;
    }
#endif

    bool check_expected_memory(long expected_memory)
    {
      long actual_memory = get_current_virtual_memory();
#ifndef _WINDOWS
      // Linux function is in kB.
      actual_memory /= 1024;
#endif
      if(actual_memory > expected_memory)
      {
        double difference = (actual_memory - expected_memory) / 1048576.;
        printf("Memory difference = %f MB!\n", difference);
        return false;
      }

      return true;
    }

    bool test_value(double obtained_value, double expected_value, const char* identifier, double absolute_precision)
    {
      if(std::abs(expected_value - obtained_value) < absolute_precision) 
        return true;
      else
      {
        std::cout << "Failed test: " << identifier << std::endl;
        std::cout << "Difference: " << expected_value << " - " << obtained_value << " higher than the precision (" << absolute_precision << ")." << std::endl;
        return false;
      }
    }

    bool compare_files(const char* filename_1, const char* filename_2)
    {
      FILE* f1 = fopen(filename_1, "r");
      FILE* f2 = fopen(filename_2, "r");

      int N = 10000;
      char* buf1 = new char[N];
      char* buf2 = new char[N];

      do
      {
        size_t r1 = fread(buf1, 1, N, f1);
        size_t r2 = fread(buf2, 1, N, f2);

        if (r1 != r2 || memcmp(buf1, buf2, r1))
        {
          fclose(f1);
          fclose(f2);
          std::cout << "Failed test: files not the same: " << filename_1 << ", " << filename_2 << std::endl;
          return false;
        }
      }
      while (!feof(f1) && !feof(f2));

      fclose(f1);
      fclose(f2);
      return true;
    }
  }
}