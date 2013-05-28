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
  }
}