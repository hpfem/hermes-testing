#include "testing-core.h"

namespace Hermes
{
  namespace Testing
  {
    Hermes::Mixins::Loggable logger(true);

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


      while (fgets(line, 128, file) != nullptr){
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
        std::cout.precision(15);
        std::cout << "Failed test: " << identifier << std::endl;
        std::cout << "Difference: " << expected_value << " - " << obtained_value << " higher than the precision (" << absolute_precision << ")." << std::endl;
        return false;
      }
    }

    static bool compare_files_internal(const char* filename_1, const char* filename_2)
    {
        std::ifstream lFile(filename_1, std::ios::in | std::ios::binary);
        std::ifstream rFile(filename_2, std::ios::in | std::ios::binary);


        if(!lFile.good() || !rFile.good())
        {
            return false;
        }

        std::streamsize lReadBytesCount = 0;
        std::streamsize rReadBytesCount = 0;

        int bufferSize = 100;
        char p_lBuffer[bufferSize];
        char p_rBuffer[bufferSize];

        do {
            lFile.read(p_lBuffer, bufferSize);
            rFile.read(p_rBuffer, bufferSize);
            lReadBytesCount = lFile.gcount();
            rReadBytesCount = rFile.gcount();

            if (lReadBytesCount != rReadBytesCount || std::memcmp(p_lBuffer, p_rBuffer, lReadBytesCount) != 0)
            {
                return false;
            }
        } while (lFile.good() || rFile.good());

        return true;
    }

    bool compare_files(const char* filename_1, const char* filename_2)
    {
      if(compare_files_internal(filename_1, filename_2))
        return true;
      else
      {
          std::cout << "Failed test: files not the same: " << filename_1 << ", " << filename_2 << std::endl;
          return false;
        }
    }
  }
}
