#include <Log.h>

#include <cstdlib>

Log *Log::s_instance = 0;

void Log::open(char const *filename)
{
  if (Log::s_instance)
  {
    std::cerr << "[ERROR] Double opening of log!" << std::endl;
    exit(1);
  }
  Log::s_instance = new Log(filename);
}

void Log::shut()
{
  if (s_instance)
    s_instance->close();
  delete Log::s_instance;
  Log::s_instance = 0;
}

std::ofstream &Log::put()
{
  if (!Log::s_instance)
  {
    std::cerr << "[ERROR] Attempted to use log without opening!" << std::endl;
    exit(1);
  }
  return *s_instance;
}