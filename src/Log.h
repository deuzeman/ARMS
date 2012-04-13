#pragma once

#include <fstream>
#include <iostream>

class Log: public std::ofstream
{
  static Log *s_instance;

  Log(char const *filename);
  
  public:
    static bool ionode;
    
    ~Log();
    static void open(char const *filename, int rank = 0);
    static void shut();
    
    static std::ofstream &put();
};

inline Log::~Log()
{}

inline Log::Log(char const *filename)
  : std::ofstream(filename, std::ofstream::app)
{}

inline std::ofstream &log()
{
  return Log::put();
}