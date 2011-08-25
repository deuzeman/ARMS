#include <fstream>
#include <iostream>

#include <Params.h>

void Params::parseInput(char const *filename)
{
  std::ifstream input(filename);
  char line[256];

  while(input.getline(line, 256))
  {
    std::string sline(line);

    if (sline.find("N = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      N = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nu = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nu = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("m = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      m = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a6 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a6 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a7 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a7 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a8 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a8 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nEigs = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nEigs = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nDet = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nDet = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("iter = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      iter = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("output = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      output = (sline.c_str() + idx + 2);
      continue;
    }
  }
}

void FitParams::parseInput(char const *filename)
{
  std::ifstream input(filename);
  char line[256];

  while(input.getline(line, 256))
  {
    std::string sline(line);

    if (sline.find("N = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      N = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nu = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nu = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("m = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      m = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a6 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a6 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a7 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a7 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a8 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a8 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nEigs = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nEigs = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nDet = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nDet = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("iter = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      iter = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("output = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      output = (sline.c_str() + idx + 2);
      continue;
    }
  }
}
