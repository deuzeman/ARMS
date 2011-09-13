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
      N = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nu = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nu = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("m = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      m = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("a6 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a6 = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("a7 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a7 = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("a8 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      a8 = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("scale = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      scale = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nEig_min = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nEig_min = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nEig_max = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nEig_max = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nDet = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nDet = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("iter = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      iter = XLat(sline.c_str() + idx + 2);
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

    if (sline.find("data = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      data = std::string(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("N = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      N = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("nu = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      nu = XLat(sline.c_str() + idx + 2);
      continue;
    }

    if (sline.find("m = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      m[0] = istr;
      m[1] = istr;
      m[2] = istr;
      continue;
    }

    if (sline.find("a6 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a6[0] = istr;
      a6[1] = istr;
      a6[2] = istr;
      continue;
    }

    if (sline.find("a7 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a7[0] = istr;
      a7[1] = istr;
      a7[2] = istr;
      continue;
    }

    if (sline.find("a8 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      a8[0] = istr;
      a8[1] = istr;
      a8[2] = istr;
      continue;
    }

    if (sline.find("scale = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      scale[0] = istr;
      scale[1] = istr;
      scale[2] = istr;
      continue;
    }

    if (sline.find("iter = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      XLat istr(sline.c_str() + idx + 2);
      iter[0] = istr;
      iter[1] = istr;
      continue;
    }

    if (sline.find("tol = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      tol = XLat(sline.c_str() + idx + 2);
      continue;
    }
  }
}
