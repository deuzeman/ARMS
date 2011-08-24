#include <fstream>
#include <iostream>

#include <Params.h>

Params *parseInput(char const *filename)
{
  std::ifstream input(filename);
  char line[256];
  Params *result = new Params();

  while(input.getline(line, 256))
  {
    std::string sline(line);

    if (sline.find("N = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->N = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nu = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->nu = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("m = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->m = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a6 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->a6 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a7 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->a7 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("a8 = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->a8 = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nEigs = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->nEigs = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("nDet = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->nDet = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("iter = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->iter = XLat(sline.c_str() + idx + 1);
      continue;
    }

    if (sline.find("output = ") != sline.npos)
    {
      size_t idx = sline.find("=");
      result->output = (sline.c_str() + idx + 2);
      continue;
    }
  }

  return result;
}
