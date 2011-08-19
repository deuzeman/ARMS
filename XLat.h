#pragma once

#include <sstream>

class XLat
{
  std::istringstream d_data;

  public:
    XLat(char const *input)
      : d_data(input)
    {};

    template< typename Type >
    operator Type()
    {
      Type result;
      d_data >> result;
      return result;
    }
};
