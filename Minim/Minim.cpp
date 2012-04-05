void Minim::reduce()
{
  // Calculate the values at the simplices.
  for (size_t idx = 0; idx < d_simplex.cols(); ++idx)
  {
    d_engine.reset();
    size_t required = 1000;
    do
    {
      d_engine.calculate(d_simplex.col(idx), required);
      double diff = kolmogorov_difference(d_data, d_engine.results());
      // NOTE How do I calculate this value? Do I just want a relative difference in the data? How to define this?
      
    }
    while (d_values[idx])
  }
}