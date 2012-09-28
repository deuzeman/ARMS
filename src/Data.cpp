#include <algorithm>
#include <iostream>
#include <mpi.h>
#include <cstdio>
#include <Data.h>
#include <sys/stat.h>

Data::Data(char const *filename)
{
  int mpirank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  // Broken LUSTRE forces us to use POSIX, then broadcast
  unsigned long int filesize = 0;
  char *buffer;

  // Detect if the data file even exists (otherwise things crash in a nasty way)
  size_t error = 0;
  if (mpirank == 0)
  {
    struct stat sbuf;
    if (stat(filename, &sbuf))
    {
      std::cerr << "Data file appears to be non-existent!" << std::endl;
      error = 1;
    }
  }
  MPI_Bcast(&error, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  if (error)
  {
    MPI_Finalize();
    exit(1);
  }

  if (mpirank == 0)
  {
    std::ifstream instream(filename);
    instream.seekg(0, std::ios::end);
    filesize = instream.tellg();
    instream.seekg(0, std::ios::beg);
    MPI_Bcast(&filesize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    buffer = new char[filesize + 1]; // To avoid weirdness, assign a little extra memory
    instream.read(buffer, filesize);
    instream.close();
    buffer[filesize] = '\0';
  }
  else
  {
    MPI_Bcast(&filesize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    buffer = new char[filesize + 1]; // To avoid weirdness, assign a little extra memory
  }
  MPI_Bcast(buffer, filesize, MPI_CHAR, 0, MPI_COMM_WORLD);

  // Now buffer should be equal for everybody -- we can continue as before
  std::istringstream input(buffer + 1); // Skip comment #
  input >> d_nSamp >> d_minEv >> d_maxEv;

  d_data.resize(d_nSamp, d_maxEv - d_minEv); // Assuming minEv < 0 and the 0 is not included!

  for (size_t row = 0; row < d_nSamp; ++row)
    for (size_t col = 0; col < (d_maxEv - d_minEv); ++col)
      input >> d_data(row, col);
  delete[] buffer;
  
  d_average.resize(d_data.cols());
  d_error.resize(d_data.cols());

  d_average = d_data.colwise().mean();
  d_error = d_data.square().colwise().mean();
  d_error -= d_average.square();
  d_error = d_error.sqrt();

  std::cout << "[DEBUG] d_average: " << d_average << std::endl;
  std::cout << "[DEBUG] d_error: " << d_error << std::endl;

  d_cumulant.resizeLike(d_data);
  Eigen::ArrayXd tempVec(d_data.rows()); // We need contiguous data, so a temporary vector
  
  for (size_t col = 0; col < d_data.cols(); ++col)
  {
    tempVec = d_data.col(col);
    std::sort(&tempVec[0], &tempVec[d_data.rows() - 1] + 1);
    d_cumulant.col(col) = tempVec;
  }
}

double *Data::average() const
{
  double *res = new double[d_data.cols()];
  for (size_t idx = 0; idx < d_data.cols(); ++idx)
    res[idx] = d_average[idx];
  return res;
}

double *Data::error() const
{
  double *res = new double[d_data.cols()];
  for (size_t idx = 0; idx < d_data.cols(); ++idx)
    res[idx] = d_error[idx];
  return res;
}

double *Data::flatPerColumn() const
{
  double *pdata = new double[d_data.rows() * d_data.cols()];
  for (size_t col = 0; col < d_data.cols(); ++col)
  {
    for (size_t row = 0; row < d_data.rows(); ++row)
      pdata[row + col * d_data.rows()] = d_data.coeff(row, col);
    std::sort(pdata + col * d_data.rows(), pdata + (col + 1) * d_data.rows());
  }
  return pdata;
}

double *Data::flat() const
{
  double *pdata = new double[d_data.rows() * d_data.cols()];
  for (size_t row = 0; row < d_data.rows(); ++row)
    for (size_t col = 0; col < d_data.cols(); ++col)
      pdata[row * d_data.cols() + col] = d_data.coeff(row, col);
  std::sort(pdata, pdata + d_data.rows() * d_data.cols());
  return pdata;
}

double *Data::flat(size_t const &col) const
{
  double *pdata = new double[d_data.rows()];
  for (size_t row = 0; row < d_data.rows(); ++row)
    pdata[row] = d_data.coeff(row, col);
  std::sort(pdata, pdata + d_data.rows());
  return pdata;
}
