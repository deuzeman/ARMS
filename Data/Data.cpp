#include <algorithm>
#include <iostream>
#include <mpi.h>
#include <cstdio>
#include <Data.h>
#include <sys/stat.h>

Data::Data(char const *filename, int bootSeed, bool const parallel)
  : d_cur(NONE), d_bootSeed(bootSeed)
{
  if (parallel)
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
  }
  else
  {
    // We assume a dataframe structure, with some meta-info prefixed.
    std::ifstream datastream(filename);

    char *buffer = new char[8192];

    datastream.getline(buffer, 8192);
    std::istringstream meta(buffer + 1); // Skip comment #
    meta >> d_nSamp >> d_minEv >> d_maxEv;

    d_data.resize(d_nSamp, d_maxEv - d_minEv); // Assuming minEv < 0 and the 0 is not included!

    for (size_t row = 0; row < d_nSamp; ++row)
      for (size_t col = 0; col < (d_maxEv - d_minEv); ++col)
        datastream >> d_data(row, col);
    delete[] buffer;
  }

  d_normalization = 1.0 / (d_data.cols() * (0.5 - (Eigen::ArrayXd::LinSpaced(d_data.rows(), 1.0, static_cast< double >(d_data.rows())) / static_cast< double >(d_data.rows()))).square().sum());
  
  if (!bootSeed)
    return;
  
  // Add a random sampling here, to enable bootstrapping fits
  CRandomSFMT sampler(bootSeed);
  Eigen::ArrayXXd bootSamp(d_data.rows(), d_data.cols());
  for (size_t ctr = 0; ctr < d_data.rows(); ++ctr)
    bootSamp.row(ctr) = d_data.row(sampler.IRandomX(0, d_data.rows() - 1));
  
  // We want to treat the resampled dataset as the true one
  // So what we want to do, is swapping the data in bootSamp with the original d_data,
  // then forget about d_data all together
  d_data.swap(bootSamp);
  
}

void Data::bootstrap(size_t const nBoot, Eigen::ArrayXXd const &data) const
{
  // Bootstrap the data
  Eigen::ArrayXXd bootSamp(data.rows(), data.cols());
  Eigen::ArrayXXd bootHist(nBoot, data.cols());
  CRandomSFMT sampler(time(0));

  for (size_t bootCtr = 0; bootCtr < nBoot; ++bootCtr)
  {
    for (size_t sampCtr = 0; sampCtr < data.rows(); ++sampCtr)
      bootSamp.row(sampCtr) = data.row(sampler.IRandomX(0, data.rows() - 1));
    bootHist.row(bootCtr) = (bootSamp.colwise().mean());
  }

  // Calculate the standard deviation from this.
  d_summary.row(1) = bootHist.square().colwise().mean();
  d_summary.row(1) -= bootHist.colwise().mean().square();
  d_summary.row(1) = std::sqrt(d_summary.row(1));
}

Eigen::ArrayXXd const &Data::average(size_t const nBoot) const
{
  if (d_cur != AVE)
  {
    d_summary.resize(2, d_data.cols());
    d_summary.row(0) = d_data.colwise().mean();

    bootstrap(nBoot, d_data);

    d_cur = AVE;
  }
  return d_summary;
}

Eigen::ArrayXXd const &Data::ratios(size_t const nBoot) const
{
  if (d_cur != RAT)
  {
    size_t ratDim = (d_data.cols() * (d_data.cols() - 1)) / 2;

    d_summary.resize(2, ratDim);
    Eigen::ArrayXXd ratios(d_data.rows(), ratDim);

    size_t ctr = 0;
    for (size_t num = 0; num < (d_data.cols() - 1); ++num)
      for (size_t den = num + 1; den < d_data.cols(); ++den)
        ratios.col(ctr++) =  d_data.col(num) / d_data.col(den);

    d_summary.row(0) = ratios.colwise().mean();

    bootstrap(nBoot, ratios);

    d_cur = RAT;
  }
  return d_summary;
}

Eigen::ArrayXXd const &Data::cumulant() const
{
  if (d_cur != CUM)
  {
    d_summary.resizeLike(d_data);
    Eigen::ArrayXd tempVec(d_data.rows()); // We need contiguous data, so a temporary vector

    for (size_t col = 0; col < d_data.cols(); ++col)
    {
      tempVec = d_data.col(col);
      std::sort(&tempVec[0], &tempVec[d_data.rows() - 1] + 1);
      d_summary.col(col) = tempVec;
    }

    d_cur = CUM;
  }
  return d_summary;
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
