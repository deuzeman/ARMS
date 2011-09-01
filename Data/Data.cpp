#include <algorithm>

#include <Data.h>

Data::Data(char const *filename)
  : d_cur(NONE)
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