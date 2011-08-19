double Minim::chiSq(Eigen::VectorXd const &pars)
{
    static Eigen::ArrayXd eval;

    d_generator.changeParams(params(pars));
    d_generator.calculate(iter);

    size_t resCtr = 0;

    for (size_t num = 0; num < nEigs - 1; ++num)
        for (size_t den = num; den < nEigs; ++den)
            eval.coeff(resCtr++) = d_generator.avRatio(num, den);

    return ((eval - ratData) / sdRatData).squaredNorm();
}

double Minim::brent(size_t const idx)
{
    // Performs Brent's line search in the direction d_dirs(idx) with start d_pars(idx)
    // It is effectively one-dimensional, but we fit a scaling factor alpha.
    // Starting position is alpha = 0, by construction.
    // Limits are imposed from the overal constraints on the parameters -- we're not
    // allowed to use values bringing any of the components out of its Cartesian box.
    double const cgold = 0.3819660;
    double const zeps = 2E-8;

    // Determine brackets on the scaling factor...
    Eigen::MatrixXd limits(dim, 2);
    limits.row(0) = (brackets.row(0) - d_pars.col(idx)) / d_dirs.col(idx);
    limits.row(1) = (brackets.row(1) - d_pars.col(idx)) / d_dirs.col(idx); // Can be better, perhaps.

    double min = limits.coeff(0, limits.row(0).minCoeff());
    double max = limits.coeff(1, limits.row(1).maxCoeff());
    double mid = 0.0;

    double d = 0.0;
    double w = mid;
    double v = mid;

    double fMid = chiSq(d_pars.col(idx));
    double fv = fMid;
    double fw = fMid;

    for (int brIter = 0; brIter < maxBrIter; ++brIter)
    {
        double boundMid =   0.5 * (max - min);
        double tol1 = tol * std::abs(mid) + zeps;
        double tol2 = 2.0 * tol1;

        if (std::abs(mid - boundMid) <= (tol2 - 0.5 * (max - min)))
        {
//             d_fmin = fMid; // NOTE How to return this?
//             d_xmin = mid;
            return mid;
        }

        if (std::abs(e) > tol1)
        {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);

            if (q > 0.0)
                p = -p;
            q = std::abs(q);
            double e_temp = e;
            e = d;
            if (abs(p) >= abs(0.5 * q * e_temp) || p <= q * (a - x) || p >= q * (b - x))
            {
                e = (x >= xm) ? a - x : b - x;
                d = cgold * x;
                d = cgold * e;
            }
            else
            {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = (xm >= x) ? tol : -tol;
            }
        }
        else
        {
            e = (x >= xm) ? a - x : b - x);
            d = cgold * e;
        }

        u = (abs(d) >= tol1) ? x + d : x + ((d >= 0) ? tol1, -tol1);
        fu = func(u);

        if (fu <= fx)
        {
            if (u >= x)
                a = x;
            else
                b = x;
            v  =  w;
            w  =  x;
            x  =  u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            if (u < x)
                a = u;
            else
                b = u;
            if (fu <= fw || w == x)
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else
            {
                if (fu <= fv || v == x || v == w)
                {
                    v = u;
                    fv = fu;
                }
            }
        }
    }
    d_error = 1;
}

