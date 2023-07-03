#include "types.hpp"
#include "geocoef.hpp"
#include "orientation.hpp"


size_t factorial(size_t n)
{
    return std::tgamma(n + 1);
}

double normalize_geofoef(size_t n, size_t m, double C_nm)
{
    return std::sqrt(factorial(n + m) / (2 * factorial(n - m))) * C_nm / std::sqrt(2*n + 1); 
}

double legandrePnm(size_t n, size_t m, double z)
{
    if (m == 0 && n < 2)
    {
        if (n == 0)
        {
            return 1;
        }
        else if (n == 1)
        {
            return z;
        }
    }
    else if (m == n)
    {
        return (2 * n - 1) * std::sqrt(1 - std::pow(z, 2)) * legandrePnm(n - 1, n - 1, z);
    }
    if (n - 1 < m)
    {
        return 0;
    }
    else if (n - 2 < m)
    {
        return (2 * n - 1) * z * legandrePnm(n - 1, m, z) / (n - m);
    }
    else
    {
        return ((2 * n - 1) * z * legandrePnm(n - 1, m, z) - (n - 1 + m) * legandrePnm(n - 2, m, z)) / (n - m);
    }
}

Vector<double> func(double t, const Vector<double> &X)
{
    namespace tools = boost::math::tools;
    Orientation ori(t);
    size_t n = 10;

    Vector<double> out(6);
    Vector<double> r = ublas::subrange(X, 0, 3);
    Vector<double> v = ublas::subrange(X, 3, 6);

    ublas::subrange(out, 0, 3) = v;
    r = ori.GCRS_to_ITRS(r);
    Vector<double> a_geo(3);
    double lambda = 0, phi = 0;
    for(size_t i = 0; i <= n; ++i)
    {
        for(size_t j = 0; j <= i; ++j)
        {
            double C_ij = normalize_geofoef(i, j, ggm03c_terms[i*(i+1)/2 + j][0]);
            double S_ij = normalize_geofoef(i, j, ggm03c_terms[i*(i+1)/2 + j][1]);
            a_geo += (i + 1) * std::pow(R/tools::l2_norm(r), i + 3) * legandrePnm(i, j, std::sin(phi)) * (C_ij*std::cos(j*lambda) + S_ij*std::sin(j*lambda))*r;
        }
    }
    a_geo *= -mu / std::pow(R, 3);
    a_geo = ori.ITRS_to_GCRS(a_geo);
    ublas::subrange(out, 3, 6) = a_geo;
    return out;
}

void Rk4(
    double &t,
    Vector<double> &X,
    std::function<Vector<double>(double, const Vector<double> &)> func)
{
    Vector<double> k1, k2, k3, k4;
    k1 = func(t, X);
    k2 = func(t + h / 2, X + k1 * h / 2);
    k3 = func(t + h / 2, X + k2 * h / 2);
    k4 = func(t + h / 2, X + k3 * h);
    X += (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
    t += h;
}



int main(int argc, char const *argv[])
{
    std::ofstream out("output.txt");

    Vector<double> P(3), V(3);
    double t = 0;
    P(0) = 1702631.521;
    P(1) = 126415.744;
    P(2) = 6769207.534;
    V(0) = -5734.531;
    V(1) = -4667.074;
    V(2) = 1528.123;
    Vector<double> X(6);
    ublas::subrange(X, 0, 3) = P;
    ublas::subrange(X, 3, 6) = V;
    out << P << '\t' << V << '\n';
    for (size_t i = 0; i < size_t(T/h); ++i)
    {
        Rk4(t, X, func);
        out << ublas::subrange(X, 0, 3) << '\t'
            << ublas::subrange(X, 3, 6) << '\n';
    }
    return 0;
}
