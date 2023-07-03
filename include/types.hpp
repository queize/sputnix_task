#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/tools/norms.hpp>

const double mu = 3.98e14;
const double h = 1;
const size_t T = 86400;
const double R = 6373e3;

template <typename T>
using Vector = boost::numeric::ublas::vector<T>;

template <typename T>
using Matrix = boost::numeric::ublas::matrix<T>;

namespace ublas = boost::numeric::ublas;
