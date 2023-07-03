#pragma once

#include "types.hpp"

class Orientation
{
    double t;
    double JD(double t);
    Matrix<double> orientaion = Matrix<double>{3, 3};
    Matrix<double> create_Pt();
    Matrix<double> create_St();
    Matrix<double> create_Nt();
    Matrix<double> create_Pr();
    Matrix<double> create_BT();

public:
    Orientation(double curr_time);
    Vector<double> GCRS_to_ITRS(const Vector<double> &vec);
    Vector<double> ITRS_to_GCRS(const Vector<double> &vec);
};