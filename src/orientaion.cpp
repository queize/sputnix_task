#include "orientation.hpp"
double Orientation::JD(double t)
{
    return 2451545.1;
}
Matrix<double> Orientation::create_Pt()
{
    double xp = 0, yp = 0;
    Matrix<double> m(3, 3);
    m(0, 0) =   1; m(0, 1) = 0;     m(0, 2) = xp;
    m(1, 0) =   0; m(1, 1) = 1;     m(1, 2) = -yp;
    m(2, 0) = -xp; m(2, 1) = yp;    m(2, 2) = 1;
    return m;
}
Matrix<double> Orientation::create_St()
{
    double S = 0;
    Matrix<double> m(3, 3);
    m(0, 0) =  std::cos(S);    m(0, 1) = std::sin(S);  m(0, 2) = 0;
    m(1, 0) = -std::sin(S);    m(1, 1) = std::cos(S);  m(1, 2) = 0;
    m(2, 0) =            0;    m(2, 1) =           0;  m(2, 2) = 1;
    return m;
}
Matrix<double> Orientation::create_Nt()
{
    double epsilon_0 = 23*M_PI/180;
    double delta_epsilon = 0, d_epsilon = 0;
    double delta_psi = 0, d_psi = 0;
    double A = -epsilon_0 - delta_epsilon - d_epsilon;
    double B = -delta_psi - d_psi;
    double C =  epsilon_0;
    Matrix<double> m(3, 3);

    m(0, 0) =  std::cos(B);
    m(0, 1) =  std::sin(B)*std::cos(C);
    m(0, 2) = -std::sin(B)*std::cos(C);

    m(1, 0) = -std::sin(B)*std::cos(A);
    m(1, 1) =  std::cos(A)*std::cos(B)*std::cos(C) - std::sin(A)*std::sin(C);
    m(1, 2) =  std::sin(A)*std::sin(B);

    m(2, 0) =  std::sin(A)*std::sin(B);
    m(2, 1) = -std::sin(A)*std::cos(B)*std::cos(C) - std::cos(A)*std::sin(C);
    m(2, 2) = -std::sin(A)*std::cos(B)*std::sin(C) + std::cos(A)*std::cos(C);
    return m;
}
Matrix<double> Orientation::create_Pr()
{
    double T = (JD(t) - 2451545.0) / 36525;
    Matrix<double> m(3, 3);
    double ksi = (2.650545 + (2306.083227 + (0.2988499 + (0.01801828 - (0.000005791 +
                    + 0.0000003173 * T) * T)*T) * T) * T) / 206264.806247097;
    double Z = (-2.650545 + (2306.077181 + (1.0927348 + (0.0182683 - (0.000028596 +
                    + 0.0000002904 * T) * T) * T) * T) * T) / 206264.806247097;
    double theta = (2004.191903 - (0.4294934 + (0.04182264 + (0.000007089 +
                    + 0.0000001274 * T) * T) * T) * T) * T / 206264.806247097;
    m(0, 0) =  std::cos(ksi)*std::cos(Z)*std::cos(theta) - std::sin(ksi)*std::sin(Z);
    m(0, 1) = -std::sin(ksi)*std::cos(Z)*std::cos(theta) - std::cos(ksi)*std::sin(Z);
    m(0, 2) = -std::cos(Z)*std::sin(theta);
    
    m(1, 0) =  std::cos(ksi)*std::sin(Z)*std::cos(theta) + std::sin(ksi)*std::cos(Z);
    m(1, 1) = -std::sin(ksi)*std::sin(Z)*std::cos(theta) + std::cos(ksi)*std::cos(Z);
    m(1, 2) = -std::sin(Z)*std::sin(theta);

    m(2, 0) =  std::cos(ksi)*std::sin(theta);
    m(2, 1) = -std::sin(ksi)*std::sin(theta);
    m(2, 2) =  std::cos(theta);
    return m;
}
Matrix<double> Orientation::create_BT()
{
    Matrix<double> m(3, 3);
    m(0, 0) = +0.9999999999999942;
    m(0, 1) = -0.0000000707827974;
    m(0, 2) = +0.0000000805621715;
    
    m(1, 0) = +0.0000000707827948;
    m(1, 1) = +0.9999999999999969;
    m(1, 2) = +0.0000000330604145;
    
    m(2, 0) = -0.0000000805621738;
    m(2, 1) = -0.0000000330604088;
    m(2, 2) = +0.9999999999999962; 
    return m;
}


Orientation::Orientation(double curr_time) : t(curr_time)  
{
    Matrix<double> A = create_Pt();
    Matrix<double> B = create_St();
    Matrix<double> C = create_Nt();
    Matrix<double> D = create_Pr();
    Matrix<double> E = create_BT();
    orientaion = ublas::prod(A, orientaion);
    orientaion = ublas::prod(B, orientaion);
    orientaion = ublas::prod(C, orientaion);
    orientaion = ublas::prod(D, orientaion);
    orientaion = ublas::prod(E, orientaion);
}
Vector<double> Orientation::GCRS_to_ITRS(const Vector<double> &vec)
{
    Vector<double> result = ublas::prod(orientaion, vec);
    return result;
}
Vector<double> Orientation::ITRS_to_GCRS(const Vector<double> &vec)
{
    Vector<double> result = ublas::prod(ublas::trans(orientaion), vec);
    return result;
}
