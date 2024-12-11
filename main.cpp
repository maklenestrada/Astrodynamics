#include "MathOperations.h"
#include "AstroCalculations.h"
#include "AttitudeRepresentation.h"
#include <iostream>

using namespace std;
int main() {

    AttitudeRepresentation AttRep;
    double R[3][3] = {
            {0.21822,0.48507,-0.84681},
            {0.43644,0.72761,0.52926},
            {0.87287,-0.48507,-0.05293}
    };
    double eps[3] = {0};
    double eta = 0;
    AttRep.QuaternionParameter(R,&eps,&eta);
    cout << "Quaternion Vector Part " << endl;
    for(int i = 0; i < 3; i++)
    {
        cout << eps[i] << endl;
    }
    cout << "Quaternion Scalar Part: " << eta << endl;

    AttitudeRepresentation attitude;
    double nR[3][3] = {0};
    string seq = "ZYZ";  // Rotation sequence
    double phi = .56;
    double theta = .78;
    double psi = .123;

    // Compute the rotation matrix
    attitude.EulerAngle_3S_RotationMatrix(&nR, seq, phi, theta, psi);

    // Print the resulting matrix R
    cout << "Matrix R:" << endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cout << R[i][j] << " ";
        }
        cout << endl;
    }

    double R1[3][3] = {0};
    string sseq = "Z";
    attitude.EulerAngle_1S_RotationMatrix(&R1,sseq,theta);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cout << R1[i][j] << " ";
        }
        cout << endl;
    }
//    double Axis[3] = {0};
//    double Angle = 0;
//    AttRep.AxisAngleParameters(R,&Axis,&Angle);
//    cout << "Axis of Rotation" << endl;
//    for(int i = 0; i < 3; i++)
//    {
//        cout << Axis[i] << endl;
//    }
//    cout << "Angle of Rotation: " << Angle << endl;

//    double r[3];
//    double v[3];
//    r[0] = 0.159321004;
//    r[1] = 0.579266185;
//    r[2] = 0.052359609;
//
//    v[0] = -9.30360303251;
//    v[1] = 3.018641330;
//    v[2] = 1.536362143;
//
//    double a = 0;
//    double e = 0;
//    double p = 0;
//    double Omega = 0;
//    double i = 0;
//    double omega = 0;
//    double f = 0;
//    AstroCalculations AstroCalc;
//    AstroCalc.rvtoOrbitElements(r,v,&a,&e,&p,&Omega,&i,&omega,&f);
//
//    cout << "Semi-major axis (a): " << a << endl;
//    cout << "Eccentricity (e): " << e << endl;
//    cout << "Semi-latus rectum (p): " << p << endl;
//    cout << "RAAN (\u03A9): " << Omega << "\u00B0" << endl;  // \u03A9 is Greek capital letter Omega (Ω)
//    cout << "Inclination (i): " << i << "\u00B0" << endl;
//    cout << "Argument of periapsis (\u03C9): " << omega << "\u00B0" << endl;  // \u03C9 is Greek small letter omega (ω)
//    cout << "True anomaly (f): " << f << "\u00B0" << endl;
//
//    double r_new[3] = {0};
//    double t1 = 0.0105076712;
//    double t2 = 0.021370777;
//    AstroCalc.F_and_G_Series(r_new,r,v,t1, t2);
//    cout << "Position at time T2 using F&G series" << endl;
//    cout << r_new[0] << " i + " << endl;
//    cout << r_new[1] << " j +" << endl;
//    cout << r_new[2] << " k" <<endl;

    return 0;
}
