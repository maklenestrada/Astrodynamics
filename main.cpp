#include "MathOperations.h"
#include "AstroCalculations.h"
#include <iostream>

using namespace std;
int main() {

    double r[3];
    double v[3];
    r[0] = 0.159321004;
    r[1] = 0.579266185;
    r[2] = 0.052359609;

    v[0] = -9.30360303251;
    v[1] = 3.018641330;
    v[2] = 1.536362143;

    double a = 0;
    double e = 0;
    double p = 0;
    double Omega = 0;
    double i = 0;
    double omega = 0;
    double f = 0;
    AstroCalculations AstroCalc;
    AstroCalc.rvtoOrbitElements(r,v,&a,&e,&p,&Omega,&i,&omega,&f);

    cout << "Semi-major axis (a): " << a << endl;
    cout << "Eccentricity (e): " << e << endl;
    cout << "Semi-latus rectum (p): " << p << endl;
    cout << "RAAN (\u03A9): " << Omega << "\u00B0" << endl;  // \u03A9 is Greek capital letter Omega (Ω)
    cout << "Inclination (i): " << i << "\u00B0" << endl;
    cout << "Argument of periapsis (\u03C9): " << omega << "\u00B0" << endl;  // \u03C9 is Greek small letter omega (ω)
    cout << "True anomaly (f): " << f << "\u00B0" << endl;

    return 0;
}
