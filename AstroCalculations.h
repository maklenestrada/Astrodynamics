//
// Created by Maklen Estrada on 10/16/24.
//

#ifndef ASTRODYNAMICS_ASTROCALCULATIONS_H
#define ASTRODYNAMICS_ASTROCALCULATIONS_H
#include <iostream>
using namespace std;
class AstroCalculations {
public:
    //Constructor
    AstroCalculations() {}

    //Function to Calculate Orbit Elements Given Position and Velocity
    void rvtoOrbitElements(double r[3],double v[3], double* a, double* e, double* p, double* Omega, double *i, double* omega, double* f);

};


#endif //ASTRODYNAMICS_ASTROCALCULATIONS_H
