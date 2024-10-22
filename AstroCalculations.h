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

    //Function to calculate new position vector at a later time using F&G series approach
    void F_and_G_Series(double r_new[3],double ro[3],double v0[3],double t1, double t2);
};


#endif //ASTRODYNAMICS_ASTROCALCULATIONS_H
