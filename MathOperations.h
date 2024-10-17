//
// Created by Maklen Estrada on 10/16/24.
//

#ifndef ASTRODYNAMICS_MATHOPERATIONS_H
#define ASTRODYNAMICS_MATHOPERATIONS_H


class MathOperations {
public:
    //Constructor
    MathOperations() {}

    //Function for cross product
    void CrossProduct(double ans[3],double X[3],double Y[3]);

    //Function for dot product
    void DotProduct(double* ans,double X[3],double Y[3]);

    //Function to calculate norm of vector
    void Norm(double* ans,double vec[3]);

    //Function to multiply vector times scalar
    void MultVecByScalar(double ans[3], double vec[3], double c);

    //Function to add two vectors
    void AddVectors(double ans[3], double x[3],double y[3]);

    //Function to subtract two vectors
    void SubtractVectors(double ans[3], double x[3],double y[3]);

    //Function to check if value is near a desired tolerance
    int ValueNear(double val, double goal, double tol);
};


#endif //ASTRODYNAMICS_MATHOPERATIONS_H
