#include "MathOperations.h"
#include "MathOperations.h"
#include <vector>
#include <cmath>

//Computes cross product
void MathOperations::CrossProduct(double ans[3],double X[3],double Y[3])
{
    ans[0] = X[1] * Y[2] - X[2] * Y[1];
    ans[1] = X[2] * Y[0] - X[0] * Y[2];
    ans[2] = X[0] * Y[1] - X[1] * Y[0];
}

//Computes dot product
void MathOperations::DotProduct(double* ans,double X[3],double Y[3])
{
    if (ans != nullptr){
        *ans = X[0] * Y[0] + X[1] * Y[1] + X[2] * Y[2];
    }

}

//Calculates norm
void MathOperations::Norm(double* ans,double vec[3])
{
    if (ans != nullptr){
        *ans = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    }
}

//Multiply vector by a scalor
void MathOperations::MultVecByScalar(double ans[3], double vec[3], double c)
{
    for(int i = 0; i <3; i++){
        ans[i] = vec[i] * c;
    }
}

//Add two vectors together
void MathOperations::AddVectors(double ans[3], double x[3],double y[3])
{
    for(int i = 0; i < 3; i++){
        ans[i] = x[i] + y[i];
    }
}

//Subtract two vectors together
void MathOperations::SubtractVectors(double ans[3], double x[3],double y[3])
{
    for(int i = 0; i < 3; i++){
        ans[i] = x[i] - y[i];
    }
}

//Check if value is near desired tolerance
int MathOperations::ValueNear(double val, double goal, double tol)
{
    if ((val > goal - tol) && (val < goal + tol))
        return 1;
    else
        return 0;
}