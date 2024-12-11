#include "AttitudeRepresentation.h"
#include "MathOperations.h"
#include<math.h>
#include <cmath>
#include <iostream>
using namespace std;

//This function forms a rotation matrix given the Quaternion vector (epsilon) and scalar (eta)
void AttitudeRepresentation::QuaternionRotationMatrix(double (*R)[3][3], double epsilon[3], double eta)
{
    double R11,R12,R13,R21,R22,R23,R31,R32,R33;
    double eps1,eps2,eps3;
    eps1 = epsilon[0];
    eps2 = epsilon[1];
    eps3 = epsilon[2];

    R11 = eta*eta + eps1*eps1 - eps2*eps2 - eps3*eps3;
    R12 = 2*(eps1*eps2 + eta*eps3);
    R13 = 2*(eps1*eps3 - eta*eps2);
    R21 = 2*(eps1*eps2 - eta*eps3);
    R22 = eta*eta - eps1*eps1 + eps2*eps2 - eps3*eps3;
    R23 = 2*(eps2*eps3 + eta*eps1);
    R31 = 2*(eps1*eps3 + eta*eps2);
    R32 = 2*(eps2*eps3 - eta*eps1);
    R33 = eta*eta - eps1*eps1 - eps2*eps2 + eps3*eps3;

    (*R)[0][0] = R11;
    (*R)[0][1] = R12;
    (*R)[0][2] = R13;
    (*R)[1][0] = R21;
    (*R)[1][1] = R22;
    (*R)[1][2] = R23;
    (*R)[2][0] = R31;
    (*R)[2][1] = R32;
    (*R)[2][2] = R33;
}

//Function to extract Quaternion vector (epsilon) and constant (eta) from a Rotation Matrix
//Extracts +constant
void AttitudeRepresentation::QuaternionParameter(double R[3][3], double (*epsilon)[3], double* eta)
{
    double eps0,eps1,eps2,eps3;
    eps0 = (0.5)*sqrt(R[0][0] + R[1][1] + R[2][2] + 1);
    eps1 = (R[1][2] - R[2][1])/(4*eps0);
    eps2 = (R[2][0] - R[0][2])/(4*eps0);
    eps3 = (R[0][1] - R[1][0])/(4*eps0);

    *eta = eps0;
    (*epsilon)[0] = eps1;
    (*epsilon)[1] = eps2;
    (*epsilon)[2] = eps3;
}

//Function to form rotation matrix from AxisAngle
void AttitudeRepresentation::AxisAngleRotationMatrix(double (*R)[3][3], double Axis[3], double Angle)
{
    // Ensuring axis is a unit vector
    MathOperations MathOps;
    double Anorm;
    double a[3] = {0};
    double ax,ay,az;
    MathOps.Norm(&Anorm,Axis);
    MathOps.MultVecByScalar(a,Axis,1/Anorm);
    ax = a[0];
    ay = a[1];
    az = a[2];

    double c,s,v;
    c = cos(Angle);
    s = sin(Angle);
    v = 1 - c;

    double R11,R12,R13,R21,R22,R23,R31,R32,R33;
    R11 = ax*ax*v + c;
    R12 = ax*ay*v - az*s;
    R13 = ax*az*v + ay*s;
    R21 = ax*ay*v + az*s;
    R22 = ay*ay*v + c;
    R23 = ay*az*v - ax*s;
    R31 = ax*az*v - ay*s;
    R32 = ay*az*v + ax*s;
    R33 = az*az*v + c;

    (*R)[0][0] = R11;
    (*R)[0][1] = R12;
    (*R)[0][2] = R13;
    (*R)[1][0] = R21;
    (*R)[1][1] = R22;
    (*R)[1][2] = R23;
    (*R)[2][0] = R31;
    (*R)[2][1] = R32;
    (*R)[2][2] = R33;
}

//This function extract the axis and angle of a rotation matrix
void AttitudeRepresentation::AxisAngleParameters(double R[3][3], double (*Axis)[3], double* Angle)
{
    double theta = acos((R[0][0] + R[1][1] + R[2][2] - 1)/2);
    double ax,ay,az;
    ax = (R[2][1] - R[1][2])/(2*sin(theta));
    ay = (R[0][2] - R[2][0])/(2*sin(theta));
    az = (R[1][0] - R[0][1])/(2*sin(theta));

    *Angle = theta;
    (*Axis)[0] = ax;
    (*Axis)[1] = ay;
    (*Axis)[2] = az;
}

//Function to construct Euler Angle Rotation Matrix
void AttitudeRepresentation::EulerAngleRotationMatrix(double (*R)[3][3], const string& seq, double phi, double theta, double psi)
{
    double c1,s1,c2,s2,c3,s3;
    c1 = cos(phi);
    s1 = sin(phi);
    c2 = cos(theta);
    s2 = sin(theta);
    c3 = cos(psi);
    s3 = sin(psi);

    double RX[3][3] = {
            {1,0,0},
            {0,c1,s1},
            {0,-s1,c1}
    };

    double RY[3][3] = {
            {c2,0,-s2},
            {0,1,0},
            {s2,0,c2}
    };

    double RZ[3][3] = {
            {c3,s3,0},
            {-s3,c3,0},
            {0,0,1}
    };
    MathOperations MathOps;
    double ans[3][3] = {0};
    // Determine the sequence and multiply the matrices accordingly
    //Antisymmetric
    if (seq == "ZYX") {
        MathOps.multiplyThreeMatrices(RZ, RY, RX, &ans);
    } else if (seq == "XYZ") {
        MathOps.multiplyThreeMatrices(RX, RY, RZ, &ans);
    } else if (seq == "YXZ") {
        MathOps.multiplyThreeMatrices(RY, RX, RZ, &ans);
    } else if (seq == "XZY") {
        MathOps.multiplyThreeMatrices(RX, RZ, RY, &ans);
    } else if (seq == "ZXY") {
        MathOps.multiplyThreeMatrices(RZ, RX, RY, &ans);
    } else if (seq == "YZX")
    //Symmetric
    {
        MathOps.multiplyThreeMatrices(RY, RZ, RX, &ans);
    } else if (seq == "XYX") {
        MathOps.multiplyThreeMatrices(RX, RY, RX, &ans);
    } else if (seq == "YZY") {
        MathOps.multiplyThreeMatrices(RY, RZ, RY, &ans);
    } else if (seq == "ZXZ") {
        MathOps.multiplyThreeMatrices(RZ, RX, RZ, &ans);
    } else if (seq == "XZX") {
        MathOps.multiplyThreeMatrices(RX, RZ, RX, &ans);
    } else if (seq == "YXY") {
        MathOps.multiplyThreeMatrices(RY, RX, RY, &ans);
    } else if (seq == "ZYZ") {
        MathOps.multiplyThreeMatrices(RZ, RY, RZ, &ans);
    } else {
        throw invalid_argument("Unsupported rotation sequence.");
    }

    // Assign values from ans to R using a for loop
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            (*R)[i][j] = ans[i][j];
        }
    }
}