#ifndef ASTRODYNAMICS_ATTITUDEREPRESENTATION_H
#define ASTRODYNAMICS_ATTITUDEREPRESENTATION_H
#include <iostream>

using namespace std;
class AttitudeRepresentation {
public:
    //Constructor
    AttitudeRepresentation() {}

    //Function to construct rotation matrix given Quaternion parameters
    void QuaternionRotationMatrix(double (*R)[3][3],double epsilon[3],double eta);

    //Function to extract Quaternion vector (epsilon) and constant (eta) from a Rotation Matrix
    void QuaternionParameter(double R[3][3], double (*epsilon)[3],double* eta);

    //Function to form rotation matrix from AxisAngle
    void AxisAngleRotationMatrix(double (*R)[3][3], double Axis[3], double Angle);

    //Function to extract axis and angle from rotation matrix
    void AxisAngleParameters(double R[3][3],double (*Axis)[3], double* Angle);

    //Function to construct a 3 sequence Euler Angle Rotation Matrix
    void EulerAngle_3S_RotationMatrix(double (*R)[3][3], const string& seq, double phi, double theta,double psi);

    //Function to construct a 1 sequence  Euler Angle Rotation Matrix
    void EulerAngle_1S_RotationMatrix(double (*R)[3][3],const string& seq, double theta);
};


#endif //ASTRODYNAMICS_ATTITUDEREPRESENTATION_H
