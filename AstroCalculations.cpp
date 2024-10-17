#include "AstroCalculations.h"
#include "MathOperations.h"
#include<math.h>
#include <cmath>

#define mu  4*M_PI*M_PI
#define mu1 1/(4*M_PI*M_PI)
#define D2R M_PI/180.0
#define R2D 180/M_PI

//Function to Calculate Orbit Elements Given Position and Velocity
void AstroCalculations::rvtoOrbitElements(double r[3],double v[3], double* a, double* e, double* p, double* Omega, double* i, double* omega, double* f)
{
    MathOperations MathOps;

    //Step 1 Calculate the norm of r and v
    double r_norm,v_norm;
    MathOps.Norm(&r_norm, r);
    MathOps.Norm(&v_norm, v);

    //Step 2 Calculate angular momentum
    double h[3],h_norm;
    MathOps.CrossProduct(h,r,v);
    MathOps.Norm(&h_norm,h);

    //Step 3 Calculate inclination
    *i = acos(h[2]/h_norm)*R2D;

    //Step 4 Calculate line of nodes vector N
    double k[3] = {0,0,1};
    double N[3],N_norm;
    MathOps.CrossProduct(N,k,h);
    MathOps.Norm(&N_norm,N);

    //Step 5 Calculate the Ecentricity
    double S1[3], S2, S3[3], e_norm;
    double S4[3];
    double e_vec[3];
    MathOps.CrossProduct(S1, v, h);
    S2 = mu / r_norm;
    MathOps.MultVecByScalar(S3, r, S2);
    MathOps.SubtractVectors(S4, S1, S3);
    MathOps.MultVecByScalar(e_vec, S4, mu1);
    MathOps.Norm(&e_norm, e_vec);
    *e = e_norm;

    //Step 6 Calculate RAAN
    if (N[1] > 0) {
        *Omega = acos(N[0] / N_norm)*R2D;
    }else {
        *Omega = 360 - acos(N[0] / N_norm)*R2D;
    }

    //Step 7 Calculate Argument of Periapse
    double Ne,Ntimese;
    MathOps.DotProduct(&Ne,N,e_vec);
    Ntimese = N_norm * e_norm;
    if(e_vec[3] > 0){
        *omega = acos(Ne/Ntimese)*R2D;
    }else {
        *omega = 360 - acos(Ne/Ntimese)*R2D;
    }

    //Step 8 Calculate True Anomoly
        //Substep Calculate Radial Vel
        double rdotv,v_r;
        MathOps.DotProduct(&rdotv,r,v);
        v_r = rdotv/r_norm;
    double edotr,etimesr;
    MathOps.DotProduct(&edotr,e_vec,r);
    etimesr = e_norm * r_norm;
    if (v_r > 0) {
        *f = acos(edotr / etimesr) * R2D;
    } else
        *f = 360 - acos(edotr / etimesr) * R2D;


    //Step 9 Calculate semi-major axis
    double epsilon; //specific orbital energy
    epsilon = (v_norm*v_norm)/2 - mu/r_norm;
    *a = -mu/(2*epsilon);

    //Step 10 Calculate Parameter P
    double a_val;
    a_val = *a;
    *p = a_val*(1-e_norm*e_norm);


}