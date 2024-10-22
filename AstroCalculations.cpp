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

void AstroCalculations::F_and_G_Series(double r_new[3],double ro[3],double vo[3],double t1,double t2)
{
    //Define MathOperations Object
    MathOperations MathOps;
    // Define epsilon, lambda, psi
    double r;
    MathOps.Norm(&r,ro);
    double eps = -mu / (pow(r, 3));
    double lam1;
    MathOps.DotProduct(&lam1, ro, vo);
    double lam = lam1/ (r * r);
    double psi1;
    MathOps.DotProduct(&psi1,vo,vo);
    double psi = psi1 / (r * r);

    double tau = t2 - t1;

    // Tau/factorial
    double T2 = pow(tau, 2) / tgamma(2 + 1);
    double T3 = pow(tau, 3) / tgamma(3 + 1);
    double T4 = pow(tau, 4) / tgamma(4 + 1);
    double T5 = pow(tau, 5) / tgamma(5 + 1);
    double T6 = pow(tau, 6) / tgamma(6 + 1);

    // r calc
    double r2 = eps * T2;
    double r3 = 3 * eps * lam * T3;
    double r4 = (-15 * eps * pow(lam, 2) + 3 * eps * psi - 2 * pow(eps, 2)) * T4;
    double r5 = (105 * eps * pow(lam, 3) - 45 * eps * lam * psi + 30 * pow(eps, 2) * lam) * T5;
    double r6 = (-945 * eps * pow(lam, 4) + 630 * eps * pow(lam, 2) * psi - 420 * pow(eps, 2) *
                pow(lam, 2) + 66 * pow(eps, 2) * psi - 45 * eps * pow(psi, 2) + 30 * pow(eps, 3)) * T6;

    double rseries = 1 - r2 + r3 + r4 + r5 + r6;

    // v calc
    double v3 = eps * T3;
    double v4 = (6 * eps * lam) * T4;
    double v5 = (-45 * eps * pow(lam, 2) + 9 * eps * psi - 8 * pow(eps, 2)) * T5;
    double v6 = (-180 * eps * lam * psi + 150 * pow(eps, 2) * lam + 315 * eps * pow(lam, 3)) * T6;

    double vseries = tau - v3 + v4 + v5 + v6;

    // New r calculation
    for (int i = 0; i < 3; ++i) {
        r_new[i] = rseries * ro[i] + vseries * vo[i];
    }
}