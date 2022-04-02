//Example1.cpp

//This program generates the initial conditions for SM 
//parameters at the high-scale Lambda in the down basis and 
//then computes the running down to muLow.

#include <iostream>
#include "RGESolver.h"

int main() {
    using namespace std;
    //Initialize an instance of the class
    RGESolver S;

    double Lambda = 10000; //Lambda = 10 TeV
    double muLow = 250.;
    int i, j;
    i = 0;
    j = 1;




    //Generate SM initial conditions using the default input parameters. 
    // S.GenerateSMInitialConditions(Lambda, "DOWN", "Numeric");

    //Generate SM initial conditions using a custom low-energy input.
    double Muin[3] = {.002, 1.2, 170.};
    double Mdin[3] = {.006, .05, 5.2};
    double Mein[3] = {.005, .1, 1.2};

    double g1in = .31;
    double g2in = .6;
    double g3in = 1.2;
    double lambdain = 0.2;
    double mh2in = 126. * 126.;
    double t12in = 0.2;
    double t13in = .1;
    double t23in = .3;
    double deltain = 3.14 / 4.;

    S.GenerateSMInitialConditions(91., Lambda, "DOWN", "Numeric",
            g1in, g2in, g3in, lambdain, mh2in,
            Muin, Mdin, Mein,
            t12in, t13in, t23in, deltain);

    //Initial conditions for SMEFT coefficients
    S.SetCoefficient("CG", 1. / (Lambda * Lambda));
    S.SetCoefficient("CHdR", 2. / (Lambda * Lambda), i, j);

    cout << "---- INITIAL CONDITIONS AT HIGHSCALE ----" << endl;
    cout << "-----------------------------------------" << endl;
    cout << "g3(" << Lambda << " GeV) : " << S.GetCoefficient("g3") << endl;
    cout << "CG(" << Lambda << " GeV) : " << S.GetCoefficient("CG") << endl;
    cout << "CHdR[" << i << "," << j << "](" << Lambda << " GeV) : "
            << S.GetCoefficient("CHdR", i, j) << endl;
    cout << "-----------------------------------------\n" << endl;

    //Evolution from Lambda to muLow (solving numerically the RGEs)
    S.EvolveToBasis("Numeric", Lambda, muLow, "UP");

    cout << "-------- AFTER THE EVOLUTION --------" << endl;
    cout << "g3(" << muLow << " GeV) : " << S.GetCoefficient("g3") << endl;
    cout << "CG(" << muLow << " GeV) : " << S.GetCoefficient("CG") << endl;
    cout << "CHdR[" << i << "," << j << "](" << muLow << " GeV) : "
            << S.GetCoefficient("CHdR", i, j) << endl;

    //Compute the masses of the fermions and the CKM matrix
    //parameters.


    cout << "md(" << muLow << " GeV) : " << S.GetFermionMass("md") << " GeV" << endl;
    cout << "CKM_theta13(" << muLow << " GeV) : " << S.GetCKMAngle("CKM_theta13") << " rad" << endl;

    cout << "-----------------------------------------" << endl;


    cout << "yd(" << muLow << " GeV) : " << endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            cout << S.GetCoefficient("YdR", i, j) << " + I * " <<
                    S.GetCoefficient("YdI", i, j) << "\t";
        }
        cout << "" << endl;
    }

    cout << "yu(" << muLow << " GeV) : " << endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            cout << S.GetCoefficient("YuR", i, j) << " + I * " <<
                    S.GetCoefficient("YuI", i, j) << "\t";
        }
        cout << "" << endl;
    }
    
    
    S.Reset();

    return 0;
}
