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



    /*
    S.SetSMInputScale(250.);

    S.SetCoefficient("g1", 0.3);
    S.SetCoefficient("g2", 0.6);
    S.SetCoefficient("g3", 1.2);
    S.SetCoefficient("lambda", 0.2);
    S.SetCoefficient("mh2", 16000.);

    S.SetCKMAngle("CKM_theta12", 0.1);
    S.SetCKMAngle("CKM_theta13", 0.2);
    S.SetCKMAngle("CKM_theta23", 0.3);
    S.SetCKMPhase(0.1);

    S.SetFermionMass("mu", 0.002);
    S.SetFermionMass("mc", 1.2);
    S.SetFermionMass("mt", 1.2);

    S.SetFermionMass("md", 0.002);
    S.SetFermionMass("ms", 1.2);
    S.SetFermionMass("mb", 1.2);

    S.SetFermionMass("mel", 0.006;);
    S.SetFermionMass("mmu", 0.100);
    S.SetFermionMass("mtau", 1.2);
     */

    //Generate SM initial conditions using the default input parameters. 
    S.GenerateSMInitialConditions(Lambda, "DOWN", "Numeric");


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
    S.Evolve("Numeric", Lambda, muLow);

    cout << "-------- AFTER THE EVOLUTION --------" << endl;
    cout << "g3(" << muLow << " GeV) : " << S.GetCoefficient("g3") << endl;
    cout << "CG(" << muLow << " GeV) : " << S.GetCoefficient("CG") << endl;
    cout << "CHdR[" << i << "," << j << "](" << muLow << " GeV) : "
            << S.GetCoefficient("CHdR", i, j) << endl;

    //Compute the masses of the fermions and the CKM matrix
    //parameters.
    S.ComputeCKMAndFermionMasses();


    cout << "md(" << muLow << " GeV) : " << S.GetFermionMass("md") << " GeV" << endl;
    cout << "CKM_theta13(" << muLow << " GeV) : " << S.GetCKMAngle("CKM_theta13") << " rad" << endl;
    cout << "-----------------------------------------" << endl;



    S.Reset();




    return 0;
}
