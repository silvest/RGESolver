
#include <iostream>
#include "RGESolver.h"
#include <chrono>
//#include <array>
using namespace std::chrono;

int main(int argc, char **argv) {
    using namespace std;

    RGESolver S;



    //[xxx0] ok
    //[xx00] ok
    //[x0x0] ok
    //[0xx0] ok


    //[xxxx] ok
    //[xx0x] not ok
    //[x0xx] not ok
    //[0xxx] ok



    S.SetCoefficient("g3", 1.2);
    S.SetCoefficient("g2", .6);
    S.SetCoefficient("g1", .3);

    S.SetCoefficient("lambda", .2);
    S.SetCoefficient("mh2", 126. * 126.);


    S.SetCKMAngle("CKM_theta12", 0.4);
    S.SetCKMAngle("CKM_theta13", 0.32);
    S.SetCKMAngle("CKM_theta23", 0.1);

    S.SetCKMPhase(.004);
    cout << "Angles and phase (input): "
            << S.GetCKMAngle("CKM_theta12") << "\t"
            << S.GetCKMAngle("CKM_theta13") << "\t"
            << S.GetCKMAngle("CKM_theta23") << "\t"
            << S.GetCKMPhase() << endl;

    cout << "Masses (input) : \n"
            << S.GetFermionMass("mu") << "\t"
            << S.GetFermionMass("mc") << "\t"
            << S.GetFermionMass("mt") << "\n"
            << S.GetFermionMass("md") << "\t"
            << S.GetFermionMass("ms") << "\t"
            << S.GetFermionMass("mb") << "\n"
            << S.GetFermionMass("mel") << "\t"
            << S.GetFermionMass("mmu") << "\t"
            << S.GetFermionMass("mtau") << endl;

    double muIn = 200.;
    double muFin = 200.;
    S.SetSMInputScale(muIn);
    S.GenerateSMInitialConditions(muFin,
            "DOWN", "Leading-Log");


    cout << "Angles and phase (output): "
            << S.GetCKMAngle("CKM_theta12") << "\t"
            << S.GetCKMAngle("CKM_theta13") << "\t"
            << S.GetCKMAngle("CKM_theta23") << "\t"
            << S.GetCKMPhase() << endl;
    cout << "Masses (output) : \n"
            << S.GetFermionMass("mu") << "\t"
            << S.GetFermionMass("mc") << "\t"
            << S.GetFermionMass("mt") << "\n"
            << S.GetFermionMass("md") << "\t"
            << S.GetFermionMass("ms") << "\t"
            << S.GetFermionMass("mb") << "\n"
            << S.GetFermionMass("mel") << "\t"
            << S.GetFermionMass("mmu") << "\t"
            << S.GetFermionMass("mtau") << endl;




    S.Evolve("Leading-Log", muFin, muFin *1.0002);
    S.ComputeCKMAndFermionMasses();
    cout << "Angles and phase (after evolve): "
            << S.GetCKMAngle("CKM_theta12") << "\t"
            << S.GetCKMAngle("CKM_theta13") << "\t"
            << S.GetCKMAngle("CKM_theta23") << "\t"
            << S.GetCKMPhase() << endl;
    cout << "Masses (after evolve) : \n"
            << S.GetFermionMass("mu") << "\t"
            << S.GetFermionMass("mc") << "\t"
            << S.GetFermionMass("mt") << "\n"
            << S.GetFermionMass("md") << "\t"
            << S.GetFermionMass("ms") << "\t"
            << S.GetFermionMass("mb") << "\n"
            << S.GetFermionMass("mel") << "\t"
            << S.GetFermionMass("mmu") << "\t"
            << S.GetFermionMass("mtau") << endl;




    return 0;
}
