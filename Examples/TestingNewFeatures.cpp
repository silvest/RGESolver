
#include <iostream>
#include "RGESolver.h"
#include <chrono>
//#include <array>
using namespace std::chrono;

int main(int argc, char **argv) {
    using namespace std;
    RGESolver S;




    double t12, t13, t23, dc;
    t12 = 0.007;
    t13 = 0.52;
    t23 = 0.11;
    dc = 0.0013;

    S.SetCoefficient("g3", 1.2);
    S.SetCoefficient("g2", .62);
    S.SetCoefficient("g1", .3);



    S.SetCoefficient("CKM_theta12", t12);
    S.SetCoefficient("CKM_theta13", t13);
    S.SetCoefficient("CKM_theta23", t23);
    S.SetCoefficient("CKM_delta", dc);
    cout << "Angles and phase (input): "
            << S.GetCoefficient("CKM_theta12") << "\t"
            << S.GetCoefficient("CKM_theta13") << "\t"
            << S.GetCoefficient("CKM_theta23") << "\t"
            << S.GetCoefficient("CKM_delta") << endl;

    double mutmp = 124.;
    S.SetSMInputScale(mutmp);
    S.GenerateSMInitialConditions(mutmp * 1.1,
            "UP", "Leading-Log");


    cout << "Angles and phase (output): "
            << S.GetCoefficient("CKM_theta12") << "\t"
            << S.GetCoefficient("CKM_theta13") << "\t"
            << S.GetCoefficient("CKM_theta23") << "\t"
            << S.GetCoefficient("CKM_delta") << endl;


    return 0;
}
