
#include <iostream>
#include "RGESolver.h"
#include <chrono>
//#include <array>
using namespace std::chrono;

int main(int argc, char **argv) {
    using namespace std;
    RGESolver S;





    S.SetCoefficient("g3", 1.2);
    S.SetCoefficient("g2", .6);
    S.SetCoefficient("g1", .3);

    S.SetCoefficient("lambda", .2);
    S.SetCoefficient("mh2", 126. * 126.);



    S.SetCoefficient("CKM_theta12",1.);
    S.SetCoefficient("CKM_theta13", 1.);
    S.SetCoefficient("CKM_theta23", 1.);
    S.SetCoefficient("CKM_delta",0.35);
    cout << "Angles and phase (input): "
            << S.GetCoefficient("CKM_theta12") << "\t"
            << S.GetCoefficient("CKM_theta13") << "\t"
            << S.GetCoefficient("CKM_theta23") << "\t"
            << S.GetCoefficient("CKM_delta") << endl;

    double muIn = 200.;
    double muFin = 200.;
    S.SetSMInputScale(muIn);
    S.GenerateSMInitialConditions(muFin,
            "DOWN", "Leading-Log");


    cout << "Angles and phase (output): "
            << S.GetCoefficient("CKM_theta12") << "\t"
            << S.GetCoefficient("CKM_theta13") << "\t"
            << S.GetCoefficient("CKM_theta23") << "\t"
            << S.GetCoefficient("CKM_delta") << endl;
    int i, j;
    cout << "Yd : " << endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            cout << S.GetCoefficient("YdR", i, j)
                    << " + i*" << S.GetCoefficient("YdI", i, j)
                    << "\t";
        }
        cout << " " << endl;
    }
    cout << "Yu : " << endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            cout << S.GetCoefficient("YuR", i, j)
                    << " + i*" << S.GetCoefficient("YuI", i, j)
                    << "\t";
        }
        cout << " " << endl;
    }
    return 0;
}
