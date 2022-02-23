
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



    S.SetCoefficient("CKM_theta12", 0.6);
    S.SetCoefficient("CKM_theta13", 0.5);
    S.SetCoefficient("CKM_theta23", 0.2);
    S.SetCoefficient("CKM_delta", 1.);
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
    cout << "\n\n----------------------\n\n" << endl;


    RGESolver S3;

    double L = 10000.;
    double L2 = L*L;
    S3.SetCoefficient("CG", 1. / L2);
    S3.SetCoefficient("g3", 1.4);
    cout << "CG,g3 (after setting) " << S3.GetCoefficient("CG")
            << ",\t" << S3.GetCoefficient("g3")
            << ",\t" << endl;
    S3.Evolve("Leading-Log", S3.GetSMInputScale(), L);
    cout << "CG,g3 (after evolution) " << S3.GetCoefficient("CG")
            << ",\t" << S3.GetCoefficient("g3")
            << ",\t" << endl;
    S3.Reset();
    cout << "CG,g3 (after reset) " << S3.GetCoefficient("CG")
            << ",\t" << S3.GetCoefficient("g3")
            << ",\t" << endl;


    return 0;
}
