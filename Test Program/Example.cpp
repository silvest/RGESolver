
#include <iostream>
#include "RGESolver.h"
#include <chrono>
//#include <array>
using namespace std::chrono;

int main(int argc, char **argv) {

    RGESolver S;
    double Lambda = 10000.;
    double a, b, c, d, e, f, g, h;
    auto start = high_resolution_clock::now();
    int Nrep = atoi(argv[1]);
    std::string method = argv[2];

    for (int k = 0; k < Nrep; k ++) {
        std::cout << "Rep n#: " << k << std::endl;
        auto part1 = high_resolution_clock::now();
        for (int i = 0; i < 300; i ++) {
            S.SetCoefficient("g3", 0.8);
            S.SetCoefficient("CG", 1. / (Lambda * Lambda));

            S.SetCoefficient("YuR", 0.95, 2, 2);
            S.SetCoefficient("CeHR", 0.45 / (Lambda * Lambda), 1, 2);
            S.SetCoefficient("CeeR", 0.5 / (Lambda * Lambda), 1, 1, 2, 1);
            S.SetCoefficient("Cqq1R", 0.7 / (Lambda * Lambda), 2, 2, 2, 1);
            S.SetCoefficient("Clequ1R", 0.6 / (Lambda * Lambda), 1, 2, 0, 1);
            S.SetCoefficient("CeuI", .9 / (Lambda * Lambda), 2, 2, 0, 1);

        }

        S.Evolve(method, Lambda, 250.);

        for (int i = 0; i < 300; i ++) {
            a = S.GetCoefficient("CG");
            b = S.GetCoefficient("g3");
            c = S.GetCoefficient("YuR", 2, 2);
            d = S.GetCoefficient("CeHR", 1, 2);
            e = S.GetCoefficient("CeeR", 1, 1, 2, 1);
            f = S.GetCoefficient("Cqq1R", 2, 2, 2, 1);
            g = S.GetCoefficient("Clequ1R", 1, 2, 0, 1);
            g = S.GetCoefficient("CeuI", 2, 2, 0, 1);

        }
        S.Reset();
        auto part2 = high_resolution_clock::now();
        auto partial = duration_cast<microseconds>(part2 - part1);
        std::cout << "Time for the " << k << "-th repetition : "
                << partial.count()*0.001 << " ms" << std::endl;
    }

    using namespace std;
    cout << a << endl;
    cout << b << endl;
    cout << c << endl;
    cout << d << endl;
    cout << e << endl;
    cout << f << endl;
    cout << g << endl;
    cout << h << endl;


    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);


   // S.SaveOutputFile("OutputSMEFTEvolution.txt", "SLHA");
    std::cout << "Time taken by setting initial conditions and performing "
            "the evolution (mean over " << Nrep << " reps) : "
            << duration.count()*0.001 / Nrep << " ms" << std::endl;


    return 0;
}
