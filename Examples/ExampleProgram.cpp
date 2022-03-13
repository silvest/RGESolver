//ExampleProgram.cpp
#include <iostream>
#include "RGESolver.h"

int main() {
    RGESolver S; 
    double Lambda = 10000.;
    
    S.SetCoefficient("CG", 1. / (Lambda * Lambda));
    S.SetCoefficient("g3", 0.8);
    S.SetCoefficient("YuR", 0.95, 2, 2);
    
    S.Evolve("Numeric", Lambda, 250.);
    
    std::cout << "CG : " << S.GetCoefficient("CG") << std::endl;
    std::cout << "g3 : " << S.GetCoefficient("g3") << std::endl;
    std::cout << "Re[Yu(3,3)] : "
              << S.GetCoefficient("YuR", 2, 2) << std::endl;
	
    S.SaveOutputFile("OutputSMEFTEvolution.txt", "SLHA");
	
	return 0;
}
