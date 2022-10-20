//ExampleEvolution.cpp

//This program is a basic example of the usage of the 
//library RGESolver.h to perform the renormalization 
//group evolution in the context of the SMEFT from an 
//high-energy scale Lambda to another scale muLow. 

//For the sake of simplicity, this program considers only 
//strong coupling and top Yukawa couplings. Only one of the 
//dimension-six operators is set =/= 0 at the highscale.

#include <iostream>
#include "RGESolver.h"

int main() {

    RGESolver S; //Initialize an instance of the class

    //Energy scales must be given in GeV 
    double Lambda = 10000.;
    double muLow = 250.;

    //Setting the initial conditions at mu = Lambda 
    S.SetCoefficient("CHG", 1. / (Lambda * Lambda));

    std::cout << "---- INITIAL CONDITIONS ----" << std::endl;
    std::cout << "g3(" << Lambda << " GeV) : "
            << S.GetCoefficient("g3") << std::endl;
    std::cout << "yt(" << Lambda << " GeV) : "
            << S.GetCoefficient("YuR", 2, 2) << std::endl;
    std::cout << "CHG(" << Lambda << " GeV) : "
            << S.GetCoefficient("CHG") << std::endl;

   
    //Evolution from Lambda to muLow via numerical solution 
    //of the RGEs
    S.Evolve("Numeric", Lambda, muLow);


    std::cout << "--- AFTER THE EVOLUTION ---" << std::endl;
    std::cout << "g3(" << muLow << " GeV) : "
            << S.GetCoefficient("g3") << std::endl;
    std::cout << "yt(" << muLow << " GeV) : "
            << S.GetCoefficient("YuR", 2, 2) << std::endl;
    std::cout << "CHG(" << muLow << " GeV) : "
            << S.GetCoefficient("CHG") << std::endl;

    return 0;
}
