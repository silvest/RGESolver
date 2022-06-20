//ExampleEvolutionBackRotation.cpp

//This program is a realistic example of the usage of the 
//library RGESolver.h to perform the renormalization 
//group evolution in the context of the SMEFT from an 
//high-energy scale Lambda to another scale muLow. 

//The initial conditions at the high-energy scale Lambda 
//are generated using pure SM RGEs. 
//A flavour rotation is then performed to go into a well-definite 
//flavour basis.


#include <iostream>
#include "RGESolver.h"

int main() {
    //Initialize an instance of the class
    RGESolver S;

    //Energy scales must be given in GeV 
    double Lambda = 1001.; //Lambda = 10 TeV
    double muLow = 1000.;

    //Generate SM initial conditions using a custom low-energy input.

    //SM fermions masses
    double Muin[3] = {.002, 1.2, 170.};
    double Mdin[3] = {.006, .05, 5.2};
    double Mein[3] = {.005, .1, 1.2};

    //CKM angles and phase
    double t12in = asin(0.225);
    double t13in = asin(0.042);
    double t23in = asin(0.003675);
    double deltain = 1.167625;

    //gauge couplings and Higgs sector parameters
    double g1in = .35;
    double g2in = .65;
    double g3in = 1.2;
    double lambdain = 0.14;
    double mh2in = 15625.;

    //Energy scale at which the input is generated 
    double SMinputScale = 91.0;

    //Generating the initial conditions for the SM 
    //parameters at mu = Lambda in the up basis.
    /*S.GenerateSMInitialConditions(SMinputScale, Lambda, "DOWN", "Numeric",
            g1in, g2in, g3in, lambdain, mh2in,
            Muin, Mdin, Mein,
            t12in, t13in, t23in, deltain);*/
    //Generate SM initial conditions using the default input parameters. 
    S.GenerateSMInitialConditions(Lambda, "DOWN", "Numeric");


    //Initial conditions for SMEFT coefficients
    S.SetCoefficient("CHdR", 1. / (Lambda * Lambda), 1, 2);
    S.SetCoefficient("Cqu1R", 2. / (Lambda * Lambda), 1, 2, 2, 2);
    S.SetCoefficient("CllR", 2. / (Lambda * Lambda), 1, 1, 1, 2);

    S.SetCoefficient("CeHR", 1. / (Lambda * Lambda), 1, 1);
    S.SetCoefficient("CeHR", 1. / (Lambda * Lambda), 0, 0);
    S.SetCoefficient("CeHR", 1. / (Lambda * Lambda), 2, 2);

    S.SetCoefficient("CHG", 1.2 / (Lambda * Lambda));


    std::cout << "---- SM INITIAL CONDITIONS ----" << std::endl;
    std::cout << "g1(" << Lambda << " GeV) : "
            << S.GetCoefficient("g1") << std::endl;
    std::cout << "g2(" << Lambda << " GeV) : "
            << S.GetCoefficient("g2") << std::endl;
    std::cout << "g3(" << Lambda << " GeV) : "
            << S.GetCoefficient("g3") << std::endl;
    std::cout << "lambda(" << Lambda << " GeV) : "
            << S.GetCoefficient("lambda") << std::endl;
    std::cout << "mh^2(" << Lambda << " GeV) : "
            << S.GetCoefficient("mh2") << " GeV^2" << std::endl;

    int i, j;

    std::cout << "Yu(" << Lambda << " GeV) : " << std::endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            std::cout << S.GetCoefficient("YuR", i, j) << " + I * " <<
                    S.GetCoefficient("YuI", i, j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Yd(" << Lambda << " GeV) : " << std::endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            std::cout << S.GetCoefficient("YdR", i, j) << " + I * " <<
                    S.GetCoefficient("YdI", i, j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Ye(" << Lambda << " GeV) : " << std::endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            std::cout << S.GetCoefficient("YeR", i, j) << " + I * " <<
                    S.GetCoefficient("YeI", i, j) << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << "---- SMEFT INITIAL CONDITIONS ----" << std::endl;
    std::cout << "CHdR(1,2)(" << Lambda << " GeV) : "
            << S.GetCoefficient("CHdR", 1, 2) << std::endl;
    std::cout << "Cqu1R(1,2,2,2)(" << Lambda << " GeV) : "
            << S.GetCoefficient("Cqu1R", 1, 2, 2, 2) << std::endl;
    std::cout << "CllR(1,1,1,2)(" << Lambda << " GeV) : "
            << S.GetCoefficient("CllR", 1, 1, 1, 2) << std::endl;
    std::cout << "CHG(" << Lambda << " GeV) : "
            << S.GetCoefficient("CHG") << std::endl;

    std::cout << "CeHR(0,0)(" << Lambda << " GeV) : "
            << S.GetCoefficient("CeHR", 0, 0) << std::endl;
    std::cout << "CeHR(1,1)(" << Lambda << " GeV) : "
            << S.GetCoefficient("CeHR", 1, 1) << std::endl;
    std::cout << "CeHR(2,2)(" << Lambda << " GeV) : "
            << S.GetCoefficient("CeHR", 2, 2) << std::endl;

    //Evolution from Lambda to muLow (solving numerically the RGEs), where 
    //a back-rotation is performed to go back into the up basis
    S.EvolveToBasis("Numeric", Lambda, muLow, "DOWN");
    // S.Evolve("Numeric", Lambda, muLow);




    std::cout << "---------------------------------" << std::endl;
    std::cout << "---------------------------------" << std::endl;

    std::cout << "---- SM AFTER THE EVOLUTION ----" << std::endl;
    std::cout << "g1(" << muLow << " GeV) : "
            << S.GetCoefficient("g1") << std::endl;
    std::cout << "g2(" << muLow << " GeV) : "
            << S.GetCoefficient("g2") << std::endl;
    std::cout << "g3(" << muLow << " GeV) : "
            << S.GetCoefficient("g3") << std::endl;
    std::cout << "lambda(" << muLow << " GeV) : "
            << S.GetCoefficient("lambda") << std::endl;
    std::cout << "mh^2(" << muLow << " GeV) : "
            << S.GetCoefficient("mh2") << " GeV^2" << std::endl;

    std::cout << "Yu(" << muLow << " GeV) : " << std::endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            std::cout << S.GetCoefficient("YuR", i, j) << " + I * " <<
                    S.GetCoefficient("YuI", i, j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Yd(" << muLow << " GeV) : " << std::endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            std::cout << S.GetCoefficient("YdR", i, j) << " + I * " <<
                    S.GetCoefficient("YdI", i, j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Ye(" << muLow << " GeV) : " << std::endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            std::cout << S.GetCoefficient("YeR", i, j) << " + I * " <<
                    S.GetCoefficient("YeI", i, j) << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << "---- SMEFT AFTER THE EVOLUTION ----" << std::endl;
    std::cout << "CHdR(1,2)(" << muLow << " GeV) : "
            << S.GetCoefficient("CHdR", 1, 2) << std::endl;
    std::cout << "Cqu1R(1,2,2,2)(" << muLow << " GeV) : "
            << S.GetCoefficient("Cqu1R", 1, 2, 2, 2) << std::endl;
    std::cout << "CllR(1,1,1,2)(" << muLow << " GeV) : "
            << S.GetCoefficient("CllR", 1, 1, 1, 2) << std::endl;
    std::cout << "CHG(" << muLow << " GeV) : "
            << S.GetCoefficient("CHG") << std::endl;


    std::cout << "CeHR(0,0)(" << Lambda << " GeV) : "
            << S.GetCoefficient("CeHR", 0, 0) << std::endl;
    std::cout << "CeHR(1,1)(" << Lambda << " GeV) : "
            << S.GetCoefficient("CeHR", 1, 1) << std::endl;
    std::cout << "CeHR(2,2)(" << Lambda << " GeV) : "
            << S.GetCoefficient("CeHR", 2, 2) << std::endl;
    return 0;
}
