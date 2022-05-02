//#include <gsl/gsl_linalg.h>

//#include "gsl/gsl_complex.h"
//#include "gsl/gsl_complex_math.h"

/*
void RGESolver::SetCKMPhase(double val) {
    if (val <= 3.141592653589793
            && val>- 3.141592653589793) {
        CKM_delta = val;
    } else {
        std::cout << "ERROR: CKM PHASE SHOULD BE IN THE INTERVAL (-pi,pi]"
                << std::endl;
    }
}

void RGESolver::SetCKMAngle(std::string name, double val) {
    if (val <= 3.141592653589793 * 0.5
            && val >= 0.) {
 *(CKMAngles.at(name)) = val;
    } else {
        std::cout << "ERROR: CKM ANGLES SHOULD BE IN THE INTERVAL [0,pi/2]"
                << std::endl;
    }
}*/

double RGESolver::GetCKMAngle(std::string name) {
    return * (CKMAngles.at(name));
}

double RGESolver::GetCKMPhase() {
    return CKM_delta;
}

/*
void RGESolver::SetFermionMass(std::string name, double val) {
    if (val >= 0.) {
 *(FermionMasses.at(name)) = val;
    } else {
        std::cout << "ERROR: FERMION MASSES MUST BE NON-NEGATIVE"
                << std::endl;
    }
}*/

double RGESolver::GetFermionMass(std::string name) {
    return * (FermionMasses.at(name));
}



//Extracts the 4 parameters from the CKM 
//affected by unphysical phases

void RGESolver::ExtractParametersFromCKM() {

    s13 = CKM(0, 2).abs();
    c13 = (sqrt(gslpp::complex(1. - s13 * s13, 0.))).real();
    c12 = CKM(0, 0).abs() / c13;
    s12 = (sqrt(gslpp::complex(1. - c12 * c12, 0.))).real();
    s23 = CKM(1, 2).abs() / c13;
    c23 = (sqrt(gslpp::complex(1. - s23 * s23, 0.))).real();
    CKM_theta12 = (arcsin(gslpp::complex(s12, 0.))).real();
    CKM_theta13 = (arcsin(gslpp::complex(s13, 0.))).real();
    CKM_theta23 = (arcsin(gslpp::complex(s23, 0.))).real();

    //See http://www.utfit.org/UTfit/Formalism
    //delta in [-Pi,Pi]
    double gamma = (- (CKM(0, 0)*(CKM(0, 2)).conjugate()) /
            (CKM(1, 0)*(CKM(1, 2)).conjugate())).arg();

    double a = (c12 * s13 * s23) / (s12 * c23);
    /* std::cout << "s12 and c12 : "
            << s12 << "  " << c12 << std::endl;
    std::cout << "s13 and c13 : "
            << s13 << "  " << c13 << std::endl;
    std::cout << "s23 and c23 : "
            << s23 << "  " << c23 << std::endl;
    std::cout << "a : " << a << std::endl;
    std::cout << "gamma : " << gamma << std::endl;*/

    double pi = 3.14159265358979323846;

    if (abs(gamma) == 0 || gamma == pi ||
            s12 == 0 || c12 == 0 ||
            s13 == 0 || c13 == 0 ||
            s23 == 0 || c23 == 0) {
        CKM_delta = 0.;
    } else {
        double tan_g = tan(gamma);
        if (gamma < pi * 0.5 && gamma > - pi * 0.5) {
            CKM_delta = 2. * atan(
                    (1. - sqrt(1. - (a * a - 1.) * tan_g * tan_g)) /
                    ((a - 1.) * tan_g)
                    );
        } else {
            CKM_delta = 2. * atan(
                    (1. + sqrt(1. - (a * a - 1.) * tan_g * tan_g)) /
                    ((a - 1.) * tan_g)
                    );
        }
    }


}

void RGESolver::FromMassesToYukawas(std::string basis) {
    double v = sqrt(0.5 * mh2 / lambda);
    UpdateCKM();

    double yudiag[3], yddiag[3], yediag[3];
    double sqrt2 = sqrt(2);

    yudiag[0] = sqrt2 * mu / v;
    yudiag[1] = sqrt2 * mc / v;
    yudiag[2] = sqrt2 * mt / v;

    yddiag[0] = sqrt2 * md / v;
    yddiag[1] = sqrt2 * ms / v;
    yddiag[2] = sqrt2 * mb / v;

    yediag[0] = sqrt2 * mel / v;
    yediag[1] = sqrt2 * mmu / v;
    yediag[2] = sqrt2 * mtau / v;
    int i, j;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            yuR[i][j] = 0;
            yuI[i][j] = 0;
            ydR[i][j] = 0;
            ydI[i][j] = 0;
            yeR[i][j] = 0;
            yeI[i][j] = 0;
        }
    }
    if (basis == "UP") {
        gslpp::matrix<gslpp::complex> CKMdag = CKM.hconjugate();
        for (i = 0; i < 3; i ++) {
            yuR[i][i] = yudiag[i];
            yeR[i][i] = yediag[i];
            for (j = 0; j < 3; j ++) {
                ydR[i][j] = yddiag[i]*(CKMdag(i, j)).real();
                ydI[i][j] = yddiag[i]*(CKMdag(i, j)).imag();
            }
        }


    }
    if (basis == "DOWN") {
        for (i = 0; i < 3; i ++) {
            ydR[i][i] = yddiag[i];
            yeR[i][i] = yediag[i];
            for (j = 0; j < 3; j ++) {
                yuR[i][j] = yudiag[i]*(CKM(i, j)).real();
                yuI[i][j] = yudiag[i]*(CKM(i, j)).imag();
            }
        }
    }



}

void RGESolver::UpdateCKM() {
    using namespace std;
    gslpp::complex PhaseFactor(1., CKM_delta,
            true);
    //First row
    CKM.assign(0, 0, c12 * c13);
    CKM.assign(0, 1, s12 * c13);
    CKM.assign(0, 2, s13 * PhaseFactor.conjugate());

    //Second row 
    CKM.assign(1, 0, - s12 * c23 - c12 * s23 * s13 * PhaseFactor);
    CKM.assign(1, 1, c12 * c23 - s12 * s23 * s13 * PhaseFactor);
    CKM.assign(1, 2, s23 * c13);

    //Third row
    CKM.assign(2, 0, s12 * s23 - c12 * c23 * s13 * PhaseFactor);
    CKM.assign(2, 1, - c12 * s23 - s12 * c23 * s13 * PhaseFactor);
    CKM.assign(2, 2, c23 * c13);
}

void RGESolver::GenerateSMInitialConditions(
        double mu, std::string basis,
        std::string method) {
    /*if (method != "Numeric" && method != "Leading-Log") {
        std::cout << "WARNING : invalid method\n"
                "Available methods: Numeric, Leading-Log"
                << std::endl;
    }*/
    SetSMDefaultInput();
    //Before evolving, eventual changes 
    //in the input values of CKM angles must be 
    //translated 
    c12 = cos(CKM_theta12);
    s12 = sin(CKM_theta12);
    c13 = cos(CKM_theta13);
    s13 = sin(CKM_theta13);
    c23 = cos(CKM_theta23);
    s23 = sin(CKM_theta23);


    //if (inputCKM == true) {
    UpdateCKM();
    FromMassesToYukawas(basis);
    // }

    EvolveSMOnly(method, InputScale_SM, mu);
    GoToBasisSMOnly(basis);
}

void RGESolver::GenerateSMInitialConditions(double muIn, double muFin, std::string basis, std::string method,
        double g1in, double g2in, double g3in, double lambdain, double mh2in,
        double Muin[3], double Mdin[3], double Mein[3],
        double t12in, double t13in, double t23in, double deltain) {


    double pihalf = 3.141592653589793;

    if (Muin[0] > 0. && Muin[1] > 0. && Muin[2] > 0. &&
            Mdin[0] > 0. && Mdin[1] > 0. && Mdin[2] > 0. &&
            Mein[0] > 0. && Mein[1] > 0. && Mein[2] > 0. &&
            t12in <= pihalf && t12in >= 0. &&
            t13in <= pihalf && t13in >= 0. &&
            t23in <= pihalf && t23in >= 0. &&
            deltain <= 3.141592653589793 && deltain > - 3.141592653589793
            ) {

        g1 = g1in;
        g2 = g2in;
        g3 = g3in;

        mh2 = mh2in;
        lambda = lambdain;

        mu = Muin[0];
        mc = Muin[1];
        mt = Muin[2];

        md = Mdin[0];
        ms = Mdin[1];
        mb = Mdin[2];

        mel = Mein[0];
        mmu = Mein[1];
        mtau = Mein[2];

        CKM_theta12 = t12in;
        CKM_theta13 = t13in;
        CKM_theta23 = t23in;
        CKM_delta = deltain;


        c12 = cos(CKM_theta12);
        s12 = sin(CKM_theta12);
        c13 = cos(CKM_theta13);
        s13 = sin(CKM_theta13);
        c23 = cos(CKM_theta23);
        s23 = sin(CKM_theta23);


        UpdateCKM();
        FromMassesToYukawas(basis);

        EvolveSMOnly(method, muIn, muFin);
        GoToBasisSMOnly(basis);
    } else {
        std::cout << "WARNING : Invalid parameters. Generation of initial conditions not completed." << std::endl;
    }
}

void RGESolver::EvolveSMOnly(std::string method, double muI, double muF) {
    if (method != "Numeric" && method != "Leading-Log") {
        std::cout << "WARNING : invalid method\n"
                "Available methods: Numeric, Leading-Log"
                << std::endl;
    }
    InitSMOnly();
    //Initial conditions are inserted 
    //in the array x


    //Numeric resolution 

    if (method == "Numeric") {
        double tI = log(muI);
        double tF = log(muF);
        double ttmp = tI;

        //Initial step 
        step_ = (tF - tI)*0.1;

        if (tF > tI) {

            while (ttmp < tF) {
                int status = gsl_odeiv2_evolve_apply(evoSMOnly_, conSMOnly_,
                        sSMOnly_, &sysSMOnly_, &ttmp, tF, &step_, x);
                if (status != GSL_SUCCESS) {
                    printf("error in adaptive integration, return value=%d\n", status);
                }
            }

        } else {

            while (ttmp > tF) {
                int status = gsl_odeiv2_evolve_apply(evoSMOnly_, conSMOnly_,
                        sSMOnly_, &sysSMOnly_, &ttmp, tF, &step_, x);
                if (status != GSL_SUCCESS) {
                    printf("error in adaptive integration, return value=%d\n", status);
                }

            }
        }


        gsl_odeiv2_evolve_reset(evoSMOnly_);
        gsl_odeiv2_step_reset(sSMOnly_);

    }

    //Leading-log resolution 
    //-------------------------------------------
    if (method == "Leading-Log") {
        double beta[59] = {0.};
        double Log_muF_over_muI = log(muF / muI);
        /*int status = */funcSMOnly(10., x, beta, NULL);
        for (int i = 0; i < 59; i ++) {
            x[i] += beta[i] * Log_muF_over_muI;
        }
    }
    //-------------------------------------------

    UpdateSMOnly();
    //Evolved values from x are put 
    //back in the coefficients 

}

void RGESolver::GoToBasisSMOnly(std::string basis) {

    gslpp::matrix<gslpp::complex> Uu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yuDiag(3, 3, 0.);
    gslpp::vector<double> Su(3, 0.);

    gslpp::matrix<gslpp::complex> Ud(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ydDiag(3, 3, 0.);
    gslpp::vector<double> Sd(3, 0.);

    gslpp::matrix<gslpp::complex> Re(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rl(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Redag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rldag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yeDiag(3, 3, 0.);
    gslpp::vector<double> Se(3, 0.);

    int i, j;
    gslpp::matrix<gslpp::complex> yu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ye(3, 3, 0.);

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            yu.assign(i, j, gslpp::complex(yuR[i][j], yuI[i][j], false));
            yd.assign(i, j, gslpp::complex(ydR[i][j], ydI[i][j], false));
            ye.assign(i, j, gslpp::complex(yeR[i][j], yeI[i][j], false));
        }
    }

    using namespace std;

    yu.singularvalue(Uu, Vu, Su);
    yd.singularvalue(Ud, Vd, Sd);
    ye.singularvalue(Re, Rl, Se);

    //Updating fermion masses

    //mf = yfdiag * v /sqrt(2)
    //v = sqrt(mh2 / (2 lambda))
    //v/sqrt(2) =  sqrt(mh2 / (4 lambda))
    double vOverSqrt2 = sqrt(0.25 * mh2 / lambda);


    mu = Su(0) * vOverSqrt2;
    mc = Su(1) * vOverSqrt2;
    mt = Su(2) * vOverSqrt2;

    md = Sd(0) * vOverSqrt2;
    ms = Sd(1) * vOverSqrt2;
    mb = Sd(2) * vOverSqrt2;

    mel = Se(0) * vOverSqrt2;
    mmu = Se(1) * vOverSqrt2;
    mtau = Se(2) * vOverSqrt2;


    //Matrix to rotate fields
    gslpp::matrix<gslpp::complex> Ru(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rudag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rddag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rq(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rqdag(3, 3, 0.);

    //Computing the CKM 
    CKM = (Vu.hconjugate()) * Vd;
    //Extract the 4 parameters from the raw CKM
    ExtractParametersFromCKM();
    //Build the CKM with the 4 parameters
    UpdateCKM();

    //Reset the Yukawas to 0
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            yuR[i][j] = 0;
            yuI[i][j] = 0;
            ydR[i][j] = 0;
            ydI[i][j] = 0;
            yeR[i][j] = 0;
            yeI[i][j] = 0;
        }
    }

    if (basis == "UP") {
        //In UP basis: 
        //yu = diag(.,.,.)
        //yd = diag(.,.,.) CKM^dagger 
        Ru = Uu;
        Rd = Ud;
        Rq = Vu;

        for (i = 0; i < 3; i ++) {
            yuR[i][i] = Su(i);
            yeR[i][i] = Se(i);
            for (j = 0; j < 3; j ++) {
                ydR[i][j] = Sd(i)*(CKM(j, i)).real();
                ydI[i][j] = - Sd(i)*(CKM(j, i)).imag();
            }
        }


    }
    if (basis == "DOWN") {
        //In DOWN basis: 
        //yd = diag(.,.,.)
        //yu = diag(.,.,.) CKM 
        Ru = Uu;
        Rd = Ud;
        Rq = Vd;

        for (i = 0; i < 3; i ++) {
            ydR[i][i] = Sd(i);
            yeR[i][i] = Se(i);
            for (j = 0; j < 3; j ++) {
                yuR[i][j] = Su(i)*(CKM(i, j)).real();
                yuI[i][j] = Su(i)*(CKM(i, j)).imag();
            }
        }


    }
}

void RGESolver::InitSMOnly() {

    int n = 0;
    int i, j;
    x[0] = g2;
    x[1] = g1;
    x[2] = g3;
    n += 3;
    x[n] = lambda;
    n ++;
    x[n] = mh2;
    n ++;

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            int count = 0;
            x[n + count * DF] = yuR[i][j];
            count ++;
            x[n + count * DF] = yuI[i][j];
            count ++;
            x[n + count * DF] = ydR[i][j];
            count ++;
            x[n + count * DF] = ydI[i][j];
            count ++;
            x[n + count * DF] = yeR[i][j];
            count ++;
            x[n + count * DF] = yeI[i][j];
            count ++;

            n ++;
        }
    }
    n += (2. * Nyukawa - 1.) * DF;

}

void RGESolver::UpdateSMOnly() {
    //Here we return to original structures

    int n = 0;
    int a, i, j;
    g2 = x[0];
    g1 = x[1];
    g3 = x[2];
    n += 3;
    lambda = x[n];
    n ++;
    mh2 = x[n];
    n ++;

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            yuR[i][j] = x[n + a * DF];
            a ++;
            yuI[i][j] = x[n + a * DF];
            a ++;
            ydR[i][j] = x[n + a * DF];
            a ++;
            ydI[i][j] = x[n + a * DF];
            a ++;
            yeR[i][j] = x[n + a * DF];
            a ++;
            yeI[i][j] = x[n + a * DF];
            a ++;
            n ++;
        }
    }
    n += (2. * Nyukawa - 1.) * DF;
}

int RGESolver::funcSMOnly(double logmu, const double y[], double f[], void* params) {
    int i, j, a, b;
    double loop_factor = 1.
            / (16. * 3.1415926535 * 3.1415926535);
    int c = 0;
    //counter initialized at 0. This index reads inside the array y, where all 
    //independent parameters are stored. In this part of the function the 
    //variables are organized in the correct structures.

    double g1 = y[c + 1]; //gauge couplings 
    double g2 = y[c ];
    double g3 = y[c + 2];

    c += Ngauge;

    double lambda = y[c]; //Higgs sector
    double mh2 = y[c + 1];
    c += Nh;

    double yeR[NG][NG], yeI[NG][NG], //yukawas
            ydR[NG][NG], ydI[NG][NG],
            yuR[NG][NG], yuI[NG][NG];
    double yedagR[NG][NG], yedagI[NG][NG],
            yddagR[NG][NG], yddagI[NG][NG],
            yudagR[NG][NG], yudagI[NG][NG];

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            yuR[i][j] = y[c];
            yudagR[j][i] = y[c];
            a ++;
            yuI[i][j] = y[c + DF];
            yudagI[j][i] = - y[c + DF];
            a ++;
            ydR[i][j] = y[c + 2 * DF];
            yddagR[j][i] = y[c + 2 * DF];
            a ++;
            ydI[i][j] = y[c + 3 * DF];
            yddagI[j][i] = - y[c + 3 * DF];
            a ++;
            yeR[i][j] = y[c + 4 * DF];
            yedagR[j][i] = y[c + 4 * DF];
            a ++;
            yeI[i][j] = y[c + 5 * DF];
            yedagI[j][i] = - y[c + 5 * DF];
            a ++;
            c ++;
        }
    }
    c += (2. * Nyukawa - 1.) * DF;

    //Auxiliary quantities: 

    //Powers and products of gauge couplings
    double g12 = g1*g1;
    double g22 = g2*g2;
    double g32 = g3*g3;
    //double g1g3 = g1*g3;
    //double g1g2 = g1*g2;
    //double g2g3 = g2*g3;

    double gammaH = 0.; //Higgs wavefunction normalization  
    double H = 0.; //scalar yukawa-trace dependent appearing in lambda RGE
    /*double gqR[NG][NG] = {{0.}}; //wavefunction ren. const. of q 
    double gqI[NG][NG] = {{0.}};
    double glR[NG][NG] = {{0.}}; //wavefunction ren. const. of l 
    double glI[NG][NG] = {{0.}};
    double guR[NG][NG] = {{0.}}; //wavefunction ren. const. of u
    double guI[NG][NG] = {{0.}};
    double gdR[NG][NG] = {{0.}}; //wavefunction ren. const. of d 
    double gdI[NG][NG] = {{0.}};
    double geR[NG][NG] = {{0.}}; //wavefunction ren. const. of e 
    double geI[NG][NG] = {{0.}};*/
    double yudyuR[NG][NG] = {
        {0.}
    }; //yu^dag yu 
    double yudyuI[NG][NG] = {
        {0.}
    };
    double yddydR[NG][NG] = {
        {0.}
    }; //yd^dag yd 
    double yddydI[NG][NG] = {
        {0.}
    };
    double yedyeR[NG][NG] = {
        {0.}
    }; //ye^dag ye 
    double yedyeI[NG][NG] = {
        {0.}
    };
    double ydyudR[NG][NG] = {
        {0.}
    }; //yd yu^dag
    double ydyudI[NG][NG] = {
        {0.}
    };
    //double yuyddR[NG][NG] = {{0.}}; //yu yd^dag
    //double yuyddI[NG][NG] = {{0.}};

    //3 Yukawa matrices product
    /*
    double yuyudyuR[NG][NG] = {{0.}}; //yu yu^dag yu
    double yuyudyuI[NG][NG] = {{0.}};
    double yudyuyudR[NG][NG] = {{0.}}; //yu^dag yu yu^dag
    double yudyuyudI[NG][NG] = {{0.}};
    double ydyddydR[NG][NG] = {{0.}}; //yd yd^dag yd
    double ydyddydI[NG][NG] = {{0.}};
    double yddydyddR[NG][NG] = {{0.}}; //yd^dag yd yd^dag
    double yddydyddI[NG][NG] = {{0.}};
    double yeyedyeR[NG][NG] = {{0.}}; //ye ye^dag ye
    double yeyedyeI[NG][NG] = {{0.}};
    double yedyeyedR[NG][NG] = {{0.}}; //ye^dag ye ye^dag
    double yedyeyedI[NG][NG] = {{0.}};

    double yddydyudR[NG][NG] = {{0.}}; // yd^dag yd yu^dag
    double yddydyudI[NG][NG] = {{0.}};
    double yudyuyddR[NG][NG] = {{0.}}; // yu^dag yu yd^dag
    double yudyuyddI[NG][NG] = {{0.}};*/


    //Products of Yukawas 
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            for (a = 0; a < NG; a ++) {

                yudyuR[i][j] += yudagR[i][a] * yuR[a][j] - yudagI[i][a] * yuI[a][j];
                yudyuI[i][j] += yudagI[i][a] * yuR[a][j] + yudagR[i][a] * yuI[a][j];


                yddydR[i][j] += yddagR[i][a] * ydR[a][j] - yddagI[i][a] * ydI[a][j];
                yddydI[i][j] += yddagI[i][a] * ydR[a][j] + yddagR[i][a] * ydI[a][j];


                yedyeR[i][j] += yedagR[i][a] * yeR[a][j] - yedagI[i][a] * yeI[a][j];
                yedyeI[i][j] += yedagI[i][a] * yeR[a][j] + yedagR[i][a] * yeI[a][j];

                ydyudR[i][j] += ydR[i][a] * yudagR[a][j] - ydI[i][a] * yudagI[a][j];
                ydyudI[i][j] += ydI[i][a] * yudagR[a][j] + ydR[i][a] * yudagI[a][j];
            }
        }
    }




    //gammaH
    for (i = 0; i < NG; i ++) {
        gammaH += yedyeR[i][i] + NC * (yudyuR[i][i] + yddydR[i][i]);
        for (j = 0; j < NG; j ++) {
            H += yedyeR[i][j] * yedyeR[j][i] - yedyeI[i][j] * yedyeI[j][i]
                    + NC * (yudyuR[i][j] * yudyuR[j][i] - yudyuI[i][j] * yudyuI[j][i]
                    + yddydR[i][j] * yddydR[j][i] - yddydI[i][j] * yddydI[j][i]);
        }
    }



    //------------------------------------
    //---------------RGE------------------
    //------------------------------------

    c = 0; //counter restarts from 0 

    //---------RGE GAUGE/HIGGS --------

    //SMEFT contributes to SM beta functions are proportional to mh2.
    //They are in RGE 1.
    //SM beta functions for mh2 and lambda are in 
    //https://arxiv.org/abs/hep-ph/0207271 
    //mh2,lambda and Yukawas follow the conventions of RGE 1

    {
        //g2
        f[c] = (- b02 * g22 //SM
                ) * g2 * loop_factor;
        c ++;
        //g1
        f[c] = (- b01 * g12) * g1 * loop_factor;
        c ++;
        //g3
        f[c] = (- b03 * g32) * g3 * loop_factor;
        c ++;

        //lambda
        f[c] = lambda * (24. * lambda - 3. * g12 - 9. * g22 + 4. * gammaH)
                + 0.375 * g12 * g12 + 0.75 * g12 * g22
                + 1.125 * g22 * g22 - 2. * H //SM
                ;
        f[c] *= loop_factor;
        c ++;
        //mh2
        f[c] = mh2 * (12. * lambda + 2. * gammaH - 1.5 * g12 - 4.5 * g22 //SM  
                ) * loop_factor;
        c ++;
    }

    //---------RGE YUKAWA --------
    //SM beta functions can be found in
    // https://arxiv.org/abs/hep-ph/0207271
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            //Entries without matrix products:      
            //yuR  
            f[c ] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yuR[i][j]
                    ;
            //yuI
            f[c + DF] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yuI[i][j]
                    ;
            //ydR     
            f[c + 2 * DF] = (gammaH - (5. / 12.) * g12 - 2.25 * g22
                    - 8. * g32) * ydR[i][j]
                    ;
            //ydI  
            f[c + 3 * DF] = (gammaH - (5. / 12.) * g12 - 2.25 * g22
                    - 8. * g32) * ydI[i][j]
                    ;
            //yeR
            f[c + 4 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * yeR[i][j]
                    ;
            //yeI
            f[c + 5 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * yeI[i][j]
                    ;
            //Entries with 1 matrix product (1 summed index)
            for (b = 0; b < NG; b ++) {
                //yuR
                f[c ] +=
                        + 1.5 * (yuR[i][b]*(yudyuR[b][j] - yddydR[b][j])
                        - yuI[i][b]*(yudyuI[b][j] - yddydI[b][j]))
                        ;
                //yuI
                f[c + DF] += 1.5 * (yuI[i][b]*(yudyuR[b][j] - yddydR[b][j])
                        + yuR[i][b]*(yudyuI[b][j] - yddydI[b][j]))
                        ;
                //ydR       
                f[c + 2 * DF] += 1.5 * (ydR[i][b]*(yddydR[b][j] - yudyuR[b][j])
                        - ydI[i][b]*(yddydI[b][j] - yudyuI[b][j]))
                        ;
                //ydI 
                f[c + 3 * DF] += 1.5 * (
                        ydI[i][b]*(yddydR[b][j] - yudyuR[b][j])
                        + ydR[i][b]*(yddydI[b][j] - yudyuI[b][j]))
                        ;
                //yeR
                f[c + 4 * DF] += 1.5 * (yeR[i][b]*(yedyeR[b][j])
                        - yeI[i][b]*(yedyeI[b][j]))
                        ;
                //yeI 
                f[c + 5 * DF] += 1.5 * (yeI[i][b]*(yedyeR[b][j])
                        + yeR[i][b]*(yedyeI[b][j]))
                        ;

            }


            for (a = 0; a < (2 * Nyukawa); a ++) {
                f[c + a * DF] *= loop_factor;
            }
            c ++;
        }
    }
    c += (2 * Nyukawa - 1) * DF;
    return GSL_SUCCESS;
}


//Default input for SM 

void RGESolver::SetSMDefaultInput() {
    //Default input scale for SM in GeV
    InputScale_SM = 91.;


    //Fermions masses in GeV
    mu = 0.002;
    mc = 1.2;
    mt = 170.;

    md = 0.006;
    ms = 0.05;
    mb = 5.2;

    mel = 0.0005;
    mmu = 0.100;
    mtau = 1.2;


    //CKM parameters in radians 
    CKM_theta12 = asin(0.225);
    CKM_theta13 = asin(0.042);
    CKM_theta23 = asin(0.003675);
    CKM_delta = 1.167625;
    c12 = cos(CKM_theta12);
    s12 = sin(CKM_theta12);
    c13 = cos(CKM_theta13);
    s13 = sin(CKM_theta13);
    c23 = cos(CKM_theta23);
    s23 = sin(CKM_theta23);



    //By default, Yukawas are aligned with CKM input
    //in up basis
    FromMassesToYukawas("UP");
    //Higgs and Gauge sector
    mh2 = 126. * 126.;
    lambda = 0.2;
    g3 = 1.22;
    g2 = 0.6516;
    g1 = .3576;


}

/*
void RGESolver::ComputeCKMAndFermionMasses() {

    gslpp::matrix<gslpp::complex> Uu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yuDiag(3, 3, 0.);
    gslpp::vector<double> Su(3, 0.);

    gslpp::matrix<gslpp::complex> Ud(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ydDiag(3, 3, 0.);
    gslpp::vector<double> Sd(3, 0.);

    gslpp::matrix<gslpp::complex> Re(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rl(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Redag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rldag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yeDiag(3, 3, 0.);
    gslpp::vector<double> Se(3, 0.);

    int i, j;
    gslpp::matrix<gslpp::complex> yu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ye(3, 3, 0.);

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            yu.assign(i, j, gslpp::complex(yuR[i][j], yuI[i][j], false));
            yd.assign(i, j, gslpp::complex(ydR[i][j], ydI[i][j], false));
            ye.assign(i, j, gslpp::complex(yeR[i][j], yeI[i][j], false));
        }
    }

    using namespace std;

    yu.singularvalue(Uu, Vu, Su);
    yd.singularvalue(Ud, Vd, Sd);
    ye.singularvalue(Re, Rl, Se);

    //Updating fermion masses

    //mf = yfdiag * v /sqrt(2)
    //v = sqrt(mh2 / (2 lambda))
    //v/sqrt(2) =  sqrt(mh2 / (4 lambda))
    double vOverSqrt2 = sqrt(0.25 * mh2 / lambda);


    mu = Su(0) * vOverSqrt2;
    mc = Su(1) * vOverSqrt2;
    mt = Su(2) * vOverSqrt2;

    md = Sd(0) * vOverSqrt2;
    ms = Sd(1) * vOverSqrt2;
    mb = Sd(2) * vOverSqrt2;

    mel = Se(0) * vOverSqrt2;
    mmu = Se(1) * vOverSqrt2;
    mtau = Se(2) * vOverSqrt2;


    //Computing the CKM 
    CKM = (Vu.hconjugate()) * Vd;
    //Extract the 4 parameters from the raw CKM
    ExtractParametersFromCKM();
    //Build the CKM with the 4 parameters
    UpdateCKM();

} */




void RGESolver::Reset() {

    int n = 0;
    int a, i, j, k, l;
    g2 = 0;
    g1 = 0;
    g3 = 0;
    n += 3;
    lambda = 0;
    n ++;
    mh2 = 0;
    n ++;

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            yuR[i][j] = 0;
            a ++;
            yuI[i][j] = 0;
            a ++;
            ydR[i][j] = 0;
            a ++;
            ydI[i][j] = 0;
            a ++;
            yeR[i][j] = 0;
            a ++;
            yeI[i][j] = 0;
            a ++;
            n ++;
        }
    }
    n += (2. * Nyukawa - 1.) * DF;

    cG = 0;
    n ++;
    cGT = 0;
    n ++;
    cW = 0;
    n ++;
    cWT = 0;
    n ++;
    cH = 0;
    n ++;
    cHBOX = 0;
    n ++;
    cHD = 0;
    n ++;

    cHG = 0;
    n ++;
    cHB = 0;
    n ++;
    cHW = 0;
    n ++;
    cHWB = 0;
    n ++;
    cHGT = 0;
    n ++;
    cHBT = 0;
    n ++;
    cHWT = 0;
    n ++;
    cHWBT = 0;
    n ++;

    //class 5
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            WC1_set(cuHR, i, j, 0);
            WC1_set(cuHI, i, j, 0);
            a ++;
            WC1_set(cdHR, i, j, 0);
            WC1_set(cdHI, i, j, 0);
            a ++;
            WC1_set(ceHR, i, j, 0);
            WC1_set(ceHI, i, j, 0);
            a ++;
            n ++;
        }
    }
    n += (N5 * 2 - 1) * DF;


    //Class 6
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            WC1_set(ceWR, i, j, 0);
            a ++;
            WC1_set(ceWI, i, j, 0);
            a ++;
            WC1_set(ceBR, i, j, 0);
            a ++;
            WC1_set(ceBI, i, j, 0);
            a ++;
            WC1_set(cuGR, i, j, 0);
            a ++;
            WC1_set(cuGI, i, j, 0);
            a ++;
            WC1_set(cuWR, i, j, 0);
            a ++;
            WC1_set(cuWI, i, j, 0);
            a ++;
            WC1_set(cuBR, i, j, 0);
            a ++;
            WC1_set(cuBI, i, j, 0);
            a ++;
            WC1_set(cdGR, i, j, 0);
            a ++;
            WC1_set(cdGI, i, j, 0);
            a ++;
            WC1_set(cdWR, i, j, 0);
            a ++;
            WC1_set(cdWI, i, j, 0);
            a ++;
            WC1_set(cdBR, i, j, 0);
            a ++;
            WC1_set(cdBI, i, j, 0);
            a ++;
            n ++;
        }
    }
    n += (N6 * 2 - 1) * DF;

    //class 7
    {
        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(cHl1R, WC2R_indices[i][0], WC2R_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(cHl1I, WC2I_indices[i][0], WC2I_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(cHl3R, WC2R_indices[i][0], WC2R_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(cHl3I, WC2I_indices[i][0], WC2I_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(cHeR, WC2R_indices[i][0], WC2R_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(cHeI, WC2I_indices[i][0], WC2I_indices[i][1], 0);
            n ++;
        }

        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(cHq1R, WC2R_indices[i][0], WC2R_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(cHq1I, WC2I_indices[i][0], WC2I_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(cHq3R, WC2R_indices[i][0], WC2R_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(cHq3I, WC2I_indices[i][0], WC2I_indices[i][1], 0);
            n ++;
        }

        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(cHuR, WC2R_indices[i][0], WC2R_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(cHuI, WC2I_indices[i][0], WC2I_indices[i][1], 0);
            n ++;
        }

        for (i = 0; i < DWC2R; i ++) {
            WC2R_set(cHdR, WC2R_indices[i][0], WC2R_indices[i][1], 0);
            n ++;
        }
        for (i = 0; i < DWC2I; i ++) {
            WC2I_set(cHdI, WC2I_indices[i][0], WC2I_indices[i][1], 0);
            n ++;
        }

        for (i = 0; i < NG; i ++) {
            for (j = 0; j < NG; j ++) {
                WC1_set(cHudR, i, j, 0);
                n ++;
            }
        }
        for (i = 0; i < NG; i ++) {
            for (j = 0; j < NG; j ++) {
                WC1_set(cHudI, i, j, 0);
                n ++;
            }
        }


    }
    //class 8_LLLL
    {
        for (a = 0; a < DWC6R; a ++) {
            WC6R_set(cllR, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6I; a ++) {
            WC6I_set(cllI, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6R; a ++) {
            WC6R_set(cqq1R, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6I; a ++) {
            WC6I_set(cqq1I, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6R; a ++) {
            WC6R_set(cqq3R, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6I; a ++) {
            WC6I_set(cqq3I, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
            n ++;
        }

        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(clq1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(clq1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(clq3R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(clq3I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }

    }
    //Class 8_RRRR
    {
        for (a = 0; a < DWC8R; a ++) {
            WC8R_set(ceeR, WC8R_indices[a][0], WC8R_indices[a][1],
                    WC8R_indices[a][2], WC8R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC8I; a ++) {
            WC8I_set(ceeI, WC8I_indices[a][0], WC8I_indices[a][1],
                    WC8I_indices[a][2], WC8I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6R; a ++) {
            WC6R_set(cuuR, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6I; a ++) {
            WC6I_set(cuuI, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6R; a ++) {
            WC6R_set(cddR, WC6R_indices[a][0], WC6R_indices[a][1],
                    WC6R_indices[a][2], WC6R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC6I; a ++) {
            WC6I_set(cddI, WC6I_indices[a][0], WC6I_indices[a][1],
                    WC6I_indices[a][2], WC6I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(ceuR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(ceuI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }


        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cedR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cedI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }

        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cud1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cud1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cud8R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cud8I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }


    }

    //Class 8_LLRR
    {
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cleR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cleI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cluR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cluI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }

        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cldR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cldI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cqeR, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cqeI, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }

        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cqu1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cqu1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cqu8R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cqu8I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }



        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cqd1R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cqd1I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7R; a ++) {
            WC7R_set(cqd8R, WC7R_indices[a][0], WC7R_indices[a][1],
                    WC7R_indices[a][2], WC7R_indices[a][3], 0);
            n ++;
        }
        for (a = 0; a < DWC7I; a ++) {
            WC7I_set(cqd8I, WC7I_indices[a][0], WC7I_indices[a][1],
                    WC7I_indices[a][2], WC7I_indices[a][3], 0);
            n ++;
        }

    }

    //Class 8_LRRL

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            for (k = 0; k < NG; k ++) {
                for (l = 0; l < NG; l ++) {
                    WC5_set(cledqR, i, j, k, l, 0);
                    WC5_set(cledqI, i, j, k, l, 0);
                    n ++;
                }
            }
        }
    }

    n += DF*DF;

    //Class 8_LRLR
    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            for (k = 0; k < NG; k ++) {
                for (l = 0; l < NG; l ++) {
                    a = 0;
                    WC5_set(cquqd1R, i, j, k, l, 0);
                    a ++;
                    WC5_set(cquqd1I, i, j, k, l, 0);
                    a ++;
                    WC5_set(cquqd8R, i, j, k, l, 0);
                    a ++;
                    WC5_set(cquqd8I, i, j, k, l, 0);
                    a ++;
                    WC5_set(clequ1R, i, j, k, l, 0);
                    a ++;
                    WC5_set(clequ1I, i, j, k, l, 0);
                    a ++;
                    WC5_set(clequ3R, i, j, k, l, 0);
                    a ++;
                    WC5_set(clequ3I, i, j, k, l, 0);
                    a ++;
                    n ++;
                }
            }
        }
    }
    n += NG * NG * NG * NG * (2 * N8_LRLR - 1);



}

void RGESolver::GoToBasis(std::string basis) {

    gslpp::matrix<gslpp::complex> Uu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yuDiag(3, 3, 0.);
    gslpp::vector<double> Su(3, 0.);

    gslpp::matrix<gslpp::complex> Ud(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Vd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ydDiag(3, 3, 0.);
    gslpp::vector<double> Sd(3, 0.);

    gslpp::matrix<gslpp::complex> Re(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rl(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Redag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rldag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yeDiag(3, 3, 0.);
    gslpp::vector<double> Se(3, 0.);

    int i, j;
    gslpp::matrix<gslpp::complex> yu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ye(3, 3, 0.);

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            yu.assign(i, j, gslpp::complex(yuR[i][j], yuI[i][j], false));
            yd.assign(i, j, gslpp::complex(ydR[i][j], ydI[i][j], false));
            ye.assign(i, j, gslpp::complex(yeR[i][j], yeI[i][j], false));
        }
    }

    using namespace std;

    yu.singularvalue(Uu, Vu, Su);
    yd.singularvalue(Ud, Vd, Sd);
    ye.singularvalue(Re, Rl, Se);

    //Updating fermion masses

    //mf = yfdiag * v /sqrt(2)
    //v = sqrt(mh2 / (2 lambda))
    //v/sqrt(2) =  sqrt(mh2 / (4 lambda))
    double vOverSqrt2 = sqrt(0.25 * mh2 / lambda);


    mu = Su(0) * vOverSqrt2;
    mc = Su(1) * vOverSqrt2;
    mt = Su(2) * vOverSqrt2;

    md = Sd(0) * vOverSqrt2;
    ms = Sd(1) * vOverSqrt2;
    mb = Sd(2) * vOverSqrt2;

    mel = Se(0) * vOverSqrt2;
    mmu = Se(1) * vOverSqrt2;
    mtau = Se(2) * vOverSqrt2;


    //Matrix to rotate fields
    gslpp::matrix<gslpp::complex> Ru(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rudag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rddag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rq(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rqdag(3, 3, 0.);

    //Computing the CKM 
    CKM = (Vu.hconjugate()) * Vd;
    //Extract the 4 parameters from the raw CKM
    ExtractParametersFromCKM();
    //Build the CKM with the 4 parameters
    UpdateCKM();

    //Reset the Yukawas to 0
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            yuR[i][j] = 0;
            yuI[i][j] = 0;
            ydR[i][j] = 0;
            ydI[i][j] = 0;
            yeR[i][j] = 0;
            yeI[i][j] = 0;
        }
    }

    if (basis == "UP") {
        //In UP basis: 
        //yu = diag(.,.,.)
        //yd = diag(.,.,.) CKM^dagger 
        Ru = Uu;
        Rd = Ud;
        Rq = Vu;

        for (i = 0; i < 3; i ++) {
            yuR[i][i] = Su(i);
            yeR[i][i] = Se(i);
            for (j = 0; j < 3; j ++) {
                ydR[i][j] = Sd(i)*(CKM(j, i)).real();
                ydI[i][j] = - Sd(i)*(CKM(j, i)).imag();
            }
        }


    }
    if (basis == "DOWN") {
        //In DOWN basis: 
        //yd = diag(.,.,.)
        //yu = diag(.,.,.) CKM 
        Ru = Uu;
        Rd = Ud;
        Rq = Vd;

        for (i = 0; i < 3; i ++) {
            ydR[i][i] = Sd(i);
            yeR[i][i] = Se(i);
            for (j = 0; j < 3; j ++) {
                yuR[i][j] = Su(i)*(CKM(i, j)).real();
                yuI[i][j] = Su(i)*(CKM(i, j)).imag();
            }
        }


    }


    //I save the hermitian conjugates of the rotation 
    //matrix to be more efficient. 
    Rudag = Ru.hconjugate();
    Rddag = Rd.hconjugate();
    Rqdag = Rq.hconjugate();

    Redag = Re.hconjugate();
    Rldag = Rl.hconjugate();





    int a, b, c, d, p, r, s, t, n;

    //Auxiliary objects to perform the rotation
    gslpp::complex z, w, x, y, h, f;
    gslpp::complex le, lu, ld, qe, qu1, qu8, qd1, qd8;



    //Coefficients in the new basis

    //Class 5
    double cuHRp[3 * 3] = {0.};
    double cuHIp[3 * 3] = {0.};
    double cdHRp[3 * 3] = {0.};
    double cdHIp[3 * 3] = {0.};
    double ceHRp[3 * 3] = {0.};
    double ceHIp[3 * 3] = {0.};

    //Class 6
    double ceWRp[3 * 3] = {0.};
    double ceWIp[3 * 3] = {0.};
    double ceBRp[3 * 3] = {0.};
    double ceBIp[3 * 3] = {0.};

    double cuWRp[3 * 3] = {0.};
    double cuWIp[3 * 3] = {0.};
    double cuBRp[3 * 3] = {0.};
    double cuBIp[3 * 3] = {0.};
    double cuGRp[3 * 3] = {0.};
    double cuGIp[3 * 3] = {0.};

    double cdWRp[3 * 3] = {0.};
    double cdWIp[3 * 3] = {0.};
    double cdBRp[3 * 3] = {0.};
    double cdBIp[3 * 3] = {0.};
    double cdGRp[3 * 3] = {0.};
    double cdGIp[3 * 3] = {0.};

    //Class 7
    double cHl1Rp[3 * 3] = {0.};
    double cHl1Ip[3 * 3] = {0.};
    double cHl3Rp[3 * 3] = {0.};
    double cHl3Ip[3 * 3] = {0.};
    double cHeRp[3 * 3] = {0.};
    double cHeIp[3 * 3] = {0.};
    double cHq1Rp[3 * 3] = {0.};
    double cHq1Ip[3 * 3] = {0.};
    double cHq3Rp[3 * 3] = {0.};
    double cHq3Ip[3 * 3] = {0.};
    double cHdRp[3 * 3] = {0.};
    double cHdIp[3 * 3] = {0.};
    double cHuRp[3 * 3] = {0.};
    double cHuIp[3 * 3] = {0.};
    double cHudRp[3 * 3] = {0.};
    double cHudIp[3 * 3] = {0.};


    //Class 8
    double cuuRp[81] = {0.};
    double cuuIp[81] = {0.};
    double cddRp[81] = {0.};
    double cddIp[81] = {0.};
    double ceeRp[81] = {0.};
    double ceeIp[81] = {0.};

    double cllRp[81] = {0.};
    double cllIp[81] = {0.};
    double cqq1Rp[81] = {0.};
    double cqq1Ip[81] = {0.};
    double cqq3Rp[81] = {0.};
    double cqq3Ip[81] = {0.};

    double clq1Rp[81] = {0.};
    double clq1Ip[81] = {0.};
    double clq3Rp[81] = {0.};
    double clq3Ip[81] = {0.};

    double ceuRp[81] = {0.};
    double ceuIp[81] = {0.};
    double cedRp[81] = {0.};
    double cedIp[81] = {0.};
    double cud1Rp[81] = {0.};
    double cud1Ip[81] = {0.};
    double cud8Rp[81] = {0.};
    double cud8Ip[81] = {0.};

    double cleRp[81] = {0.};
    double cleIp[81] = {0.};
    double cluRp[81] = {0.};
    double cluIp[81] = {0.};
    double cldRp[81] = {0.};
    double cldIp[81] = {0.};
    double cqeRp[81] = {0.};
    double cqeIp[81] = {0.};

    double cqu1Rp[81] = {0.};
    double cqu1Ip[81] = {0.};
    double cqu8Rp[81] = {0.};
    double cqu8Ip[81] = {0.};
    double cqd1Rp[81] = {0.};
    double cqd1Ip[81] = {0.};
    double cqd8Rp[81] = {0.};
    double cqd8Ip[81] = {0.};


    double cledqRp[81] = {0.};
    double cledqIp[81] = {0.};

    double clequ1Rp[81] = {0.};
    double clequ1Ip[81] = {0.};
    double clequ3Rp[81] = {0.};
    double clequ3Ip[81] = {0.};

    double cquqd1Rp[81] = {0.};
    double cquqd1Ip[81] = {0.};
    double cquqd8Rp[81] = {0.};
    double cquqd8Ip[81] = {0.};

    //Class 5
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cuHR, a, b), WC1(cuHI, a, b)) *
                            Ru(b, j);
                }
            }
            WC1_set(cuHRp, i, j, z.real());
            WC1_set(cuHIp, i, j, z.imag());
        }
    }
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cdHR, a, b), WC1(cdHI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cdHRp, i, j, z.real());
            WC1_set(cdHIp, i, j, z.imag());
        }
    }
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rldag(i, a) *
                            gslpp::complex(
                            WC1(ceHR, a, b), WC1(ceHI, a, b)) *
                            Re(b, j);
                }
            }
            WC1_set(ceHRp, i, j, z.real());
            WC1_set(ceHIp, i, j, z.imag());
        }
    }



    //The coefficients are updated with the rotated ones
    std::copy(std::begin(cuHRp), std::end(cuHRp), std::begin(cuHR));
    std::copy(std::begin(cuHIp), std::end(cuHIp), std::begin(cuHI));
    std::copy(std::begin(cdHRp), std::end(cdHRp), std::begin(cdHR));
    std::copy(std::begin(cdHIp), std::end(cdHIp), std::begin(cdHI));
    std::copy(std::begin(ceHRp), std::end(ceHRp), std::begin(ceHR));
    std::copy(std::begin(ceHIp), std::end(ceHIp), std::begin(ceHI));


    //Class 6
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rldag(i, a) *
                            gslpp::complex(
                            WC1(ceWR, a, b), WC1(ceWI, a, b)) *
                            Re(b, j);
                }
            }
            WC1_set(ceWRp, i, j, z.real());
            WC1_set(ceWIp, i, j, z.imag());
        }
    }
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rldag(i, a) *
                            gslpp::complex(
                            WC1(ceBR, a, b), WC1(ceBI, a, b)) *
                            Re(b, j);
                }
            }
            WC1_set(ceBRp, i, j, z.real());
            WC1_set(ceBIp, i, j, z.imag());
        }
    }


    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cuWR, a, b), WC1(cuWI, a, b)) *
                            Ru(b, j);
                }
            }
            WC1_set(cuWRp, i, j, z.real());
            WC1_set(cuWIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cuBR, a, b), WC1(cuBI, a, b)) *
                            Ru(b, j);
                }
            }
            WC1_set(cuBRp, i, j, z.real());
            WC1_set(cuBIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cuGR, a, b), WC1(cuGI, a, b)) *
                            Ru(b, j);
                }
            }
            WC1_set(cuGRp, i, j, z.real());
            WC1_set(cuGIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cdWR, a, b), WC1(cdWI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cdWRp, i, j, z.real());
            WC1_set(cdWIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cdBR, a, b), WC1(cdBI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cdBRp, i, j, z.real());
            WC1_set(cdBIp, i, j, z.imag());
        }
    }

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rqdag(i, a) *
                            gslpp::complex(
                            WC1(cdGR, a, b), WC1(cdGI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cdGRp, i, j, z.real());
            WC1_set(cdGIp, i, j, z.imag());
        }
    }


    std::copy(std::begin(ceWRp), std::end(ceWRp), std::begin(ceWR));
    std::copy(std::begin(ceWIp), std::end(ceWIp), std::begin(ceWI));
    std::copy(std::begin(ceBRp), std::end(ceBRp), std::begin(ceBR));
    std::copy(std::begin(ceBIp), std::end(ceBIp), std::begin(ceBI));

    std::copy(std::begin(cuWRp), std::end(cuWRp), std::begin(cuWR));
    std::copy(std::begin(cuWIp), std::end(cuWIp), std::begin(cuWI));
    std::copy(std::begin(cuBRp), std::end(cuBRp), std::begin(cuBR));
    std::copy(std::begin(cuBIp), std::end(cuBIp), std::begin(cuBI));
    std::copy(std::begin(cuGRp), std::end(cuGRp), std::begin(cuGR));
    std::copy(std::begin(cuGIp), std::end(cuGIp), std::begin(cuGI));

    std::copy(std::begin(cdWRp), std::end(cdWRp), std::begin(cdWR));
    std::copy(std::begin(cdWIp), std::end(cdWIp), std::begin(cdWI));
    std::copy(std::begin(cdBRp), std::end(cdBRp), std::begin(cdBR));
    std::copy(std::begin(cdBIp), std::end(cdBIp), std::begin(cdBI));
    std::copy(std::begin(cdGRp), std::end(cdGRp), std::begin(cdGR));
    std::copy(std::begin(cdGIp), std::end(cdGIp), std::begin(cdGI));




    //Class 7
    for (n = 0; n < DWC2R; n ++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                z += Rldag(i, a) *
                        gslpp::complex(
                        WC2R(cHl1R, a, b), WC2I(cHl1I, a, b)) *
                        Rl(b, j);
            }
        }
        WC2R_set(cHl1Rp, i, j, z.real());
        WC2I_set(cHl1Ip, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n ++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                z += Rldag(i, a) *
                        gslpp::complex(
                        WC2R(cHl3R, a, b), WC2I(cHl3I, a, b)) *
                        Rl(b, j);
            }
        }
        WC2R_set(cHl3Rp, i, j, z.real());
        WC2I_set(cHl3Ip, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n ++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                z += Redag(i, a) *
                        gslpp::complex(
                        WC2R(cHeR, a, b), WC2I(cHeI, a, b)) *
                        Re(b, j);
            }
        }
        WC2R_set(cHeRp, i, j, z.real());
        WC2I_set(cHeIp, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n ++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                z += Rqdag(i, a) *
                        gslpp::complex(
                        WC2R(cHq1R, a, b), WC2I(cHq1I, a, b)) *
                        Rq(b, j);
            }
        }
        WC2R_set(cHq1Rp, i, j, z.real());
        WC2I_set(cHq1Ip, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n ++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                z += Rqdag(i, a) *
                        gslpp::complex(
                        WC2R(cHq3R, a, b), WC2I(cHq3I, a, b)) *
                        Rq(b, j);
            }
        }
        WC2R_set(cHq3Rp, i, j, z.real());
        WC2I_set(cHq3Ip, i, j, z.imag());
    }

    for (n = 0; n < DWC2R; n ++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                z += Rddag(i, a) *
                        gslpp::complex(
                        WC2R(cHdR, a, b), WC2I(cHdI, a, b)) *
                        Rd(b, j);
            }
        }
        WC2R_set(cHdRp, i, j, z.real());
        WC2I_set(cHdIp, i, j, z.imag());
    }


    for (n = 0; n < DWC2R; n ++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                z += Rudag(i, a) *
                        gslpp::complex(
                        WC2R(cHuR, a, b), WC2I(cHuI, a, b)) *
                        Ru(b, j);
            }
        }
        WC2R_set(cHuRp, i, j, z.real());
        WC2I_set(cHuIp, i, j, z.imag());
    }

    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            z = gslpp::complex(0., 0.);
            for (a = 0; a < 3; a ++) {
                for (b = 0; b < 3; b ++) {
                    z += Rudag(i, a) *
                            gslpp::complex(
                            WC1(cHudR, a, b), WC1(cHudI, a, b)) *
                            Rd(b, j);
                }
            }
            WC1_set(cHudRp, i, j, z.real());
            WC1_set(cHudIp, i, j, z.imag());
        }
    }

    std::copy(std::begin(cHl1Rp), std::end(cHl1Rp), std::begin(cHl1R));
    std::copy(std::begin(cHl1Ip), std::end(cHl1Ip), std::begin(cHl1I));
    std::copy(std::begin(cHl3Rp), std::end(cHl3Rp), std::begin(cHl3R));
    std::copy(std::begin(cHl3Ip), std::end(cHl3Ip), std::begin(cHl3I));
    std::copy(std::begin(cHeRp), std::end(cHeRp), std::begin(cHeR));
    std::copy(std::begin(cHeIp), std::end(cHeIp), std::begin(cHeI));
    std::copy(std::begin(cHq1Rp), std::end(cHq1Rp), std::begin(cHq1R));
    std::copy(std::begin(cHq1Ip), std::end(cHq1Ip), std::begin(cHq1I));
    std::copy(std::begin(cHq3Rp), std::end(cHq3Rp), std::begin(cHq3R));
    std::copy(std::begin(cHq3Ip), std::end(cHq3Ip), std::begin(cHq3I));
    std::copy(std::begin(cHdRp), std::end(cHdRp), std::begin(cHdR));
    std::copy(std::begin(cHdIp), std::end(cHdIp), std::begin(cHdI));
    std::copy(std::begin(cHuRp), std::end(cHuRp), std::begin(cHuR));
    std::copy(std::begin(cHuIp), std::end(cHuIp), std::begin(cHuI));
    std::copy(std::begin(cHudRp), std::end(cHudRp), std::begin(cHudR));
    std::copy(std::begin(cHudIp), std::end(cHudIp), std::begin(cHudI));







    //Class 8 (WC6)
    for (n = 0; n < DWC6R; n ++) {
        p = WC6R_indices[n][0];
        r = WC6R_indices[n][1];
        s = WC6R_indices[n][2];
        t = WC6R_indices[n][3];
        z = gslpp::complex(0., 0.);
        w = gslpp::complex(0., 0.);
        x = gslpp::complex(0., 0.);
        y = gslpp::complex(0., 0.);
        h = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                for (c = 0; c < 3; c ++) {
                    for (d = 0; d < 3; d ++) {
                        z += Rudag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC6R(cuuR, a, b, c, d), WC6I(cuuI, a, b, c, d))
                                * Ru(b, r) * Ru(d, t);
                        w += Rddag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC6R(cddR, a, b, c, d), WC6I(cddI, a, b, c, d))
                                * Rd(b, r) * Rd(d, t);
                        x += Rldag(p, a) * Rldag(s, c) *
                                gslpp::complex(
                                WC6R(cllR, a, b, c, d), WC6I(cllI, a, b, c, d))
                                * Rl(b, r) * Rl(d, t);
                        y += Rqdag(p, a) * Rqdag(s, c) *
                                gslpp::complex(
                                WC6R(cqq1R, a, b, c, d), WC6I(cqq1I, a, b, c, d))
                                * Rq(b, r) * Rq(d, t);
                        h += Rqdag(p, a) * Rqdag(s, c) *
                                gslpp::complex(
                                WC6R(cqq3R, a, b, c, d), WC6I(cqq3I, a, b, c, d))
                                * Rq(b, r) * Rq(d, t);

                    }
                }
            }
        }
        WC6R_set(cuuRp, p, r, s, t, z.real());
        WC6I_set(cuuIp, p, r, s, t, z.imag());
        WC6R_set(cddRp, p, r, s, t, w.real());
        WC6I_set(cddIp, p, r, s, t, w.imag());

        WC6R_set(cllRp, p, r, s, t, x.real());
        WC6I_set(cllIp, p, r, s, t, x.imag());
        WC6R_set(cqq1Rp, p, r, s, t, y.real());
        WC6I_set(cqq1Ip, p, r, s, t, y.imag());
        WC6R_set(cqq3Rp, p, r, s, t, h.real());
        WC6I_set(cqq3Ip, p, r, s, t, h.imag());
    }


    std::copy(std::begin(cuuRp), std::end(cuuRp), std::begin(cuuR));
    std::copy(std::begin(cuuIp), std::end(cuuIp), std::begin(cuuI));
    std::copy(std::begin(cddRp), std::end(cddRp), std::begin(cddR));
    std::copy(std::begin(cddIp), std::end(cddIp), std::begin(cddI));

    std::copy(std::begin(cllRp), std::end(cllRp), std::begin(cllR));
    std::copy(std::begin(cllIp), std::end(cllIp), std::begin(cllI));
    std::copy(std::begin(cqq1Rp), std::end(cqq1Rp), std::begin(cqq1R));
    std::copy(std::begin(cqq1Ip), std::end(cqq1Ip), std::begin(cqq1I));
    std::copy(std::begin(cqq3Rp), std::end(cqq3Rp), std::begin(cqq3R));
    std::copy(std::begin(cqq3Ip), std::end(cqq3Ip), std::begin(cqq3I));

    //Class 8 (WC7)
    for (n = 0; n < DWC7R; n ++) {
        p = WC7R_indices[n][0];
        r = WC7R_indices[n][1];
        s = WC7R_indices[n][2];
        t = WC7R_indices[n][3];

        z = gslpp::complex(0., 0.);
        w = gslpp::complex(0., 0.);
        x = gslpp::complex(0., 0.);
        y = gslpp::complex(0., 0.);
        h = gslpp::complex(0., 0.);
        f = gslpp::complex(0., 0.);

        le = gslpp::complex(0., 0.);
        lu = gslpp::complex(0., 0.);
        ld = gslpp::complex(0., 0.);
        qe = gslpp::complex(0., 0.);

        qu1 = gslpp::complex(0., 0.);
        qu8 = gslpp::complex(0., 0.);

        qd1 = gslpp::complex(0., 0.);
        qd8 = gslpp::complex(0., 0.);

        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                for (c = 0; c < 3; c ++) {
                    for (d = 0; d < 3; d ++) {
                        z += Rldag(p, a) * Rqdag(s, c) *
                                gslpp::complex(
                                WC7R(clq1R, a, b, c, d), WC7I(clq1I, a, b, c, d))
                                * Rl(b, r) * Rq(d, t);
                        w += Rldag(p, a) * Rqdag(s, c) *
                                gslpp::complex(
                                WC7R(clq3R, a, b, c, d), WC7I(clq3I, a, b, c, d))
                                * Rl(b, r) * Rq(d, t);

                        x += Redag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC7R(ceuR, a, b, c, d), WC7I(ceuI, a, b, c, d))
                                * Re(b, r) * Ru(d, t);
                        y += Redag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cedR, a, b, c, d), WC7I(cedI, a, b, c, d))
                                * Re(b, r) * Rd(d, t);

                        h += Rudag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cud1R, a, b, c, d), WC7I(cud1I, a, b, c, d))
                                * Ru(b, r) * Rd(d, t);
                        f += Rudag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cud8R, a, b, c, d), WC7I(cud8I, a, b, c, d))
                                * Ru(b, r) * Rd(d, t);

                        le += Rldag(p, a) * Redag(s, c) *
                                gslpp::complex(
                                WC7R(cleR, a, b, c, d), WC7I(cleI, a, b, c, d))
                                * Rl(b, r) * Re(d, t);
                        lu += Rldag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC7R(cluR, a, b, c, d), WC7I(cluI, a, b, c, d))
                                * Rl(b, r) * Ru(d, t);
                        ld += Rldag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cldR, a, b, c, d), WC7I(cldI, a, b, c, d))
                                * Rl(b, r) * Rd(d, t);
                        qe += Rqdag(p, a) * Redag(s, c) *
                                gslpp::complex(
                                WC7R(cqeR, a, b, c, d), WC7I(cqeI, a, b, c, d))
                                * Rq(b, r) * Re(d, t);

                        qu1 += Rqdag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC7R(cqu1R, a, b, c, d), WC7I(cqu1I, a, b, c, d))
                                * Rq(b, r) * Ru(d, t);
                        qu8 += Rqdag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC7R(cqu8R, a, b, c, d), WC7I(cqu8I, a, b, c, d))
                                * Rq(b, r) * Ru(d, t);

                        qd1 += Rqdag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cqd1R, a, b, c, d), WC7I(cqd1I, a, b, c, d))
                                * Rq(b, r) * Rd(d, t);
                        qd8 += Rqdag(p, a) * Rddag(s, c) *
                                gslpp::complex(
                                WC7R(cqd8R, a, b, c, d), WC7I(cqd8I, a, b, c, d))
                                * Rq(b, r) * Rd(d, t);

                    }
                }
            }
        }
        WC7R_set(clq1Rp, p, r, s, t, z.real());
        WC7I_set(clq1Ip, p, r, s, t, z.imag());
        WC7R_set(clq3Rp, p, r, s, t, w.real());
        WC7I_set(clq3Ip, p, r, s, t, w.imag());

        WC7R_set(ceuRp, p, r, s, t, x.real());
        WC7I_set(ceuIp, p, r, s, t, x.imag());
        WC7R_set(cedRp, p, r, s, t, y.real());
        WC7I_set(cedIp, p, r, s, t, y.imag());

        WC7R_set(cud1Rp, p, r, s, t, h.real());
        WC7I_set(cud1Ip, p, r, s, t, h.imag());
        WC7R_set(cud8Rp, p, r, s, t, f.real());
        WC7I_set(cud8Ip, p, r, s, t, f.imag());

        WC7R_set(cleRp, p, r, s, t, le.real());
        WC7I_set(cleIp, p, r, s, t, le.imag());
        WC7R_set(cluRp, p, r, s, t, lu.real());
        WC7I_set(cluIp, p, r, s, t, lu.imag());
        WC7R_set(cldRp, p, r, s, t, ld.real());
        WC7I_set(cldIp, p, r, s, t, ld.imag());
        WC7R_set(cqeRp, p, r, s, t, qe.real());
        WC7I_set(cqeIp, p, r, s, t, qe.imag());

        WC7R_set(cqu1Rp, p, r, s, t, qu1.real());
        WC7I_set(cqu1Ip, p, r, s, t, qu1.imag());
        WC7R_set(cqu8Rp, p, r, s, t, qu8.real());
        WC7I_set(cqu8Ip, p, r, s, t, qu8.imag());
        WC7R_set(cqd1Rp, p, r, s, t, qd1.real());
        WC7I_set(cqd1Ip, p, r, s, t, qd1.imag());
        WC7R_set(cqd8Rp, p, r, s, t, qd8.real());
        WC7I_set(cqd8Ip, p, r, s, t, qd8.imag());

    }


    std::copy(std::begin(clq1Rp), std::end(clq1Rp), std::begin(clq1R));
    std::copy(std::begin(clq1Ip), std::end(clq1Ip), std::begin(clq1I));
    std::copy(std::begin(clq3Rp), std::end(clq3Rp), std::begin(clq3R));
    std::copy(std::begin(clq3Ip), std::end(clq3Ip), std::begin(clq3I));
    std::copy(std::begin(ceuRp), std::end(ceuRp), std::begin(ceuR));
    std::copy(std::begin(ceuIp), std::end(ceuIp), std::begin(ceuI));
    std::copy(std::begin(cedRp), std::end(cedRp), std::begin(cedR));
    std::copy(std::begin(cedIp), std::end(cedIp), std::begin(cedI));
    std::copy(std::begin(cud1Rp), std::end(cud1Rp), std::begin(cud1R));
    std::copy(std::begin(cud1Ip), std::end(cud1Ip), std::begin(cud1I));
    std::copy(std::begin(cud8Rp), std::end(cud8Rp), std::begin(cud8R));
    std::copy(std::begin(cud8Ip), std::end(cud8Ip), std::begin(cud8I));

    std::copy(std::begin(cleRp), std::end(cleRp), std::begin(cleR));
    std::copy(std::begin(cleIp), std::end(cleIp), std::begin(cleI));
    std::copy(std::begin(cldRp), std::end(cldRp), std::begin(cldR));
    std::copy(std::begin(cldIp), std::end(cldIp), std::begin(cldI));
    std::copy(std::begin(cluRp), std::end(cluRp), std::begin(cluR));
    std::copy(std::begin(cluIp), std::end(cluIp), std::begin(cluI));
    std::copy(std::begin(cqeRp), std::end(cqeRp), std::begin(cqeR));
    std::copy(std::begin(cqeIp), std::end(cqeIp), std::begin(cqeI));

    std::copy(std::begin(cqu1Rp), std::end(cqu1Rp), std::begin(cqu1R));
    std::copy(std::begin(cqu1Ip), std::end(cqu1Ip), std::begin(cqu1I));
    std::copy(std::begin(cqu8Rp), std::end(cqu8Rp), std::begin(cqu8R));
    std::copy(std::begin(cqu8Ip), std::end(cqu8Ip), std::begin(cqu8I));
    std::copy(std::begin(cqd1Rp), std::end(cqd1Rp), std::begin(cqd1R));
    std::copy(std::begin(cqd1Ip), std::end(cqd1Ip), std::begin(cqd1I));
    std::copy(std::begin(cqd8Rp), std::end(cqd8Rp), std::begin(cqd8R));
    std::copy(std::begin(cqd8Ip), std::end(cqd8Ip), std::begin(cqd8I));

    //Class 8 (WC8)
    for (n = 0; n < DWC8R; n ++) {
        p = WC8R_indices[n][0];
        r = WC8R_indices[n][1];
        s = WC8R_indices[n][2];
        t = WC8R_indices[n][3];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                for (c = 0; c < 3; c ++) {
                    for (d = 0; d < 3; d ++) {
                        z += Redag(p, a) * Redag(s, c) *
                                gslpp::complex(
                                WC8R(ceeR, a, b, c, d), WC8I(ceeI, a, b, c, d))
                                * Re(b, r) * Re(d, t);
                    }
                }
            }
        }
        WC8R_set(ceeRp, p, r, s, t, z.real());
        WC8I_set(ceeIp, p, r, s, t, z.imag());
    }

    std::copy(std::begin(ceeRp), std::end(ceeRp), std::begin(ceeR));
    std::copy(std::begin(ceeIp), std::end(ceeIp), std::begin(ceeI));








    //Class 8 (WC5)
    for (p = 0; p < 3; p ++) {
        for (r = 0; r < 3; r ++) {
            for (s = 0; s < 3; s ++) {
                for (t = 0; t < 3; t ++) {
                    z = gslpp::complex(0., 0.);
                    x = gslpp::complex(0., 0.);
                    y = gslpp::complex(0., 0.);
                    h = gslpp::complex(0., 0.);
                    w = gslpp::complex(0., 0.);


                    for (a = 0; a < 3; a ++) {
                        for (b = 0; b < 3; b ++) {
                            for (c = 0; c < 3; c ++) {
                                for (d = 0; d < 3; d ++) {
                                    h += Rldag(p, a) * Rddag(s, c) *
                                            gslpp::complex(
                                            WC5(cledqR, a, b, c, d), WC5(cledqI, a, b, c, d))
                                            * Re(b, r) * Rq(d, t);

                                    z += Rldag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC5(clequ1R, a, b, c, d), WC5(clequ1I, a, b, c, d))
                                            * Re(b, r) * Ru(d, t);
                                    w += Rldag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC5(clequ3R, a, b, c, d), WC5(clequ3I, a, b, c, d))
                                            * Re(b, r) * Ru(d, t);

                                    x += Rqdag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC5(cquqd1R, a, b, c, d), WC5(cquqd1I, a, b, c, d))
                                            * Ru(b, r) * Rd(d, t);
                                    y += Rqdag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC5(cquqd8R, a, b, c, d), WC5(cquqd8I, a, b, c, d))
                                            * Ru(b, r) * Rd(d, t);



                                }
                            }
                        }
                    }

                    WC5_set(cledqRp, p, r, s, t, h.real());
                    WC5_set(cledqIp, p, r, s, t, h.imag());

                    WC5_set(clequ1Rp, p, r, s, t, z.real());
                    WC5_set(clequ1Ip, p, r, s, t, z.imag());
                    WC5_set(clequ3Rp, p, r, s, t, w.real());
                    WC5_set(clequ3Ip, p, r, s, t, w.imag());

                    WC5_set(cquqd1Rp, p, r, s, t, x.real());
                    WC5_set(cquqd1Ip, p, r, s, t, x.imag());
                    WC5_set(cquqd8Rp, p, r, s, t, y.real());
                    WC5_set(cquqd8Ip, p, r, s, t, y.imag());


                }
            }
        }
    }


    std::copy(std::begin(cledqRp), std::end(cledqRp), std::begin(cledqR));
    std::copy(std::begin(cledqIp), std::end(cledqIp), std::begin(cledqI));

    std::copy(std::begin(clequ1Rp), std::end(clequ1Rp), std::begin(clequ1R));
    std::copy(std::begin(clequ1Ip), std::end(clequ1Ip), std::begin(clequ1I));
    std::copy(std::begin(clequ3Rp), std::end(clequ3Rp), std::begin(clequ3R));
    std::copy(std::begin(clequ3Ip), std::end(clequ3Ip), std::begin(clequ3I));

    std::copy(std::begin(cquqd1Rp), std::end(cquqd1Rp), std::begin(cquqd1R));
    std::copy(std::begin(cquqd1Ip), std::end(cquqd1Ip), std::begin(cquqd1I));
    std::copy(std::begin(cquqd8Rp), std::end(cquqd8Rp), std::begin(cquqd8R));
    std::copy(std::begin(cquqd8Ip), std::end(cquqd8Ip), std::begin(cquqd8I));



}

void RGESolver::EvolveToBasis(
        std::string method, double muI, double muF,
        std::string basis) {
    Evolve(method, muI, muF);
    GoToBasis(basis);

}



