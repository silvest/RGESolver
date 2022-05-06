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

    gslpp::vector<double> yudiag(3), yddiag(3), yediag(3);
    double sqrt2ov = sqrt(2.)/v;

    yudiag(0) = sqrt2ov * mu;
    yudiag(1) = sqrt2ov * mc;
    yudiag(2) = sqrt2ov * mt;

    yddiag(0) = sqrt2ov * md;
    yddiag(1) = sqrt2ov * ms;
    yddiag(2) = sqrt2ov * mb;

    yediag(0) = sqrt2ov * mel;
    yediag(1) = sqrt2ov * mmu;
    yediag(2) = sqrt2ov * mtau;

    yuR = gslpp::matrix<double>(yudiag);
    yuI.reset();
    yeR = gslpp::matrix<double>(yediag);
    yeI.reset();
    ydR = gslpp::matrix<double>(yddiag);
    ydI.reset();

    if (basis == "UP") {

        gslpp::matrix<gslpp::complex> yd = ydR * CKM.hconjugate(); 
        ydR = yd.real();
        ydI = yd.imag();

    }
    else if (basis == "DOWN") {

        gslpp::matrix<gslpp::complex> yu = yuR * CKM; 
        yuR = yu.real();
        yuI = yu.imag();

    }
    else 
        std::cout << "WARNING: wrong basis choice: " << basis << ": Yukawa couplings not updated!!!" << std::endl;

}

void RGESolver::UpdateCKM() {
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
                int status = gsl_odeiv2_evolve_apply(evo_, con_, s_, &sys_, &ttmp, tF, &step_, x);
                if (status != GSL_SUCCESS) {
                    printf("error in adaptive integration, return value=%d\n", status);
                }
            }

        } else {

            while (ttmp > tF) {
                int status = gsl_odeiv2_evolve_apply(evo_, con_, s_, &sys_, &ttmp, tF, &step_, x);
                if (status != GSL_SUCCESS) {
                    printf("error in adaptive integration, return value=%d\n", status);
                }

            }
        }


        gsl_odeiv2_evolve_reset(evo_);
        gsl_odeiv2_step_reset(s_);

    }

    //Leading-log resolution 
    //-------------------------------------------
    if (method == "Leading-Log") {
        double beta[dim] = {0.};
        double Log_muF_over_muI = log(muF / muI);
        /*int status = */func(10., x, beta, NULL);
        for (int i = 0; i < dim; i ++) {
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

    gslpp::matrix<gslpp::complex> yu(3, 3, 0.);
    gslpp::matrix<gslpp::complex> yd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> ye(3, 3, 0.);

    yu = yuR + gslpp::complex::i() * yuI;
    yd = ydR + gslpp::complex::i() * ydI;
    ye = yeR + gslpp::complex::i() * yeI;
    
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
    gslpp::matrix<gslpp::complex> Rd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rq(3, 3, 0.);

    //Computing the CKM 
    CKM = (Vu.hconjugate()) * Vd;
    //Extract the 4 parameters from the raw CKM
    ExtractParametersFromCKM();
    //Build the CKM with the 4 parameters
    UpdateCKM();

    yuR = gslpp::matrix<double>(Su);
    yuI.reset();
    yeR = gslpp::matrix<double>(Se);
    yeI.reset();
    ydR = gslpp::matrix<double>(Sd);
    ydI.reset();

    if (basis == "UP") {

        Ru = Uu;
        Rd = Ud;
        Rq = Vu;

        gslpp::matrix<gslpp::complex> yd = ydR * CKM.hconjugate(); 
        ydR = yd.real();
        ydI = yd.imag();

    }
    else if (basis == "DOWN") {

        Ru = Uu;
        Rd = Ud;
        Rq = Vd;

        gslpp::matrix<gslpp::complex> yu = yuR * CKM; 
        yuR = yu.real();
        yuI = yu.imag();

    }
    else 
        std::cout << "WARNING: wrong basis choice: " << basis << ": Yukawa couplings not updated!!!" << std::endl;


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
            x[n + count * DF] = yuR(i,j);
            count ++;
            x[n + count * DF] = yuI(i,j);
            count ++;
            x[n + count * DF] = ydR(i,j);
            count ++;
            x[n + count * DF] = ydI(i,j);
            count ++;
            x[n + count * DF] = yeR(i,j);
            count ++;
            x[n + count * DF] = yeI(i,j);
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
            yuR.assign(i,j,x[n + a * DF]);
            a ++;
            yuI.assign(i,j,x[n + a * DF]);
            a ++;
            ydR.assign(i,j,x[n + a * DF]);
            a ++;
            ydI.assign(i,j,x[n + a * DF]);
            a ++;
            yeR.assign(i,j,x[n + a * DF]);
            a ++;
            yeI.assign(i,j,x[n + a * DF]);
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

    gslpp::matrix<double> yeR(NG,0.), yeI(NG,0.), //yukawas
            ydR(NG,0.), ydI(NG,0.),
            yuR(NG,0.), yuI(NG,0.);
    gslpp::matrix<double> yedagR(NG,0.), yedagI(NG,0.),
            yddagR(NG,0.), yddagI(NG,0.),
            yudagR(NG,0.), yudagI(NG,0.);

    for (i = 0; i < NG; i ++) {
        for (j = 0; j < NG; j ++) {
            a = 0;
            yuR.assign(i,j,y[c]);
            yudagR.assign(j,i,y[c]);
            a ++;
            yuI.assign(i,j,y[c + DF]);
            yudagI.assign(j,i,- y[c + DF]);
            a ++;
            ydR.assign(i,j,y[c + 2 * DF]);
            yddagR.assign(j,i,y[c + 2 * DF]);
            a ++;
            ydI.assign(i,j,y[c + 3 * DF]);
            yddagI.assign(j,i,- y[c + 3 * DF]);
            a ++;
            yeR.assign(i,j,y[c + 4 * DF]);
            yedagR.assign(j,i,y[c + 4 * DF]);
            a ++;
            yeI.assign(i,j,y[c + 5 * DF]);
            yedagI.assign(j,i,- y[c + 5 * DF]);
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
    gslpp::matrix<double> yudyuR(NG,0.),yudyuI(NG,0.) ; //yu^dag yu 
    gslpp::matrix<double> yddydR(NG,0.),yddydI(NG,0.); //yd^dag yd 
    gslpp::matrix<double> yedyeR(NG,0.),yedyeI(NG,0.); //ye^dag ye 
    gslpp::matrix<double> ydyudR(NG,0.),ydyudI(NG,0.); //yd yu^dag

    yudyuR = yudagR * yuR - yudagI * yuI;
    yudyuI = yudagI * yuR + yudagR * yuI;
    yddydR = yddagR * ydR - yddagI * ydI;
    yddydI = yddagI * ydR + yddagR * ydI;
    yedyeR = yedagR * yeR - yedagI * yeI;
    yedyeI = yedagI * yeR + yedagR * yeI;
    ydyudR = ydR * yudagR - ydI * yudagI;
    ydyudI = ydI * yudagR + ydR * yudagI;
  




    //gammaH
    for (i = 0; i < NG; i ++) {
        gammaH += yedyeR(i,i) + NC * (yudyuR(i,i) + yddydR(i,i));
        for (j = 0; j < NG; j ++) {
            H += yedyeR(i,j) * yedyeR(j,i) - yedyeI(i,j) * yedyeI(j,i)
                    + NC * (yudyuR(i,j) * yudyuR(j,i) - yudyuI(i,j) * yudyuI(j,i)
                    + yddydR(i,j) * yddydR(j,i) - yddydI(i,j) * yddydI(j,i));
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
            f[c ] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yuR(i,j)
                    ;
            //yuI
            f[c + DF] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yuI(i,j)
                    ;
            //ydR     
            f[c + 2 * DF] = (gammaH - (5. / 12.) * g12 - 2.25 * g22
                    - 8. * g32) * ydR(i,j)
                    ;
            //ydI  
            f[c + 3 * DF] = (gammaH - (5. / 12.) * g12 - 2.25 * g22
                    - 8. * g32) * ydI(i,j)
                    ;
            //yeR
            f[c + 4 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * yeR(i,j)
                    ;
            //yeI
            f[c + 5 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * yeI(i,j)
                    ;
            //Entries with 1 matrix product (1 summed index)
            for (b = 0; b < NG; b ++) {
                //yuR
                f[c ] +=
                        + 1.5 * (yuR(i,b)*(yudyuR(b,j) - yddydR(b,j))
                        - yuI(i,b)*(yudyuI(b,j) - yddydI(b,j)))
                        ;
                //yuI
                f[c + DF] += 1.5 * (yuI(i,b)*(yudyuR(b,j) - yddydR(b,j))
                        + yuR(i,b)*(yudyuI(b,j) - yddydI(b,j)))
                        ;
                //ydR       
                f[c + 2 * DF] += 1.5 * (ydR(i,b)*(yddydR(b,j) - yudyuR(b,j))
                        - ydI(i,b)*(yddydI(b,j) - yudyuI(b,j)))
                        ;
                //ydI 
                f[c + 3 * DF] += 1.5 * (
                        ydI(i,b)*(yddydR(b,j) - yudyuR(b,j))
                        + ydR(i,b)*(yddydI(b,j) - yudyuI(b,j)))
                        ;
                //yeR
                f[c + 4 * DF] += 1.5 * (yeR(i,b)*(yedyeR(b,j))
                        - yeI(i,b)*(yedyeI(b,j)))
                        ;
                //yeI 
                f[c + 5 * DF] += 1.5 * (yeI(i,b)*(yedyeR(b,j))
                        + yeR(i,b)*(yedyeI(b,j)))
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
    CKM_theta12 = 0.2;
    CKM_theta13 = 0.1;
    CKM_theta23 = 0.3;
    CKM_delta = 3.14 / 4.;
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
    g3 = 1.2;
    g2 = 0.6;
    g1 = .31;

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

    yuR.reset();
    yuI.reset();
    ydR.reset();
    ydI.reset();
    yeR.reset();
    yeI.reset();
    
//    for (i = 0; i < NG; i ++) {
//        for (j = 0; j < NG; j ++) {
//            a = 0;
//            yuR.assign(i,j,0);
//            a ++;
//            yuI.assign(i,j,0);
//            a ++;
//            ydR.assign(i,j,0);
//            a ++;
//            ydI.assign(i,j,0);
//            a ++;
//            yeR.assign(i,j,0);
//            a ++;
//            yeI.assign(i,j,0);
//            a ++;
//            n ++;
//        }
//    }
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
            yu.assign(i, j, gslpp::complex(yuR(i,j), yuI(i,j), false));
            yd.assign(i, j, gslpp::complex(ydR(i,j), ydI(i,j), false));
            ye.assign(i, j, gslpp::complex(yeR(i,j), yeI(i,j), false));
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

    yuR = gslpp::matrix<double>(Su);
    yuI.reset();
    yeR = gslpp::matrix<double>(Se);
    yeI.reset();
    ydR = gslpp::matrix<double>(Sd);
    ydI.reset();

    if (basis == "UP") {

        Ru = Uu;
        Rd = Ud;
        Rq = Vu;

        gslpp::matrix<gslpp::complex> yd = ydR * CKM.hconjugate(); 
        ydR = yd.real();
        ydI = yd.imag();

    }
    else if (basis == "DOWN") {

        Ru = Uu;
        Rd = Ud;
        Rq = Vd;

        gslpp::matrix<gslpp::complex> yu = yuR * CKM; 
        yuR = yu.real();
        yuI = yu.imag();

    }
    else 
        std::cout << "WARNING: wrong basis choice: " << basis << ": Yukawa couplings not updated!!!" << std::endl;

    //I save the hermitian conjugates of the rotation 
    //matrix to be more efficient. 
    Rudag = Ru.hconjugate();
    Rddag = Rd.hconjugate();
    Rqdag = Rq.hconjugate();

    Redag = Re.hconjugate();
    Rldag = Rl.hconjugate();





    int a, b, c, d, p, r, s, t, n;
    gslpp::complex z;



    //Coefficients in the new basis
    double cuHRp[3 * 3] = {0.};
    double cuHIp[3 * 3] = {0.};

    double cHl1Rp[3 * 3] = {0.};
    double cHl1Ip[3 * 3] = {0.};

    double cuuRp[81] = {0.};
    double cuuIp[81] = {0.};

    double clequ1Rp[81] = {0.};
    double clequ1Ip[81] = {0.};



    //Class WC1
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

    //The coefficients are updated with the rotated ones
    std::copy(std::begin(cuHRp), std::end(cuHRp), std::begin(cuHR));
    std::copy(std::begin(cuHIp), std::end(cuHIp), std::begin(cuHI));


    //Class WC2
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

    std::copy(std::begin(cHl1Rp), std::end(cHl1Rp), std::begin(cHl1R));
    std::copy(std::begin(cHl1Ip), std::end(cHl1Ip), std::begin(cHl1I));


    //Class WC6
    for (n = 0; n < DWC6R; n ++) {
        p = WC6R_indices[n][0];
        r = WC6R_indices[n][1];
        s = WC6R_indices[n][2];
        t = WC6R_indices[n][3];
        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                for (c = 0; c < 3; c ++) {
                    for (d = 0; d < 3; d ++) {
                        z += Rudag(p, a) * Rudag(s, c) *
                                gslpp::complex(
                                WC6R(cuuR, a, b, c, d), WC6I(cuuI, a, b, c, d))
                                * Ru(b, r) * Ru(d, t);
                    }
                }
            }
        }
        WC6R_set(cuuRp, p, r, s, t, z.real());
        WC6I_set(cuuIp, p, r, s, t, z.imag());
    }

    std::copy(std::begin(cuuRp), std::end(cuuRp), std::begin(cuuR));
    std::copy(std::begin(cuuIp), std::end(cuuIp), std::begin(cuuI));

    //Class WC5
    for (p = 0; p < 3; p ++) {
        for (r = 0; r < 3; r ++) {
            for (s = 0; s < 3; s ++) {
                for (t = 0; t < 3; t ++) {
                    z = gslpp::complex(0., 0.);
                    for (a = 0; a < 3; a ++) {
                        for (b = 0; b < 3; b ++) {
                            for (c = 0; c < 3; c ++) {
                                for (d = 0; d < 3; d ++) {
                                    z += Rldag(p, a) * Rqdag(s, c) *
                                            gslpp::complex(
                                            WC6R(cuuR, a, b, c, d), WC6I(cuuI, a, b, c, d))
                                            * Re(b, r) * Ru(d, t);
                                }
                            }
                        }
                    }
                    WC5_set(clequ1Rp, p, r, s, t, z.real());
                    WC5_set(clequ1Ip, p, r, s, t, z.imag());
                }
            }
        }
    }

    std::copy(std::begin(clequ1Rp), std::end(clequ1Rp), std::begin(clequ1R));
    std::copy(std::begin(clequ1Ip), std::end(clequ1Ip), std::begin(clequ1I));






}

void RGESolver::EvolveToBasis(
        std::string method, double muI, double muF,
        std::string basis) {
    Evolve(method, muI, muF);
    GoToBasis(basis);

}




/*
void RGESolver::TEST() {

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


    using namespace std;


    //Matrix to rotate fields
    gslpp::matrix<gslpp::complex> Ru(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rudag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rddag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rq(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rqdag(3, 3, 0.);





    double cuHRp[3 * 3] = {0.};
    double cuHIp[3 * 3] = {0.};

    double cHl1Rp[3 * 3] = {0.};
    double cHl1Ip[3 * 3] = {0.};


    WC2R_set(cHl1R, 0, 0, 0.4);
    WC2R_set(cHl1R, 0, 2, 0.5);
    WC2I_set(cHl1I, 0, 1, -0.2);


    cout << "Before : " << endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            cout << WC2R(cHl1R, i, j) << " + i * "
                    << WC2I(cHl1I, i, j) << "\t";
        }
        cout << "" << endl;
    }


    Ru.assign(0, 1, gslpp::complex(1., 2.));
    Ru.assign(2, 1, gslpp::complex(1., 0.));
    Ru.assign(0, 0, 2.);
    Ru.assign(2, 0, 1);
    Ru.assign(2, 2, 4);


    Rq.assign(0, 1, gslpp::complex(1., 2.));
    Rq.assign(2, 1, gslpp::complex(1., 0.));
    Rq.assign(0, 0, gslpp::complex(2., 0.));
    Rq.assign(1, 2, gslpp::complex(3., 0.));
    Rq.assign(2, 2, gslpp::complex(1., 0.));


    //I save the hermitian conjugates of the rotation 
    //matrix to be more efficient. 
    Rudag = Ru.hconjugate();
    Rddag = Rd.hconjugate();
    Rqdag = Rq.hconjugate();

    Redag = Re.hconjugate();
    Rldag = Rl.hconjugate();


    int a, b,n;
    gslpp::complex z;

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

    std::copy(std::begin(cuHRp), std::end(cuHRp), std::begin(cuHR));
    std::copy(std::begin(cuHIp), std::end(cuHIp), std::begin(cuHI));

    

    for (n = 0; n < DWC2R; n ++) {
        i = WC2R_indices[n][0];
        j = WC2R_indices[n][1];

        z = gslpp::complex(0., 0.);
        for (a = 0; a < 3; a ++) {
            for (b = 0; b < 3; b ++) {
                z += Rudag(i,a) *
                        gslpp::complex(
                        WC2R(cHl1R, a, b), WC2I(cHl1I, a, b)) *
                        Ru(b,j);
            }
        }
        WC2R_set(cHl1Rp, i, j, z.real());
        WC2I_set(cHl1Ip, i, j, z.imag());
    }

    std::copy(std::begin(cHl1Rp), std::end(cHl1Rp), std::begin(cHl1R));
    std::copy(std::begin(cHl1Ip), std::end(cHl1Ip), std::begin(cHl1I));
    cout << "After : " << endl;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            cout << WC2R(cHl1R, i, j) << " + i * "
                    << WC2I(cHl1I, i, j) << "\t";
        }
        cout << "" << endl;
    }


}
 */