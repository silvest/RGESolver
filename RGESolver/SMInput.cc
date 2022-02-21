#include <gsl/gsl_linalg.h>

#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
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

 
    //First method 
    //delta in [0,Pi]
    double temp = (CKM(1, 0).abs() * CKM(1, 0).abs() - s12 * s12 * c23 * c23
            - c12 * c12 * s23 * s23 * s13 * s13) /
            (2. * s12 * c23 * c12 * s23 * s13);
    CKM_delta = (arccos(gslpp::complex(temp, 0.))).real();

    //Second method 
    //delta in [-Pi,Pi]
    gslpp::complex R = (CKM(0, 0) * CKM(1, 1)) /
            (CKM(0, 1) * CKM(1, 0));
    gslpp::complex temp2 = (- c12 * c12 * c23 - c23 * R * s12 * s12) /
            (c12 * (R - 1.) * s12 * s13 * s23);
    CKM_delta = (gslpp::complex(0., - 1.) * log(temp2)).real();
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
        std::string method,
        bool inputCKM) {
    if (method != "Numeric" && method != "Leading-Log") {
        std::cout << "WARNING : invalid method\n"
                "Available methods: Numeric, Leading-Log"
                << std::endl;
    }


    //Before evolving, eventual changes 
    //in the input values of CKM angles must be 
    //translated 
    c12 = cos(CKM_theta12);
    s12 = sin(CKM_theta12);
    c13 = cos(CKM_theta13);
    s13 = sin(CKM_theta13);
    c23 = cos(CKM_theta23);
    s23 = sin(CKM_theta23);


    if (inputCKM == true) {
        UpdateCKM();
        FromMassesToYukawas(basis);
    }

    EvolveSMOnly(method, InputScale_SM, mu);
    GoToBasisSMOnly(basis);
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
    /*double gqR[NG][NG] = {0.}; //wavefunction ren. const. of q 
    double gqI[NG][NG] = {0.};
    double glR[NG][NG] = {0.}; //wavefunction ren. const. of l 
    double glI[NG][NG] = {0.};
    double guR[NG][NG] = {0.}; //wavefunction ren. const. of u
    double guI[NG][NG] = {0.};
    double gdR[NG][NG] = {0.}; //wavefunction ren. const. of d 
    double gdI[NG][NG] = {0.};
    double geR[NG][NG] = {0.}; //wavefunction ren. const. of e 
    double geI[NG][NG] = {0.};*/
    double yudyuR[NG][NG] = {0.}; //yu^dag yu 
    double yudyuI[NG][NG] = {0.};
    double yddydR[NG][NG] = {0.}; //yd^dag yd 
    double yddydI[NG][NG] = {0.};
    double yedyeR[NG][NG] = {0.}; //ye^dag ye 
    double yedyeI[NG][NG] = {0.};
    double ydyudR[NG][NG] = {0.}; //yd yu^dag
    double ydyudI[NG][NG] = {0.};
    //double yuyddR[NG][NG] = {0.}; //yu yd^dag
    //double yuyddI[NG][NG] = {0.};

    //3 Yukawa matrices product
    /*
    double yuyudyuR[NG][NG] = {0.}; //yu yu^dag yu
    double yuyudyuI[NG][NG] = {0.};
    double yudyuyudR[NG][NG] = {0.}; //yu^dag yu yu^dag
    double yudyuyudI[NG][NG] = {0.};
    double ydyddydR[NG][NG] = {0.}; //yd yd^dag yd
    double ydyddydI[NG][NG] = {0.};
    double yddydyddR[NG][NG] = {0.}; //yd^dag yd yd^dag
    double yddydyddI[NG][NG] = {0.};
    double yeyedyeR[NG][NG] = {0.}; //ye ye^dag ye
    double yeyedyeI[NG][NG] = {0.};
    double yedyeyedR[NG][NG] = {0.}; //ye^dag ye ye^dag
    double yedyeyedI[NG][NG] = {0.};

    double yddydyudR[NG][NG] = {0.}; // yd^dag yd yu^dag
    double yddydyudI[NG][NG] = {0.};
    double yudyuyddR[NG][NG] = {0.}; // yu^dag yu yd^dag
    double yudyuyddI[NG][NG] = {0.};*/


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


/*void RGESolver::GoToBasis(std::string basis) {

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

    int i, j, k, l;
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
    yu.singularvalue(Uu, Vu, Su);
    std::cout << Su << std::endl;
    yd.singularvalue(Ud, Vd, Sd);
    std::cout << Sd << std::endl;
    ye.singularvalue(Re, Rl, Se);
    std::cout << Se << std::endl;

    //Matrix to rotate fields
    gslpp::matrix<gslpp::complex> Ru(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rudag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rd(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rddag(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rq(3, 3, 0.);
    gslpp::matrix<gslpp::complex> Rqdag(3, 3, 0.);


    CKM = Vu * (Vd.hconjugate());
    ExtractParametersFromCKM();
    UpdateCKM();

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
        Ru = Uu;
        Rd = Ud;
        Rq = Vu;
        gslpp::matrix<gslpp::complex> CKMdag = CKM.hconjugate();

        for (i = 0; i < 3; i ++) {
            yuR[i][i] = Su(i);
            yeR[i][i] = Se(i);
            for (j = 0; j < 3; j ++) {
                ydR[i][j] = Sd(i)*(CKMdag(i, j)).real();
                ydI[i][j] = Sd(i)*(CKMdag(i, j)).imag();
            }
        }


    }
    if (basis == "DOWN") {
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

    Rqdag = Rq.hconjugate();
    Rudag = Ru.hconjugate();
    Rddag = Rd.hconjugate();

    Rldag = Rl.hconjugate();
    Redag = Re.hconjugate();




    double ceHRtmp[3 * 3];
    double ceHItmp[3 * 3];
    double cuHRtmp[3 * 3];
    double cuHItmp[3 * 3];
    double cdHRtmp[3 * 3];
    double cdHItmp[3 * 3];


    //Class 5 rotations 
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            WC1_set(ceHRtmp, i, j, WC1(ceHR, i, j));
            WC1_set(ceHItmp, i, j, WC1(ceHI, i, j));
            WC1_set(ceHR, i, j, 0.);
            WC1_set(ceHI, i, j, 0.);

            WC1_set(cuHRtmp, i, j, WC1(cuHR, i, j));
            WC1_set(cuHItmp, i, j, WC1(cuHI, i, j));
            WC1_set(cuHR, i, j, 0.);
            WC1_set(cuHI, i, j, 0.);

            WC1_set(cdHRtmp, i, j, WC1(cdHR, i, j));
            WC1_set(cdHItmp, i, j, WC1(cdHI, i, j));
            WC1_set(cdHR, i, j, 0.);
            WC1_set(cdHI, i, j, 0.);

        }
    }

    double o1r, o1i, o2r, o2i, o3r, o3i;



    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            o1r = o1i = o2r = o2i = 0;
            for (k = 0; k < 3; k ++) {
                for (l = 0; l < 3; l ++) {
                    o1r += (Rldag(i, k) * WC1(ceHRtmp, k, l) * Re(l, j)).real();
                    o1i += (Rldag(i, k) * WC1(ceHItmp, k, l) * Re(l, j)).real();
                }
            }
        }
    }


}*/