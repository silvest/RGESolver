#include <gsl/gsl_linalg.h>

#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
//Extracts the 4 parameters from the CKM 
//affected by unphysical phases

void RGESolver::ExtractParametersFromCKM() {

    double s13 = gsl_complex_abs(gsl_matrix_complex_get(CKM, 0, 2));
    CKM_theta13 = asin(s13);
    double c13 = sqrt(1. - s13 * s13);

    double c12 = gsl_complex_abs(gsl_matrix_complex_get(CKM, 0, 0)) / c13;
    CKM_theta12 = acos(c12);
    double s12 = sqrt(1. - c12 * c12);

    double s23 = gsl_complex_abs(gsl_matrix_complex_get(CKM, 1, 2)) / c13;
    double c23 = sqrt(1. - s23 * s23);
    CKM_theta23 = asin(s23);
    double temp = (gsl_complex_abs(gsl_matrix_complex_get(CKM, 1, 0)) *
            gsl_complex_abs(gsl_matrix_complex_get(CKM, 1, 0)) - s12 * s12 * c23 * c23
            - c12 * c12 * s23 * s23 * s13 * s13) /
            (2. * s12 * c23 * c12 * s23 * s13);
    CKM_delta = acos((temp));
    //Debugging:
    std::cout << "angles and phase : "
            << CKM_theta12 << "\t"
            << CKM_theta13 << "\t"
            << CKM_theta23 << "\t"
            << CKM_delta << "\t" << std::endl;

}


//Builds the CKM from the 4 parameters

void RGESolver::UpdateCKM() {
    using namespace std;
    gsl_complex c12 = gsl_complex_rect(cos(CKM_theta12), 0.);
    gsl_complex s12 = gsl_complex_rect(sin(CKM_theta12), 0.);
    gsl_complex c13 = gsl_complex_rect(cos(CKM_theta13), 0.);
    gsl_complex s13 = gsl_complex_rect(sin(CKM_theta13), 0.);
    gsl_complex c23 = gsl_complex_rect(cos(CKM_theta23), 0.);
    gsl_complex s23 = gsl_complex_rect(sin(CKM_theta23), 0.);
    gsl_complex PhaseFactor;
    GSL_REAL(PhaseFactor) = cos(CKM_delta);
    GSL_IMAG(PhaseFactor) = sin(CKM_delta);
    cout << GSL_REAL(PhaseFactor) <<
            " +I*" << GSL_IMAG(PhaseFactor) << endl;

    gsl_complex a, b;

    //First row 
    gsl_matrix_complex_set(CKM, 0, 0,
            gsl_complex_mul(c12, c13));
    gsl_matrix_complex_set(CKM, 0, 1,
            gsl_complex_mul(s12, c13));
    gsl_matrix_complex_set(CKM, 0, 2,
            gsl_complex_mul(s13, gsl_complex_conjugate(PhaseFactor)));



    //Second row
    a = gsl_complex_negative(gsl_complex_mul(s12, c23));
    b = gsl_complex_negative(gsl_complex_mul(gsl_complex_mul(c12, s23), gsl_complex_mul(s13, PhaseFactor)));
    gsl_matrix_complex_set(CKM, 1, 0,
            gsl_complex_add(a, b));

    a = gsl_complex_mul(c12, c23);
    b = gsl_complex_negative(gsl_complex_mul(gsl_complex_mul(s12, s23), gsl_complex_mul(s13, PhaseFactor)));
    gsl_matrix_complex_set(CKM, 1, 1,
            gsl_complex_add(a, b));

    gsl_matrix_complex_set(CKM, 1, 2, gsl_complex_mul(s23, c13));


    //Third row
    a = gsl_complex_mul(s12, s23);
    b = gsl_complex_negative(gsl_complex_mul(gsl_complex_mul(c12, c23), gsl_complex_mul(s13, PhaseFactor)));
    gsl_matrix_complex_set(CKM, 2, 0,
            gsl_complex_add(a, b));

    a = gsl_complex_negative(gsl_complex_mul(c12, s23));
    b = gsl_complex_negative(gsl_complex_mul(gsl_complex_mul(s12, c23), gsl_complex_mul(s13, PhaseFactor)));
    gsl_matrix_complex_set(CKM, 2, 1,
            gsl_complex_add(a, b));

    gsl_matrix_complex_set(CKM, 2, 2, gsl_complex_mul(c23, c13));

}

void RGESolver::GenerateSMInitialConditions(
        double mu, std::string basis) {
    std::cout << "muSM, mu : " << InputScale_SM << "  " << mu << std::endl;
    EvolveSMOnly("Numeric", InputScale_SM, mu);
    ExtractParametersFromCKM();
    UpdateCKM();
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

void RGESolver::GoToBasis(std::string basis) {
    gsl_matrix_complex *A = gsl_matrix_complex_alloc(3, 3);
    gsl_matrix_complex *V = gsl_matrix_complex_alloc(3, 3);
    gsl_vector_complex *S = gsl_vector_complex_alloc(3);

    gsl_vector *work = gsl_vector_alloc(3);

    srand(time(0));
    int i, j;
    double a, b;
    gsl_complex z;
    for (i = 0; i < 3; i ++) {
        for (j = 0; j < 3; j ++) {
            a = rand();
            b = rand();
            z = gsl_complex_rect(a, b);
            gsl_matrix_complex_set(A, i, j, z);
        }
    }
    //gsl_linalg_complex_SV_decomp(A, V, S, work);
    
}