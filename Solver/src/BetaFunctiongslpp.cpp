

int RGESolver::func(double logmu, const double y[], double f[], void* params) {
    int i, j, a, b, l, k, d, w, count, v;
    int p, r, s, t;
    double loop_factor = 1.
    / (16.*3.1415926535 * 3.1415926535);
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

    gslpp::matrix<gslpp::complex> ye(NG,0.), yd(NG,0.), yu(NG,0.);
    gslpp::matrix<gslpp::complex> yedag(NG,0.), yddag(NG,0.), yudag(NG,0.);
    
    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            a = 0;
            yu.assignre(i, j, y[c]);
            yudag.assignre(j, i, y[c]);
            a++;
            yu.assignim(i, j, y[c + DF]);
            yudag.assignim(j, i, -y[c + DF]);
            a++;
            yd.assignre(i, j, y[c + 2 * DF]);
            yddag.assignre(j, i, y[c + 2 * DF]);
            a++;
            yd.assignim(i, j, y[c + 3 * DF]);
            yddag.assignim(j, i, -y[c + 3 * DF]);
            a++;
            ye.assignre(i, j, y[c + 4 * DF]);
            yedag.assignre(j, i, y[c + 4 * DF]);
            a++;
            ye.assignim(i, j, y[c + 5 * DF]);
            yedag.assignim(j, i, -y[c + 5 * DF]);
            a++;
            c++;
        }
    }
    c += (2. * Nyukawa - 1.) * DF;
    
    
    
    //SMEFT class 1 
    double CG = y[c];
    double CGT = y[c + 1];
    double CW = y[c + 2];
    double CWT = y[c + 3];
    c += N1;
    
    //SMEFT class 2-3 
    double CH = y[c];
    double CHBOX = y[c + 1];
    double CHD = y[c + 2];
    c += N23;
    //SMEFT class 4 
    double CHG = y[c];
    double CHB = y[c + 1];
    double CHW = y[c + 2];
    double CHWB = y[c + 3];
    
    double CHGT = y[c + 4];
    double CHBT = y[c + 5];
    double CHWT = y[c + 6];
    double CHWBT = y[c + 7];
    c += N4;
    //SMEFT class 5
    double CeHR[NG * NG], CeHI[NG * NG], CuHR[NG * NG],
            CuHI[NG * NG], CdHR[NG * NG], CdHI[NG * NG];
    
    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            a = 0;
            WC1_set(CuHR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CuHI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CdHR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CdHI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CeHR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CeHI, i, j, y[c + a * DF]);
            a++;
            c++;
        }
    }
    c += (N5 * 2 - 1) * DF;
    
    
    //SMEFT class 6 
    double CeWR[NG * NG], CeWI[NG * NG];
    double CeBR[NG * NG], CeBI[NG * NG];
    double CuGR[NG * NG], CuGI[NG * NG];
    double CuWR[NG * NG], CuWI[NG * NG], CuBR[NG * NG], CuBI[NG * NG];
    double CdGR[NG * NG], CdGI[NG * NG], CdWR[NG * NG], CdWI[NG * NG],
            CdBR[NG * NG], CdBI[NG * NG];
    for (i = 0; i < NG; i++) {
        for (j = 0; j < NG; j++) {
            a = 0;
            WC1_set(CeWR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CeWI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CeBR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CeBI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CuGR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CuGI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CuWR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CuWI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CuBR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CuBI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CdGR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CdGI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CdWR, i, j, y[c + a * DF]);
            a++;
            WC1_set(CdWI, i, j, y[c + a * DF]);
            a++;
            WC1_set(CdBR, i, j, y[c + a * DF]);
			a++;
			WC1_set(CdBI, i, j, y[c + a * DF]);
			a++;
			c++;
		}
	}
	c += (N6 * 2 - 1) * DF;

	//SMEFT class 7
	double CHl1R[NG * NG], CHl3R[NG * NG], CHl1I[NG * NG], CHl3I[NG * NG];
	double CHeR[NG * NG], CHeI[NG * NG];
	double CHq1R[NG * NG], CHq1I[NG * NG], CHq3R[NG * NG], CHq3I[NG * NG];
	double CHuR[NG * NG], CHuI[NG * NG], CHdR[NG * NG], CHdI[NG * NG];
	double CHudR[NG * NG], CHudI[NG * NG];

	{
		for (i = 0; i < DWC2R; i++) {
			WC2R_set(CHl1R, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2I; i++) {
			WC2I_set(CHl1I, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2R; i++) {
			WC2R_set(CHl3R, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2I; i++) {
			WC2I_set(CHl3I, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2R; i++) {
			WC2R_set(CHeR, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2I; i++) {
			WC2I_set(CHeI, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
			c++;
		}

		for (i = 0; i < DWC2R; i++) {
			WC2R_set(CHq1R, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2I; i++) {
			WC2I_set(CHq1I, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2R; i++) {
			WC2R_set(CHq3R, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2I; i++) {
			WC2I_set(CHq3I, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
			c++;
		}

		for (i = 0; i < DWC2R; i++) {
			WC2R_set(CHuR, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2I; i++) {
			WC2I_set(CHuI, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
			c++;
		}

		for (i = 0; i < DWC2R; i++) {
			WC2R_set(CHdR, WC2R_indices[i][0], WC2R_indices[i][1], y[c]);
			c++;
		}
		for (i = 0; i < DWC2I; i++) {
			WC2I_set(CHdI, WC2I_indices[i][0], WC2I_indices[i][1], y[c]);
			c++;
		}

		for (i = 0; i < NG; i++) {
			for (j = 0; j < NG; j++) {
				WC1_set(CHudR, i, j, y[c]);
				c++;
			}
		}
		for (i = 0; i < NG; i++) {
			for (j = 0; j < NG; j++) {
				WC1_set(CHudI, i, j, y[c]);
				c++;
			}
		}
	}

	//SMEFT class 8_LLLL
	double CllR[NG * NG * NG * NG], CllI[NG * NG * NG * NG], Cqq1R[NG * NG * NG * NG], Cqq1I[NG * NG * NG * NG],
		Cqq3R[NG * NG * NG * NG], Cqq3I[NG * NG * NG * NG];

	double Clq1R[NG * NG * NG * NG];
	double Clq1I[NG * NG * NG * NG];
	double Clq3R[NG * NG * NG * NG];
	double Clq3I[NG * NG * NG * NG];

	for (a = 0; a < DWC6R; a++) {
		WC6R_set(CllR, WC6R_indices[a][0], WC6R_indices[a][1],
						WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC6I; a++) {
		WC6I_set(CllI, WC6I_indices[a][0], WC6I_indices[a][1],
						WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC6R; a++) {
		WC6R_set(Cqq1R, WC6R_indices[a][0], WC6R_indices[a][1],
						WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC6I; a++) {
		WC6I_set(Cqq1I, WC6I_indices[a][0], WC6I_indices[a][1],
						WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC6R; a++) {
		WC6R_set(Cqq3R, WC6R_indices[a][0], WC6R_indices[a][1],
						WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC6I; a++) {
		WC6I_set(Cqq3I, WC6I_indices[a][0], WC6I_indices[a][1],
						WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC7R; a++) {
		WC7R_set(Clq1R, WC7R_indices[a][0], WC7R_indices[a][1],
						WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC7I; a++) {
		WC7I_set(Clq1I, WC7I_indices[a][0], WC7I_indices[a][1],
						WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC7R; a++) {
		WC7R_set(Clq3R, WC7R_indices[a][0], WC7R_indices[a][1],
						WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC7I; a++) {
		WC7I_set(Clq3I, WC7I_indices[a][0], WC7I_indices[a][1],
						WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
		c++;
	}

	//SMEFT class 8_RRRR
	double CuuR[NG * NG * NG * NG];
	double CuuI[NG * NG * NG * NG];
	double CeeR[NG * NG * NG * NG];
	double CeeI[NG * NG * NG * NG];
	double CddR[NG * NG * NG * NG];
	double CddI[NG * NG * NG * NG];
	double CeuR[NG * NG * NG * NG];
	double CeuI[NG * NG * NG * NG];
	double CedR[NG * NG * NG * NG];
	double CedI[NG * NG * NG * NG];
	double Cud1R[NG * NG * NG * NG];
	double Cud1I[NG * NG * NG * NG];
	double Cud8R[NG * NG * NG * NG];
	double Cud8I[NG * NG * NG * NG];


	
	for (a = 0; a < DWC8R; a++) {
		WC8R_set(CeeR, WC8R_indices[a][0], WC8R_indices[a][1],
						WC8R_indices[a][2], WC8R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC8I; a++) {
		WC8I_set(CeeI, WC8I_indices[a][0], WC8I_indices[a][1],
						WC8I_indices[a][2], WC8I_indices[a][3], y[c]);
		c++;
	}
	 
	for (a = 0; a < DWC6R; a++) {
		WC6R_set(CuuR, WC6R_indices[a][0], WC6R_indices[a][1],
						WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC6I; a++) {
		WC6I_set(CuuI, WC6I_indices[a][0], WC6I_indices[a][1],
						WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
		c++;
	}
	
	for (a = 0; a < DWC6R; a++) {
		WC6R_set(CddR, WC6R_indices[a][0], WC6R_indices[a][1],
						WC6R_indices[a][2], WC6R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC6I; a++) {
		WC6I_set(CddI, WC6I_indices[a][0], WC6I_indices[a][1],
						WC6I_indices[a][2], WC6I_indices[a][3], y[c]);
		c++;
	}


	
	for (a = 0; a < DWC7R; a++) {
		WC7R_set(CeuR, WC7R_indices[a][0], WC7R_indices[a][1],
						WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
		c++;

	}
	for (a = 0; a < DWC7I; a++) {
		WC7I_set(CeuI, WC7I_indices[a][0], WC7I_indices[a][1],
						WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
		c++;
	}

	
	for (a = 0; a < DWC7R; a++) {
		WC7R_set(CedR, WC7R_indices[a][0], WC7R_indices[a][1],
						WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC7I; a++) {
		WC7I_set(CedI, WC7I_indices[a][0], WC7I_indices[a][1],
						WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
		c++;
	}

	
	for (a = 0; a < DWC7R; a++) {
		WC7R_set(Cud1R, WC7R_indices[a][0], WC7R_indices[a][1],
						WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC7I; a++) {
		WC7I_set(Cud1I, WC7I_indices[a][0], WC7I_indices[a][1],
						WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
		c++;
	}

	
	for (a = 0; a < DWC7R; a++) {
		WC7R_set(Cud8R, WC7R_indices[a][0], WC7R_indices[a][1],
						WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
		c++;
	}
	for (a = 0; a < DWC7I; a++) {
		WC7I_set(Cud8I, WC7I_indices[a][0], WC7I_indices[a][1],
						WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
		c++;
	}

	//SMEFT class 8_LLRR
	double CleR[NG * NG * NG * NG];
	double CleI[NG * NG * NG * NG];
	double CluR[NG * NG * NG * NG];
	double CluI[NG * NG * NG * NG];
	double CldR[NG * NG * NG * NG];
	double CldI[NG * NG * NG * NG];
	double CqeR[NG * NG * NG * NG];
	double CqeI[NG * NG * NG * NG];
	double Cqu1R[NG * NG * NG * NG];
	double Cqu1I[NG * NG * NG * NG];
	double Cqu8R[NG * NG * NG * NG];
	double Cqu8I[NG * NG * NG * NG];
	double Cqd1R[NG * NG * NG * NG];
	double Cqd1I[NG * NG * NG * NG];
	double Cqd8R[NG * NG * NG * NG];
	double Cqd8I[NG * NG * NG * NG];



	{
		for (a = 0; a < DWC7R; a++) {
			WC7R_set(CleR, WC7R_indices[a][0], WC7R_indices[a][1],
							WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7I; a++) {
			WC7I_set(CleI, WC7I_indices[a][0], WC7I_indices[a][1],
							WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7R; a++) {
			WC7R_set(CluR, WC7R_indices[a][0], WC7R_indices[a][1],
							WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7I; a++) {
			WC7I_set(CluI, WC7I_indices[a][0], WC7I_indices[a][1],
							WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
			c++;
		}

		for (a = 0; a < DWC7R; a++) {
			WC7R_set(CldR, WC7R_indices[a][0], WC7R_indices[a][1],
							WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7I; a++) {
			WC7I_set(CldI, WC7I_indices[a][0], WC7I_indices[a][1],
							WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7R; a++) {
			WC7R_set(CqeR, WC7R_indices[a][0], WC7R_indices[a][1],
							WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7I; a++) {
			WC7I_set(CqeI, WC7I_indices[a][0], WC7I_indices[a][1],
							WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
			c++;
		}

		for (a = 0; a < DWC7R; a++) {
			WC7R_set(Cqu1R, WC7R_indices[a][0], WC7R_indices[a][1],
							WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7I; a++) {
			WC7I_set(Cqu1I, WC7I_indices[a][0], WC7I_indices[a][1],
							WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7R; a++) {
			WC7R_set(Cqu8R, WC7R_indices[a][0], WC7R_indices[a][1],
							WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7I; a++) {
			WC7I_set(Cqu8I, WC7I_indices[a][0], WC7I_indices[a][1],
							WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
			c++;
		}

		for (a = 0; a < DWC7R; a++) {
			WC7R_set(Cqd1R, WC7R_indices[a][0], WC7R_indices[a][1],
							WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7I; a++) {
			WC7I_set(Cqd1I, WC7I_indices[a][0], WC7I_indices[a][1],
							WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7R; a++) {
			WC7R_set(Cqd8R, WC7R_indices[a][0], WC7R_indices[a][1],
							WC7R_indices[a][2], WC7R_indices[a][3], y[c]);
			c++;
		}
		for (a = 0; a < DWC7I; a++) {
			WC7I_set(Cqd8I, WC7I_indices[a][0], WC7I_indices[a][1],
							WC7I_indices[a][2], WC7I_indices[a][3], y[c]);
			c++;
		}
	}


	//SMEFT class 8_LRRL
	double CledqR[NG * NG * NG * NG], CledqI[NG * NG * NG * NG];

	for (i = 0; i < NG; i++) {
		for (j = 0; j < NG; j++) {
			for (k = 0; k < NG; k++) {
				for (l = 0; l < NG; l++) {
					CledqR[NG * NG * NG * i + NG * NG * j + NG * k + l] = y[c];
					CledqI[NG * NG * NG * i + NG * NG * j + NG * k + l] = y[c + NG * NG * NG * NG];
					c++;
				}
			}
		}
	}
	c += NG * NG * NG*NG;

	//SMEFT class 8_LRLR

	double Cquqd1R[NG * NG * NG * NG], Cquqd1I[NG * NG * NG * NG],
		Cquqd8R[NG * NG * NG * NG], Cquqd8I[NG * NG * NG * NG],
		Clequ1R[NG * NG * NG * NG], Clequ1I[NG * NG * NG * NG],
		Clequ3R[NG * NG * NG * NG], Clequ3I[NG * NG * NG * NG];

	for (i = 0; i < NG; i++) {
		for (j = 0; j < NG; j++) {
			for (k = 0; k < NG; k++) {
				for (l = 0; l < NG; l++) {
					a = 0;
					WC5_set(Cquqd1R, i, j, k, l, y[c + a * NG * NG * NG * NG]);
					a++;
					WC5_set(Cquqd1I, i, j, k, l, y[c + a * NG * NG * NG * NG]);
					a++;
					WC5_set(Cquqd8R, i, j, k, l, y[c + a * NG * NG * NG * NG]);
					a++;
					WC5_set(Cquqd8I, i, j, k, l, y[c + a * NG * NG * NG * NG]);
					a++;
					WC5_set(Clequ1R, i, j, k, l, y[c + a * NG * NG * NG * NG]);
					a++;
					WC5_set(Clequ1I, i, j, k, l, y[c + a * NG * NG * NG * NG]);
					a++;
					WC5_set(Clequ3R, i, j, k, l, y[c + a * NG * NG * NG * NG]);
					a++;
					WC5_set(Clequ3I, i, j, k, l, y[c + a * NG * NG * NG * NG]);
					a++;
					c++;
				}
			}
		}
	}

	c += NG * NG * NG * NG * (2 * N8_LRLR - 1);

	//Auxiliary quantities: 

	//Powers and products of gauge couplings
	double g12 = g1*g1;
	double g22 = g2*g2;
	double g32 = g3*g3;
	double g1g3 = g1*g3;
	double g1g2 = g1*g2;
	double g2g3 = g2*g3;

	double gammaH = 0.; //Higgs wavefunction normalization  
	double H = 0.; //scalar yukawa-trace dependent appearing in lambda RGE
	gslpp::matrix<gslpp::complex> gq(NG,0.); //wavefunction ren. const. of q 
	gslpp::matrix<gslpp::complex> gl(NG,0.); //wavefunction ren. const. of l 
	gslpp::matrix<gslpp::complex> gu(NG,0.); //wavefunction ren. const. of u
	gslpp::matrix<gslpp::complex> gd(NG,0.); //wavefunction ren. const. of d
	gslpp::matrix<gslpp::complex> ge(NG,0.); //wavefunction ren. const. of e 
	gslpp::matrix<gslpp::complex> yudyu(NG,0.); //yu^dag yu 
	gslpp::matrix<gslpp::complex> yddyd(NG,0.); //yd^dag yd 
	gslpp::matrix<gslpp::complex> yedye(NG,0.); //ye^dag ye 
	gslpp::matrix<gslpp::complex> ydyud(NG,0.); //yd yu^dag
	gslpp::matrix<gslpp::complex> yuydd(NG,0.); //yu yd^dag

	//3 Yukawa matrices product

	gslpp::matrix<gslpp::complex> yuyudyu(NG,0.); //yu yu^dag yu
	gslpp::matrix<gslpp::complex> yudyuyud(NG,0.); //yu^dag yu yu^dag
	gslpp::matrix<gslpp::complex> ydyddyd(NG,0.); //yd yd^dag yd
	gslpp::matrix<gslpp::complex> yddydydd(NG,0.); //yd^dag yd yd^dag
	gslpp::matrix<gslpp::complex> yeyedye(NG,0.); //ye ye^dag ye
	gslpp::matrix<gslpp::complex> yedyeyed(NG,0.); //ye^dag ye ye^dag

	gslpp::matrix<gslpp::complex> yddydyud(NG,0.); // yd^dag yd yu^dag
	gslpp::matrix<gslpp::complex> yudyuydd(NG,0.); // yu^dag yu yd^dag


	//Real numbers defined in (A.3) https://arxiv.org/abs/1310.4838 
	double eta1 = 0.; 
	double eta2 = 0.;
	double eta3 = 0.;
	double eta4 = 0.;
	double eta5 = 0.;

	//Complex matrices defined in (A.4) https://arxiv.org/abs/1310.4838 
	gslpp::matrix<gslpp::complex> xie(NG,0.);
	gslpp::matrix<gslpp::complex> xiu(NG,0.);
	gslpp::matrix<gslpp::complex> xid(NG,0.);


	//Products of Yukawas 
        gu = yu * yudag;
        yudyu = yudag * yu;
        gd = yd * yddag;
        yddyd = yddag * yd;
        ge = ye * yedag;
        yedye = yedag * ye;
        ydyud = yd * yudag;
        gq = (yudyu + yddyd) * 0.5;
        gl = yedye * 0.5;
        yuyudyu = gu * yu;
        ydyddyd = gd * yd;
        yeyedye = ge * ye;
        yddydyud = yddyd * yudag;
        yudyuydd = yudyu * yddag;

        yudyuyud = yuyudyu.hconjugate();
        yedyeyed = yeyedye.hconjugate();
        yddydydd = ydyddyd.hconjugate();
        yuydd = ydyud.hconjugate();
        
	for (r = 0; r < NG; r++) {
		for (s = 0; s < NG; s++) {
			eta1 += +0.5 * NC * (
				+WC1(CdHR, r, s) * yd(s,r).real() - WC1(CdHI, r, s) * yd(s,r).imag()
				+ WC1(CdHR, r, s) * yddag(r,s).real() + WC1(CdHI, r, s) * yddag(r,s).imag()
				)
				+ 0.5 * NC * (
				+WC1(CuHR, r, s) * yu(s,r).real() - WC1(CuHI, r, s) * yu(s,r).imag()
				+ WC1(CuHR, r, s) * yudag(r,s).real() + WC1(CuHI, r, s) * yudag(r,s).imag()
				)
				+ 0.5 * (
				+WC1(CeHR, r, s) * ye(s,r).real() - WC1(CeHI, r, s) * ye(s,r).imag()
				+ WC1(CeHR, r, s) * yedag(r,s).real() + WC1(CeHI, r, s) * yedag(r,s).imag()
				);

			eta2 += -2. * NC * (
				+WC2R(CHq3R, r, s)*(yudyu(s,r).real() + yddyd(s,r).real())
				- WC2I(CHq3I, r, s)*(yudyu(s,r).imag() + yddyd(s,r).imag())
				)
				+ NC * (
				+WC1(CHudR, r, s) * ydyud(s,r).real() - WC1(CHudI, r, s) * ydyud(s,r).imag()
				+ WC1(CHudR, s, r) * yuydd(s,r).real() + WC1(CHudI, s, r) * yuydd(s,r).imag()
				)
				- 2. * (
				+WC2R(CHl3R, r, s)*(yedye(s,r).real())
				- WC2I(CHl3I, r, s)*(yedye(s,r).imag())
				);

			eta3 += +NC * (
				+WC2R(CHq1R, r, s)*(yddyd(s,r).real() - yudyu(s,r).real())
				- WC2I(CHq1I, r, s)*(yddyd(s,r).imag() - yudyu(s,r).imag())
				)
				+ 3. * NC * (
				+WC2R(CHq3R, r, s)*(yudyu(s,r).real() + yddyd(s,r).real())
				- WC2I(CHq3I, r, s)*(yudyu(s,r).imag() + yddyd(s,r).imag())
				)
				+ NC * (
				+WC2R(CHuR, r, s) * gu(s,r).real() - WC2I(CHuI, r, s) * gu(s,r).imag())
				- NC * (
				+WC2R(CHdR, r, s) * gd(s,r).real() - WC2I(CHdI, r, s) * gd(s,r).imag())
				- NC * (
				+WC1(CHudR, r, s) * ydyud(s,r).real() - WC1(CHudI, r, s) * ydyud(s,r).imag()
				+ WC1(CHudR, s, r) * yuydd(s,r).real() + WC1(CHudI, s, r) * yuydd(s,r).imag()
				)
				+ (
				+(3. * WC2R(CHl3R, r, s) + WC2R(CHl1R, r, s)) * yedye(s,r).real()
				-(3. * WC2I(CHl3I, r, s) + WC2I(CHl1I, r, s)) * yedye(s,r).imag()
				)
				-(WC2R(CHeR, r, s) * ge(s,r).real() - WC2I(CHeI, r, s) * ge(s,r).imag());

			eta4 += 4. * NC * (
				+WC2R(CHq1R, r, s)*(yddyd(s,r).real() - yudyu(s,r).real())
				- WC2I(CHq1I, r, s)*(yddyd(s,r).imag() - yudyu(s,r).imag())
				)
				+ 4. * NC * (
				+WC2R(CHuR, r, s) * gu(s,r).real() - WC2I(CHuI, r, s) * gu(s,r).imag())
				- 4. * NC * (
				+WC2R(CHdR, r, s) * gd(s,r).real() - WC2I(CHdI, r, s) * gd(s,r).imag())
				+ 2. * NC * (
				+WC1(CHudR, r, s) * ydyud(s,r).real() - WC1(CHudI, r, s) * ydyud(s,r).imag()
				+ WC1(CHudR, s, r) * yuydd(s,r).real() + WC1(CHudI, s, r) * yuydd(s,r).imag()
				)
				+ 4. * (
				+WC2R(CHl1R, r, s) * yedye(s,r).real() - WC2I(CHl1R, r, s) * yedye(s,r).imag()
				)
				- 4. *
				(WC2R(CHeR, r, s) * ge(s,r).real() - WC2I(CHeI, r, s) * ge(s,r).imag());


			eta5 +=
				+0.5 * NC * (
				+(+WC1(CdHR, r, s) * yd(s,r).imag() + WC1(CdHI, r, s) * yd(s,r).real())
				-(+WC1(CdHR, r, s) * yddag(r,s).imag() - WC1(CdHI, r, s) * yddag(r,s).real())
				)
				- 0.5 * NC * (
				+(+WC1(CuHR, r, s) * yu(s,r).imag() + WC1(CuHI, r, s) * yu(s,r).real())
				-(+WC1(CuHR, r, s) * yudag(r,s).imag() - WC1(CuHI, r, s) * yudag(r,s).real())
				)
				+ 0.5 * (
				+(+WC1(CeHR, r, s) * ye(s,r).imag() + WC1(CeHI, r, s) * ye(s,r).real())
				-(+WC1(CeHR, r, s) * yedag(r,s).imag() - WC1(CeHI, r, s) * yedag(r,s).real())
				);
		}
	}


	//gammaH
	for (i = 0; i < NG; i++) {
		gammaH += yedye(i,i).real() + NC * (yudyu(i,i).real() + yddyd(i,i).real());
		for (j = 0; j < NG; j++) {
			H += yedye(i,j).real() * yedye(j,i).real() - yedye(i,j).imag() * yedye(j,i).imag()
				+ NC * (yudyu(i,j).real() * yudyu(j,i).real() - yudyu(i,j).imag() * yudyu(j,i).imag()
				+ yddyd(i,j).real() * yddyd(j,i).real() - yddyd(i,j).imag() * yddyd(j,i).imag());
		}
	}

	//xiB defined in (C.1) https://arxiv.org/abs/1312.2014
	double xiB = +FOUR_THIRDS * Yh * (CHBOX + CHD)
		+ EIGHT_THIRDS * (
		+2. * Yl * (WC2R(CHl1R, 0, 0) + WC2R(CHl1R, 1, 1) + WC2R(CHl1R, 2, 2))
		+ 2. * Yq * NC * (WC2R(CHq1R, 0, 0) + WC2R(CHq1R, 1, 1) + WC2R(CHq1R, 2, 2))
		+ Ye * (WC2R(CHeR, 0, 0) + WC2R(CHeR, 1, 1) + WC2R(CHeR, 2, 2))
		+ Yu * NC * (WC2R(CHuR, 0, 0) + WC2R(CHuR, 1, 1) + WC2R(CHuR, 2, 2))
		+ Yd * NC * (WC2R(CHdR, 0, 0) + WC2R(CHdR, 1, 1) + WC2R(CHdR, 2, 2))
		);

	//xiu,xid,xid
	for (p = 0; p < NG; p++) {
		for (t = 0; t < NG; t++) {
			//2 summed indices
                    double xiertmp = 0., xieitmp = 0., xidrtmp = 0., xiditmp = 0., xiurtmp = 0., xiuitmp = 0.;
			for (r = 0; r < NG; r++) {
				for (s = 0; s < NG; s++) {
					xiertmp += 2. * (WC7R(CleR, p, r, s, t) * ye(s,r).real() + WC7I(CleI, p, r, s, t) * ye(s,r).imag())
						- NC * (WC5(CledqR, p, t, s, r) * yd(s,r).real() + WC5(CledqI, p, t, s, r) * yd(s,r).imag())
						+ NC * (WC5(Clequ1R, p, t, s, r) * yu(r,s).real() - WC5(Clequ1I, p, t, s, r) * yu(r,s).imag());
					xieitmp += 2. * (-WC7R(CleR, p, r, s, t) * ye(s,r).imag() + WC7I(CleI, p, r, s, t) * ye(s,r).real())
						- NC * (-WC5(CledqR, p, t, s, r) * yd(s,r).imag() + WC5(CledqI, p, t, s, r) * yd(s,r).real())
						+ NC * (WC5(Clequ1R, p, t, s, r) * yu(r,s).imag() + WC5(Clequ1I, p, t, s, r) * yu(r,s).real());

					xidrtmp += 2. * ((WC7R(Cqd1R, p, r, s, t) + cF3 * WC7R(Cqd8R, p, r, s, t)) * yd(s,r).real()
						+(WC7I(Cqd1I, p, r, s, t) + cF3 * WC7I(Cqd8I, p, r, s, t)) * yd(s,r).imag())
						-(
						(NC * WC5(Cquqd1R, s, r, p, t) + 0.5 * WC5(Cquqd1R, p, r, s, t)
						+ 0.5 * cF3 * WC5(Cquqd8R, p, r, s, t)) * yu(r,s).real()

						-(NC * WC5(Cquqd1I, s, r, p, t) + 0.5 * WC5(Cquqd1I, p, r, s, t)
						+ 0.5 * cF3 * WC5(Cquqd8I, p, r, s, t)) * yu(r,s).imag()
						)
						-(WC5(CledqR, s, r, t, p) * ye(r,s).real() - WC5(CledqI, s, r, t, p) * ye(r,s).imag());

					xiditmp += 2. * (-(WC7R(Cqd1R, p, r, s, t) + cF3 * WC7R(Cqd8R, p, r, s, t)) * yd(s,r).imag()
						+(WC7I(Cqd1I, p, r, s, t) + cF3 * WC7I(Cqd8I, p, r, s, t)) * yd(s,r).real())
						-(
						(NC * WC5(Cquqd1R, s, r, p, t) + 0.5 * WC5(Cquqd1R, p, r, s, t)
						+ 0.5 * cF3 * WC5(Cquqd8R, p, r, s, t)) * yu(r,s).imag()
						+(NC * WC5(Cquqd1I, s, r, p, t) + 0.5 * WC5(Cquqd1I, p, r, s, t)
						+ 0.5 * cF3 * WC5(Cquqd8I, p, r, s, t)) * yu(r,s).real()
						)
						-(-WC5(CledqR, s, r, t, p) * ye(r,s).imag() - WC5(CledqI, s, r, t, p) * ye(r,s).real());

					xiurtmp += 2. * (
						(WC7R(Cqu1R, p, r, s, t) + cF3 * WC7R(Cqu8R, p, r, s, t)) * yu(s,r).real()
						+(WC7I(Cqu1I, p, r, s, t) + cF3 * WC7I(Cqu8I, p, r, s, t)) * yu(s,r).imag()
						)
						-(
						(NC * WC5(Cquqd1R, p, t, s, r) + 0.5 * WC5(Cquqd1R, s, t, p, r)
						+ 0.5 * cF3 * WC5(Cquqd8R, s, t, p, r)) * yd(r,s).real()
						-(NC * WC5(Cquqd1I, p, t, s, r) + 0.5 * WC5(Cquqd1I, s, t, p, r)
						+ 0.5 * cF3 * WC5(Cquqd8I, s, t, p, r)) * yd(r,s).imag()
						)
						+(WC5(Clequ1R, s, r, p, t) * ye(r,s).real() - WC5(Clequ1I, s, r, p, t) * ye(r,s).imag());
					xiuitmp += 2. * (-(WC7R(Cqu1R, p, r, s, t) + cF3 * WC7R(Cqu8R, p, r, s, t)) * yu(s,r).imag()
						+(WC7I(Cqu1I, p, r, s, t) + cF3 * WC7I(Cqu8I, p, r, s, t)) * yu(s,r).real())
						-((NC * WC5(Cquqd1R, p, t, s, r) + 0.5 * WC5(Cquqd1R, s, t, p, r)
						+ 0.5 * cF3 * WC5(Cquqd8R, s, t, p, r)) * yd(r,s).imag()
						+(NC * WC5(Cquqd1I, p, t, s, r) + 0.5 * WC5(Cquqd1I, s, t, p, r)
						+ 0.5 * cF3 * WC5(Cquqd8I, s, t, p, r)) * yd(r,s).real())
						+(WC5(Clequ1R, s, r, p, t) * ye(r,s).imag() + WC5(Clequ1I, s, r, p, t) * ye(r,s).real());
				}
			}
                    xid.assignre(p,t,xidrtmp);
                    xid.assignim(p,t,xiditmp);
                    xie.assignre(p,t,xiertmp);
                    xie.assignim(p,t,xieitmp);
                    xiu.assignre(p,t,xiurtmp);
                    xiu.assignim(p,t,xiuitmp);
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
		f[c] = (-b02 * g22 //SM
			- 4. * mh2 * CHW //SMEFT
			) * g2 * loop_factor;
		c++;
		//g1
		f[c] = (-b01 * g12 //SM
			- 4. * mh2 * CHB //SMEFT
			) * g1 * loop_factor;
		c++;
		//g3
		f[c] = (-b03 * g32 //SM
			- 4. * mh2 * CHG //SMEFT
			) * g3 * loop_factor;
		c++;

		//lambda
		f[c] = lambda * (24. * lambda - 3. * g12 - 9. * g22 + 4. * gammaH)
			+ 0.375 * g12 * g12 + 0.75 * g12 * g22
			+ 1.125 * g22 * g22 - 2. * H //SM
			+ mh2 * (12. * CH + (-32. * lambda + TEN_THIRDS * g22) * CHBOX
			+ (12. * lambda - 1.5 * g22 + 6. * g12 * Yh2) * CHD + 2. * (eta1 + eta2)
			+ 12. * g22 * cF2 * CHW + 12. * g12 * Yh2 * CHB + 6. * g1g2 * Yh * CHWB
			+ FOUR_THIRDS * g22 * (
			(WC2R(CHl3R, 0, 0) + WC2R(CHl3R, 1, 1) + WC2R(CHl3R, 2, 2)) +
			NC * (WC2R(CHq3R, 0, 0) + WC2R(CHq3R, 1, 1) + WC2R(CHq3R, 2, 2)))
			) //SMEFT
			;
		f[c] *= loop_factor;
		c++;
		//mh2
		f[c] = mh2 * (12. * lambda + 2. * gammaH - 1.5 * g12 - 4.5 * g22 //SM  
			+ mh2 * (2. * CHD - 4. * CHBOX) //SMEFT
			) * loop_factor;
		c++;
	}

	//---------RGE YUKAWA --------
	//SM beta functions can be found in
	// https://arxiv.org/abs/hep-ph/0207271
	for (i = 0; i < NG; i++) {
		for (j = 0; j < NG; j++) {
			//Entries without matrix products:      
			//yuR  
			f[c ] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yu(i,j).real()
				+(3. * WC1(CuHR, j, i) - CHBOX * yu(i,j).real()
				+ 0.5 * CHD * yu(i,j).real()) * mh2
				;
			//yuI
			f[c + DF] = (gammaH - (17. / 12.) * g12 - 2.25 * g22 - 8. * g32) * yu(i,j).imag()
				+(-3. * WC1(CuHI, j, i) - CHBOX * yu(i,j).imag()
				+ 0.5 * CHD * yu(i,j).imag()) * mh2
				;
			//ydR     
			f[c + 2 * DF] = (gammaH - (5./12.) * g12 - 2.25 * g22
				- 8. * g32) * yd(i,j).real()
				+(3. * WC1(CdHR, j, i) - CHBOX * yd(i,j).real()
				+ 0.5 * CHD * yd(i,j).real()) * mh2
				;
			//ydI  
			f[c + 3 * DF] = (gammaH - (5./12.) * g12 - 2.25 * g22
				- 8. * g32) * yd(i,j).imag()
				+(-3. * WC1(CdHI, j, i) - CHBOX * yd(i,j).imag()
				+ 0.5 * CHD * yd(i,j).imag()) * mh2
				;
			//yeR
			f[c + 4 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * ye(i,j).real()
				+(3. * WC1(CeHR, j, i) - CHBOX * ye(i,j).real()
				+ 0.5 * CHD * ye(i,j).real()) * mh2
				;
			//yeI
			f[c + 5 * DF] = (gammaH - 2.25 * g22 - 3.75 * g12) * ye(i,j).imag()
				+(-3. * WC1(CeHI, j, i) - CHBOX * ye(i,j).imag()
				+ 0.5 * CHD * ye(i,j).imag()) * mh2
				;
			//Entries with 1 matrix product (1 summed index)
			for (b = 0; b < NG; b++) {
				//yuR
				f[c ] +=
					+1.5 * (yu(i,b).real()*(yudyu(b,j).real() - yddyd(b,j).real())
					- yu(i,b).imag()*(yudyu(b,j).imag() - yddyd(b,j).imag()))
					+(-(yu(i,b).real()*(WC2R(CHq1R, b, j) - 3. * WC2R(CHq3R, b, j))
					- yu(i,b).imag()*(WC2I(CHq1I, b, j) - 3. * WC2I(CHq3I, b, j)))
					+(WC2R(CHuR, i, b) * yu(b,j).real() - WC2I(CHuI, i, b) * yu(b,j).imag())
					- (WC1(CHudR, i, b) * yd(b,j).real() - WC1(CHudI, i, b) * yd(b,j).imag())
					) * mh2
					;
				//yuI
				f[c + DF] += 1.5 * (yu(i,b).imag()*(yudyu(b,j).real() - yddyd(b,j).real())
					+ yu(i,b).real()*(yudyu(b,j).imag() - yddyd(b,j).imag()))
					+(-(yu(i,b).imag()*(WC2R(CHq1R, b, j) - 3. * WC2R(CHq3R, b, j))
					+ yu(i,b).real()*(WC2I(CHq1I, b, j) - 3. * WC2I(CHq3I, b, j)))
					+(WC2R(CHuR, i, b) * yu(b,j).imag() + WC2I(CHuI, i, b) * yu(b,j).real())
					- (WC1(CHudR, i, b) * yd(b,j).imag() + WC1(CHudI, i, b) * yd(b,j).real())
					) * mh2
					;
				//ydR       
				f[c + 2 * DF] += 1.5 * (yd(i,b).real()*(yddyd(b,j).real() - yudyu(b,j).real())
					- yd(i,b).imag()*(yddyd(b,j).imag() - yudyu(b,j).imag()))
					+(yd(i,b).real()*(WC2R(CHq1R, b, j) + 3. * WC2R(CHq3R, b, j))
					- yd(i,b).imag()*(WC2I(CHq1I, b, j) + 3. * WC2I(CHq3I, b, j))
					-(WC2R(CHdR, i, b) * yd(b,j).real() - WC2I(CHdI, i, b) * yd(b,j).imag())
					-(yu(b,j).real() * WC1(CHudR, b, i) + yu(b,j).imag() * WC1(CHudI, b, i))
					) * mh2
					;
				//ydI 
				f[c + 3 * DF] +=1.5 * (
					yd(i,b).imag()*(yddyd(b,j).real() - yudyu(b,j).real())
					+ yd(i,b).real()*(yddyd(b,j).imag() - yudyu(b,j).imag()))
					+(yd(i,b).imag()*(WC2R(CHq1R, b, j) + 3. * WC2R(CHq3R, b, j))
					+ yd(i,b).real()*(WC2I(CHq1I, b, j) + 3. * WC2I(CHq3I, b, j))
					-(WC2R(CHdR, i, b) * yd(b,j).imag() + WC2I(CHdI, i, b) * yd(b,j).real())
					-(yu(b,j).imag() * WC1(CHudR, b, i) - yu(b,j).real() * WC1(CHudI, b, i))
					) * mh2
					;
				//yeR
				f[c + 4 * DF] += 1.5 * (ye(i,b).real()*(yedye(b,j).real())
					- ye(i,b).imag()*(yedye(b,j).imag()))
					+(ye(i,b).real()*(WC2R(CHl1R, b, j) + 3. * WC2R(CHl3R, b, j))
					- ye(i,b).imag()*(WC2I(CHl1I, b, j) + 3. * WC2I(CHl3I, b, j))
					-(WC2R(CHeR, i, b) * ye(b,j).real() - WC2I(CHeI, i, b) * ye(b,j).imag())
					) * mh2
					;
				//yeI 
				f[c + 5 * DF] += 1.5 * (ye(i,b).imag()*(yedye(b,j).real())
					+ ye(i,b).real()*(yedye(b,j).imag()))
					+(ye(i,b).imag()*(WC2R(CHl1R, b, j) + 3. * WC2R(CHl3R, b, j))
					+ ye(i,b).real()*(WC2I(CHl1I, b, j) + 3. * WC2I(CHl3I, b, j))
					-(WC2R(CHeR, i, b) * ye(b,j).imag() + WC2I(CHeI, i, b) * ye(b,j).real())
					) * mh2
					;

				for (a = 0; a < NG; a++) {
					//yuR
					f[c] += (-2. * (
						(WC7R(Cqu1R, a, j, i, b) + cF3 * WC7R(Cqu8R, a, j, i, b)) * yu(b,a).real()
						-(WC7I(Cqu1I, a, j, i, b) + cF3 * WC7I(Cqu8I, a, j, i, b)) * yu(b,a).imag())
						-(WC5(Clequ1R, a, b, j, i) * ye(b,a).real() - WC5(Clequ1I, a, b, j, i) * ye(b,a).imag())
						+ NC * (WC5(Cquqd1R, j, i, a, b) * yd(b,a).real() - WC5(Cquqd1I, j, i, a, b) * yd(b,a).imag())
						+ 0.5 * (
						(WC5(Cquqd1R, a, i, j, b) + cF3 * WC5(Cquqd8R, a, i, j, b)) * yd(b,a).real()
						-(WC5(Cquqd1I, a, i, j, b) + cF3 * WC5(Cquqd8I, a, i, j, b)) * yd(b,a).imag()
						)
						) * mh2;
					//yuI 
					f[c + DF] += (-2. * (
						(WC7R(Cqu1R, a, j, i, b) + cF3 * WC7R(Cqu8R, a, j, i, b)) * yu(b,a).imag()
						+(WC7I(Cqu1I, a, j, i, b) + cF3 * WC7I(Cqu8I, a, j, i, b)) * yu(b,a).real())
						-(-WC5(Clequ1R, a, b, j, i) * ye(b,a).imag() - WC5(Clequ1I, a, b, j, i) * ye(b,a).real())
						+ NC * (-WC5(Cquqd1R, j, i, a, b) * yd(b,a).imag() - WC5(Cquqd1I, j, i, a, b) * yd(b,a).real())
						+ 0.5 * (
						-(WC5(Cquqd1R, a, i, j, b) + cF3 * WC5(Cquqd8R, a, i, j, b)) * yd(b,a).imag()
						-(WC5(Cquqd1I, a, i, j, b) + cF3 * WC5(Cquqd8I, a, i, j, b)) * yd(b,a).real())
						) * mh2;
					//ydR
					f[c + 2 * DF] += (-2. * (
						(WC7R(Cqd1R, a, j, i, b) + cF3 * WC7R(Cqd8R, a, j, i, b)) * yd(b,a).real()
						-(WC7I(Cqd1I, a, j, i, b) + cF3 * WC7I(Cqd8I, a, j, i, b)) * yd(b,a).imag())
						+(WC5(CledqR, a, b, i, j) * ye(b,a).real() - WC5(CledqI, a, b, i, j) * ye(b,a).imag())
						+ NC * (WC5(Cquqd1R, a, b, j, i) * yu(b,a).real() - WC5(Cquqd1I, a, b, j, i) * yu(b,a).imag())
						+ 0.5 * ((WC5(Cquqd1R, j, a, b, i) + cF3 * WC5(Cquqd8R, j, a, b, i)) * yu(a,b).real()
						-(WC5(Cquqd1I, j, a, b, i) + cF3 * WC5(Cquqd8I, j, a, b, i)) * yu(a,b).imag())
						) * mh2;
					//ydI
					f[c + 3 * DF] += (-2. * (
						(WC7R(Cqd1R, a, j, i, b) + cF3 * WC7R(Cqd8R, a, j, i, b)) * yd(b,a).imag()
						+(WC7I(Cqd1I, a, j, i, b) + cF3 * WC7I(Cqd8I, a, j, i, b)) * yd(b,a).real())
						+(WC5(CledqR, a, b, i, j) * ye(b,a).imag() + WC5(CledqI, a, b, i, j) * ye(b,a).real())
						+ NC * (-WC5(Cquqd1R, a, b, j, i) * yu(b,a).imag() - WC5(Cquqd1I, a, b, j, i) * yu(b,a).real())
						+ 0.5 * (
						-(WC5(Cquqd1R, j, a, b, i) + cF3 * WC5(Cquqd8R, j, a, b, i)) * yu(a,b).imag()
						-(WC5(Cquqd1I, j, a, b, i) + cF3 * WC5(Cquqd8I, j, a, b, i)) * yu(a,b).real())
						) * mh2;

					//yeR
					f[c + 4 * DF] += (
						-2. * (WC7R(CleR, a, j, i, b) * ye(b,a).real() - WC7I(CleI, a, j, i, b) * ye(b,a).imag())
						+ NC * (WC5(CledqR, j, i, a, b) * yd(a,b).real() + WC5(CledqI, j, i, a, b) * yd(a,b).imag())
						- NC * (WC5(Clequ1R, j, i, a, b) * yu(b,a).real() - WC5(Clequ1I, j, i, a, b) * yu(b,a).imag())
						) * mh2;
					//yeI
					f[c + 5 * DF] += (
						-2. * (WC7R(CleR, a, j, i, b) * ye(b,a).imag() + WC7I(CleI, a, j, i, b) * ye(b,a).real())
						+ NC * (WC5(CledqR, j, i, a, b) * yd(a,b).imag() - WC5(CledqI, j, i, a, b) * yd(a,b).real())
						- NC * (-WC5(Clequ1R, j, i, a, b) * yu(b,a).imag() - WC5(Clequ1I, j, i, a, b) * yu(b,a).real())
						) * mh2;
				}
			}


			for (a = 0; a < (2 * Nyukawa); a++) {
				f[c + a * DF] *= loop_factor;
			}
			c++;
		}
	}
	c += (2 * Nyukawa - 1) * DF;


	//SMEFT RGE are in 
	//RGE 1 : https://arxiv.org/abs/1308.2627
	//RGE 2 : https://arxiv.org/abs/1310.4838
	//RGE 3 : https://arxiv.org/abs/1312.2014

	//----------------------RGE SMEFT class 1-3--------------------------
	{
		f[c] = (12. * cA3 - 3. * b03) * g32 * CG * loop_factor; //CG
		f[c + 1] = (12. * cA3 - 3. * b03) * g32 * CGT * loop_factor; //CGT
		f[c + 2] = (12. * cA2 - 3. * b02) * g22 * CW * loop_factor; //CW
		f[c + 3] = (12. * cA2 - 3. * b02) * g22 * CWT * loop_factor; //CWT
		c += N1;

		//CH
		f[c] = lambda * (108. * CH + 8. * eta1 + 8. * eta2) //RGEI 
			- lambda * lambda * (160. * CHBOX - 48. * CHD)
			+ 6. * gammaH*CH;
		for (i = 0; i < NG; i++) {
			for (j = 0; j < NG; j++) {
				f[c] += -4. * NC * (
					2. * (WC1(CuHR, i, j) * yuyudyu(j,i).real() - WC1(CuHI, i, j) * yuyudyu(j,i).imag())
					+ 2. * (WC1(CdHR, i, j) * ydyddyd(j,i).real() - WC1(CdHI, i, j) * ydyddyd(j,i).imag()))
					- 4. * (2. * (WC1(CeHR, i, j) * yeyedye(j,i).real() - WC1(CeHI, i, j) * yeyedye(j,i).imag()));
			}
		}

		//RGE 3 
		f[c] += (-13.5 * g22 - 4.5 * g12) * CH +
			lambda * ((40./3.) * g22 * CHBOX + (-6. * g22 + 24. * g12 * Yh2) * CHD)
			- 0.75 * (4. * Yh2 * g12 + g22) * (4. * Yh2 * g12 + g22) * CHD
			+ 12. * lambda * (3. * g22 * CHW + 4. * g12 * Yh2 * CHB + 2 * g1g2 * Yh * CHWB)
			-(12. * g12 * g22 * Yh2 + 9. * g22 * g22) * CHW
			- g12 * Yh2 * (48. * g12 * Yh2 + 12. * g22) * CHB
			- g1g2 * Yh * (24. * g12 * Yh2 + 6. * g22) * CHWB
			+ (16./3.) * lambda * g22 * (
			WC2R(CHl3R, 0, 0) + WC2R(CHl3R, 1, 1) + WC2R(CHl3R, 2, 2)
			+ NC * (WC2R(CHq3R, 0, 0) + WC2R(CHq3R, 1, 1) + WC2R(CHq3R, 2, 2)));

		//CHBOX
		f[c + 1] = 24. * lambda * CHBOX - (4. * g22 + (16./3.) * g12 * Yh2) * CHBOX
			+ (20. / 3.) * g12 * Yh2 * CHD
			+ 2. * g22 * (WC2R(CHl3R, 0, 0) + WC2R(CHl3R, 1, 1) + WC2R(CHl3R, 2, 2)
			+ NC * (WC2R(CHq3R, 0, 0) + WC2R(CHq3R, 1, 1) + WC2R(CHq3R, 2, 2)))
			+ FOUR_THIRDS * g12 * Yh * (NC * (
			Yu * (WC2R(CHuR, 0, 0) + WC2R(CHuR, 1, 1) + WC2R(CHuR, 2, 2))
			+ Yd * (WC2R(CHdR, 0, 0) + WC2R(CHdR, 1, 1) + WC2R(CHdR, 2, 2))
			+ 2. * Yq * (WC2R(CHq1R, 0, 0) + WC2R(CHq1R, 1, 1) + WC2R(CHq1R, 2, 2)))
			+ Ye * (WC2R(CHeR, 0, 0) + WC2R(CHeR, 1, 1) + WC2R(CHeR, 2, 2))
			+ 2. * Yl * (WC2R(CHl1R, 0, 0) + WC2R(CHl1R, 1, 1) + WC2R(CHl1R, 2, 2)))
			- 2. * eta3 + 4. * gammaH*CHBOX;

		//CHD
		f[c + 2] = 12. * lambda * CHD + (80. / 3.) * g12 * Yh2 * CHBOX
			+ (4.5 * g22 - TEN_THIRDS * g12 * Yh2) * CHD
			+ (16. / 3.) * g12 * Yh * (NC * (
			+ Yu * (WC2R(CHuR, 0, 0) + WC2R(CHuR, 1, 1) + WC2R(CHuR, 2, 2))
			+ Yd * (WC2R(CHdR, 0, 0) + WC2R(CHdR, 1, 1) + WC2R(CHdR, 2, 2))
			+ 2. * Yq * (WC2R(CHq1R, 0, 0) + WC2R(CHq1R, 1, 1) + WC2R(CHq1R, 2, 2)))
			+ Ye * (WC2R(CHeR, 0, 0) + WC2R(CHeR, 1, 1) + WC2R(CHeR, 2, 2))
			+ 2. * Yl * (WC2R(CHl1R, 0, 0) + WC2R(CHl1R, 1, 1) + WC2R(CHl1R, 2, 2)))
			- 2. * eta4 + 4. * gammaH*CHD;

		f[c] *= loop_factor;
		f[c + 1] *= loop_factor;
		f[c + 2] *= loop_factor;

		c += N23;
	}
	//----------------------RGE SMEFT class 4--------------------------
	
	{
		//Entries without matrix products

		//CHG and CHGT
		f[c] = (12. * lambda + 2. * gammaH
			- 6. * Yh2 * g12 - 4.5 * g22
			- 2. * b03 * g32) * CHG;
		f[c + 4] = (12. * lambda + 2. * gammaH
			- 6. * Yh2 * g12 - 4.5 * g22
			- 2. * b03 * g32) * CHGT;

		//CHW and CHWT
		f[c + 2] = (12. * lambda + 2. * gammaH
			- 6. * Yh2 * g12 -
			g22 * (2.5 + 2. * b02)) * CHW
			- 15. * g22 * g2 * CW + 2. * g1g2 * Yh*CHWB;
		f[c + 6] = (12. * lambda + 2. * gammaH
			- 6. * Yh2 * g12 -
			g22 * (2.5 + 2. * b02)) * CHWT
			- 15. * g22 * g2 * CWT + 2. * g1g2 * Yh*CHWBT;

		//CHB and CHBT
		f[c + 1] = (12. * lambda + 2. * gammaH
			+ (2. * Yh2 - 2. * b01) * g12
			- 4.5 * g22) * CHB
			+ 6. * g1g2 * Yh * CHWB;
		f[c + 5] = (12. * lambda + 2. * gammaH
			+ (2. * Yh2 - 2. * b01) * g12
			- 4.5 * g22) * CHBT
			+ 6. * g1g2 * Yh * CHWBT;

		//CHWB and CHWBT
		f[c + 3] = (4. * lambda + 2. * gammaH
			- g12 * (2. * Yh2 + b01)
			+ g22 * (4.5 - b02)) * CHWB
			+ 6. * g1 * g22 * Yh * CW
			+ 4. * g1g2 * Yh * (CHB + CHW);

		f[c + 7] = (4. * lambda + 2. * gammaH
			- g12 * (2 * Yh2 + b01)
			+ g22 * (4.5 - b02)) * CHWBT
			+ 6. * g1 * g22 * Yh * CWT
			+ 4. * g1g2 * Yh * (CHBT + CHWT);


		//Entries with traces (2 summed indices)
		for (a = 0; a < NG; a++) {
			for (b = 0; b < NG; b++) {
				//CHG and CHGT
				f[c] += -2. * g3 * (2. * (WC1(CuGR, a, b) * yu(b,a).real() - WC1(CuGI, a, b) * yu(b,a).imag()))
					- 2. * g3 * (2. * (WC1(CdGR, a, b) * yd(b,a).real() - WC1(CdGI, a, b) * yd(b,a).imag()));
				f[c + 4] += -2. * g3 * (2. * (WC1(CuGI, a, b) * yu(b,a).real() + WC1(CuGR, a, b) * yu(b,a).imag()))
					- 2. * g3 * (2. * (WC1(CdGI, a, b) * yd(b,a).real() + WC1(CdGR, a, b) * yd(b,a).imag()));

				//CHW and CHWT
				f[c + 2 ] += -g2 * NC * (2. * (WC1(CuWR, a, b) * yu(b,a).real() - WC1(CuWI, a, b) * yu(b,a).imag()))
					- g2 * NC * (2. * (WC1(CdWR, a, b) * yd(b,a).real() - WC1(CdWI, a, b) * yd(b,a).imag()))
					- g2 * (2. * (WC1(CeWR, a, b) * ye(b,a).real() - WC1(CeWI, a, b) * ye(b,a).imag()));
				f[c + 6 ] += -g2 * NC * (2. * (WC1(CuWI, a, b) * yu(b,a).real() + WC1(CuWR, a, b) * yu(b,a).imag()))
					- g2 * NC * (2. * (WC1(CdWI, a, b) * yd(b,a).real() + WC1(CdWR, a, b) * yd(b,a).imag()))
					- g2 * (2. * (WC1(CeWI, a, b) * ye(b,a).real() + WC1(CeWR, a, b) * ye(b,a).imag()));

				//CHB and CHBT
				f[c + 1 ] += -2. * g1 * NC * (Yu + Yq) * (2. * (WC1(CuBR, a, b) * yu(b,a).real() - WC1(CuBI, a, b) * yu(b,a).imag()))
					- 2. * g1 * (Yl + Ye) * (2. * (WC1(CeBR, a, b) * ye(b,a).real() - WC1(CeBI, a, b) * ye(b,a).imag()))
					- 2. * g1 * NC * (Yd + Yq) * (2. * (WC1(CdBR, a, b) * yd(b,a).real() - WC1(CdBI, a, b) * yd(b,a).imag()));
				f[c + 5 ] +=
					-2. * g1 * NC * (Yu + Yq) *(2. * (WC1(CuBI, a, b) * yu(b,a).real() + WC1(CuBR, a, b) * yu(b,a).imag()))
					- 2. * g1 * (Yl + Ye) *(2. * (WC1(CeBI, a, b) * ye(b,a).real() + WC1(CeBR, a, b) * ye(b,a).imag()))
					- 2. * g1 * NC * (Yd + Yq) *(2. * (WC1(CdBI, a, b) * yd(b,a).real() + WC1(CdBR, a, b) * yd(b,a).imag()));

				//CHWB and CHWBT
				f[c + 3 ] += g2 * NC * (2. * (WC1(CuBR, a, b) * yu(b,a).real() - WC1(CuBI, a, b) * yu(b,a).imag()))
					- g2 * NC * (2. * (WC1(CdBR, a, b) * yd(b,a).real() - WC1(CdBI, a, b) * yd(b,a).imag()))
					- g2 * (2. * (WC1(CeBR, a, b) * ye(b,a).real() - WC1(CeBI, a, b) * ye(b,a).imag()))
					+ 2. * g1 * (Yq + Yu) * NC * (2. * (WC1(CuWR, a, b) * yu(b,a).real() - WC1(CuWI, a, b) * yu(b,a).imag()))
					- 2. * g1 * (Yq + Yd) * NC * (2. * (WC1(CdWR, a, b) * yd(b,a).real() - WC1(CdWI, a, b) * yd(b,a).imag()))
					- 2. * g1 * (Ye + Yl) * (2. * (WC1(CeWR, a, b) * ye(b,a).real() - WC1(CeWI, a, b) * ye(b,a).imag()));

				f[c + 7 ] += g2 * NC * (2. * (WC1(CuBR, a, b) * yu(b,a).imag() + WC1(CuBI, a, b) * yu(b,a).real()))
					- g2 * NC * (2. * (WC1(CdBR, a, b) * yd(b,a).imag() + WC1(CdBI, a, b) * yd(b,a).real()))
					- g2 * (2. * (WC1(CeBR, a, b) * ye(b,a).imag() + WC1(CeBI, a, b) * ye(b,a).real()))
					+ 2. * g1 * (Yq + Yu) * NC * (2. * (WC1(CuWR, a, b) * yu(b,a).imag() + WC1(CuWI, a, b) * yu(b,a).real()))
					- 2. * g1 * (Yq + Yd) * NC * (2. * (WC1(CdWR, a, b) * yd(b,a).imag() + WC1(CdWI, a, b) * yd(b,a).real()))
					- 2. * g1 * (Ye + Yl) * (2. * (WC1(CeWR, a, b) * ye(b,a).imag() + WC1(CeWI, a, b) * ye(b,a).real()));
			}
		}


		for (i = 0; i < N4; i++) {
			f[c + i] *= loop_factor;
		}
		c += N4;
	}

	//----------------------RGE SMEFT class 5--------------------------
	for (r = 0; r < NG; r++) {
		for (s = 0; s < NG; s++) {
			count = 0;
			//CuH : Re and Im
			f[c + count * DF] = lambda * (24. * WC1(CuHR, r, s)
				- 4. * yudag(r,s).real() * CHBOX + 2. * yudag(r,s).real() * CHD) //RGE 1
				+ 2. * ((eta1 + eta2) * yudag(r,s).real() + eta5 * yudag(r,s).imag())
				+ yudyuyud(r,s).real()*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CuHR, r, s)
				//RGE 2
				+ yudag(r,s).real()*(TEN_THIRDS  * g22 * CHBOX
				+ (-1.5 * g22 + 6. * g12 * Yh2) * CHD)
				-(3. * (3. * Yq2 + 3. * Yu2 - 4. * YuYq) * g12
				+ 6.75 * g22 + 6. * cF3 * g32) * WC1(CuHR, r, s)
				+ 3. * (yudag(r,s).real()*(8. * g32 * cF3 * CHG + 3. * g22 * CHW
				+ 4. * (Yh2 + 2. * YuYq) * g12 * CHB - 2. * Yq * g1g2 * CHWB)
				- yudag(r,s).imag()*(8. * g32 * cF3 * CHGT + 3. * g22 * CHWT
				+ 4. * (Yh2 + 2. * YuYq) * g12 * CHBT - 2. * Yq * g1g2 * CHWBT))
				- 6. * (4. * g12 * g1 * Yh2 * Yu + 4. * g12 * g1 * Yh2 * Yq
				- g22 * g1 * Yh) * WC1(CuBR, r, s)
				+ 3. * (4. * g12 * g2 * YhYu + 4. * g12 * g2 * YhYq - 3. * g22 * g2) * WC1(CuWR, r, s)
				; //RGE 3
			count++;
			f[c + count * DF] = lambda * (24. * WC1(CuHI, r, s)
				- 4. * yudag(r,s).imag() * CHBOX + 2. * yudag(r,s).imag() * CHD) //RGE 1
				+ 2. * ((eta1 + eta2) * yudag(r,s).imag() - eta5 * yudag(r,s).real())
				+ yudyuyud(r,s).imag()*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CuHI, r, s)
				//RGE 2
				+ yudag(r,s).imag()*(TEN_THIRDS  * g22 * CHBOX
				+ (-1.5 * g22 + 6. * g12 * Yh2) * CHD)
				-(3. * (3. * Yq2 + 3. * Yu2 - 4. * YuYq) * g12
				+ 6.75 * g22 + 6. * cF3 * g32) * WC1(CuHI, r, s)
				+ 3. * (yudag(r,s).imag()*(8. * g32 * cF3 * CHG + 3. * g22 * CHW
				+ 4. * (Yh2 + 2. * YuYq) * g12 * CHB - 2. * Yq * g1g2 * CHWB)
				+ yudag(r,s).real()*(8. * g32 * cF3 * CHGT + 3. * g22 * CHWT
				+ 4. * (Yh2 + 2. * YuYq) * g12 * CHBT - 2. * Yq * g1g2 * CHWBT))
				- 6. * (4. * g12 * g1 * Yh2 * Yu + 4. * g12 * g1 * Yh2 * Yq
				- g22 * g1 * Yh) * WC1(CuBI, r, s)
				+ 3. * (4. * g12 * g2 * YhYu + 4. * g12 * g2 * YhYq
				- 3. * g22 * g2) * WC1(CuWI, r, s)
				; //RGE 3
			count++;
			//CdH : Re and Im
			f[c + count * DF] = lambda * (24. * WC1(CdHR, r, s)
				- 4. * yddag(r,s).real() * CHBOX + 2. * yddag(r,s).real() * CHD)
				//RGE 1
				+ 2. * ((eta1 + eta2) * yddag(r,s).real() - eta5 * yddag(r,s).imag())
				+ yddydydd(r,s).real()*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CdHR, r, s)
				//RGE 2
				+ yddag(r,s).real()*(TEN_THIRDS  * g22 * CHBOX + (-1.5 * g22 + 6. * Yh2 * g12) * CHD)
				-(3. * (3. * Yq2 + 3. * Yd2 - 4. * YdYq) * g12
				+ 6.75 * g22 + 6. * cF3 * g32) * WC1(CdHR, r, s)
				+ 3. * (
				yddag(r,s).real()*(8. * cF3 * g32 * CHG + 3. * g22 * CHW
				+ 4. * (Yh2 + 2. * YdYq) * g12 * CHB + 2. * Yq * g1g2 * CHWB)
				- yddag(r,s).imag()*(8. * cF3 * g32 * CHGT + 3. * g22 * CHWT
				+ 4. * (Yh2 + 2. * YdYq) * g12 * CHBT + 2. * Yq * g1g2 * CHWBT)
				)
				- 6. * (4. * g12 * g1 * Yh2 * Yd + 4. * g12 * g1 * Yh2 * Yq
				+ g22 * g1 * Yh) * WC1(CdBR, r, s)
				- 3. * (4. * g12 * g2 * YhYd + 4. * g12 * g2 * YhYq
				+ 3. * g22 * g2) * WC1(CdWR, r, s)
				; //RGE 3
			count++;

			f[c + count * DF] = lambda * (24. * WC1(CdHI, r, s)
				- 4. * yddag(r,s).imag() * CHBOX + 2. * yddag(r,s).imag() * CHD)
				//RGE 1
				+ 2. * ((eta1 + eta2) * yddag(r,s).imag() + eta5 * yddag(r,s).real())
				+ yddydydd(r,s).imag()*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CdHI, r, s)
				//RGE 2
				+ yddag(r,s).imag()*(TEN_THIRDS  * g22 * CHBOX + (-1.5 * g22 + 6. * Yh2 * g12) * CHD)
				-(3. * (3. * Yq2 + 3. * Yd2 - 4. * YdYq) * g12
				+ 6.75 * g22 + 6. * cF3 * g32) * WC1(CdHI, r, s)
				+ 3. * (
				yddag(r,s).imag()*(8. * cF3 * g32 * CHG + 3. * g22 * CHW
				+ 4. * (Yh2 + 2. * YdYq) * g12 * CHB + 2. * Yq * g1g2 * CHWB)
				+ yddag(r,s).real()*(8. * cF3 * g32 * CHGT + 3. * g22 * CHWT
				+ 4. * (Yh2 + 2. * YdYq) * g12 * CHBT + 2. * Yq * g1g2 * CHWBT)
				)
				- 6. * (4. * g12 * g1 * Yh2 * Yd + 4. * g12 * g1 * Yh2 * Yq
				+ g22 * g1 * Yh) * WC1(CdBI, r, s)
				- 3. * (4. * g12 * g2 * YhYd + 4. * g12 * g2 * YhYq
				+ 3. * g22 * g2) * WC1(CdWI, r, s)
				; //RGE 3
			count++;
			//CeH : Re and Im
			f[c + count * DF] = lambda * (24. * WC1(CeHR, r, s)
				- 4. * yedag(r,s).real() * CHBOX + 2. * yedag(r,s).real() * CHD) //RGE 1
				+ 2. * ((eta1 + eta2) * yedag(r,s).real() - eta5 * yedag(r,s).imag())
				+ yedyeyed(r,s).real()*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CeHR, r, s)
				//RGE 2
				+ yedag(r,s).real()*(TEN_THIRDS  * g22 * CHBOX + (-1.5 * g22 + 6. * g12 * Yh2) * CHD)
				-(3. * (3. * Yl2 + 3. * Ye2 - 4. * YeYl) * g12
				+ 6.75 * g22) * WC1(CeHR, r, s)
				+ 3. * (yedag(r,s).real()*(3. * g22 * CHW
				+ 4. * (Yh2 + 2. * YeYl) * g12 * CHB
				+ 2. * g1g2 * Yl * CHWB)
				- yedag(r,s).imag()*(3. * g22 * CHWT
				+ 4. * (Yh2 + 2. * YeYl) * g12 * CHBT
				+ 2. * g1g2 * Yl * CHWBT))
				- 6. * (4. * g12 * g1 * Yh2 * Ye + 4. * g12 * g1 * Yh2 * Yl
				+ g22 * g1 * Yh) * WC1(CeBR, r, s)
				- 3. * (4. * g12 * g2 * YhYe + 4. * g12 * g2 * YhYl
				+ 3. * g22 * g2) * WC1(CeWR, r, s);
			//RGE 3
			count++;
			f[c + count * DF] = lambda * (24. * WC1(CeHI, r, s)
				- 4. * yedag(r,s).imag() * CHBOX + 2. * yedag(r,s).imag() * CHD) //RGE 1
				+ 2. * ((eta1 + eta2) * yedag(r,s).imag() + eta5 * yedag(r,s).real())
				+ yedyeyed(r,s).imag()*(CHD - 6. * CHBOX) + 3. * gammaH * WC1(CeHI, r, s)
				//RGE 2
				+ yedag(r,s).imag()*(TEN_THIRDS  * g22 * CHBOX + (-1.5 * g22 + 6. * g12 * Yh2) * CHD)
				-(3. * (3. * Yl2 + 3. * Ye2 - 4. * YeYl) * g12
				+ 6.75 * g22) * WC1(CeHI, r, s)
				+ 3. * (yedag(r,s).imag()*(3. * g22 * CHW
				+ 4. * (Yh2 + 2. * YeYl) * g12 * CHB
				+ 2. * g1g2 * Yl * CHWB)
				+ yedag(r,s).real()*(3. * g22 * CHWT
				+ 4. * (Yh2 + 2. * YeYl) * g12 * CHBT
				+ 2. * g1g2 * Yl * CHWBT))
				- 6. * (4. * g12 * g1 * Yh2 * Ye + 4. * g12 * g1 * Yh2 * Yl
				+ g22 * g1 * Yh) * WC1(CeBI, r, s)
				- 3. * (4. * g12 * g2 * YhYe + 4. * g12 * g2 * YhYl
				+ 3. * g22 * g2) * WC1(CeWI, r, s)
				;
			//RGE 3
			count++;
			//Entries with 1 summed index
			for (t = 0; t < NG; t++) {
				count = 0;
				//CuH : Re and Im
				f[c + count * DF] += lambda * (-4. * (WC2R(CHq1R, r, t) * yudag(t,s).real()
					- WC2I(CHq1I, r, t) * yudag(t,s).imag())
					+ 12. * (WC2R(CHq3R, r, t) * yudag(t,s).real() - WC2I(CHq3I, r, t) * yudag(t,s).imag())
					+ 4. * (yudag(r,t).real() * WC2R(CHuR, t, s) - yudag(r,t).imag() * WC2I(CHuI, t, s))
					- 4. * (yddag(r,t).real() * WC1(CHudR, s, t) + yddag(r,t).imag() * WC1(CHudI, s, t))) //RGE1
					- 2. * (WC2R(CHq1R, r, t) * yudyuyud(t,s).real() - WC2I(CHq1I, r, t) * yudyuyud(t,s).imag())
					+ 6. * (WC2R(CHq3R, r, t) * yddydyud(t,s).real() - WC2I(CHq3I, r, t) * yddydyud(t,s).imag())
					+ 2. * (yudyuyud(r,t).real() * WC2R(CHuR, t, s) - yudyuyud(r,t).imag() * WC2I(CHuI, t, s))
					- 2. * (yddydydd(r,t).real() * WC1(CHudR, s, t) + yddydydd(r,t).imag() * WC1(CHudI, s, t))
					+ 4. * (WC1(CuHR, r, t) * gu(t,s).real() - WC1(CuHI, r, t) * gu(t,s).imag())
					+ 5. * (yudyu(r,t).real() * WC1(CuHR, t, s) - yudyu(r,t).imag() * WC1(CuHI, t, s))
					-(WC1(CdHR, r, t) * ydyud(t,s).real() - WC1(CdHI, r, t) * ydyud(t,s).imag())
					- 2. * (yddyd(r,t).real() * WC1(CuHR, t, s) - yddyd(r,t).imag() * WC1(CuHI, t, s))
					+ (gq(r,t).real() * WC1(CuHR, t, s) - gq(r,t).imag() * WC1(CuHI, t, s))
					+(WC1(CuHR, r, t) * gu(t,s).real() - WC1(CuHI, r, t) * gu(t,s).imag())
					//RGE 2
					- 12. * g2 * (yddyd(r,t).real() * WC1(CuWR, t, s) - yddyd(r,t).imag() * WC1(CuWI, t, s))
					- 6. * g2 * (WC1(CdWR, r, t) * ydyud(t,s).real() - WC1(CdWI, r, t) * ydyud(t,s).imag())
					- 3. * ((4. * g3 * cF3 * WC1(CuGR, r, t) + g2 * WC1(CuWR, r, t)
					+(3. * Yu + Yd) * g1 * WC1(CuBR, r, t)) * gu(t,s).real()
					-(4. * g3 * cF3 * WC1(CuGI, r, t) + g2 * WC1(CuWI, r, t)
					+(3. * Yu + Yd) * WC1(CuBI, r, t)) * gu(t,s).imag())
					- 3. * (yudyu(r,t).real()*(4. * cF3 * g3 * WC1(CuGR, t, s)
					- g2 * WC1(CuWR, t, s) + 2. * (Yq + Yu) * g1 * WC1(CuBR, t, s))
					- yudyu(r,t).imag()*(4. * cF3 * g3 * WC1(CuGI, t, s)
					- g2 * WC1(CuWI, t, s) + 2. * (Yq + Yu) * g1 * WC1(CuBI, t, s)))
					-(3. * g22 - 12. * g12 * YhYq)*(
					yudag(r,t).real() * WC2R(CHuR, t, s) - yudag(r,t).imag() * WC2I(CHuI, t, s))
					+ 3. * g22 * (yddag(r,t).real() * WC1(CHudR, s, t)
					+ yddag(r,t).imag() * WC1(CHudI, s, t))
					+ 12. * g12 * YhYu * (WC2R(CHq1R, r, t) * yudag(t,s).real()
					- WC2I(CHq1I, r, t) * yudag(t,s).imag())
					- 12. * g12 * YhYu * (WC2R(CHq3R, r, t) * yudag(t,s).real()
					- WC2I(CHq3I, r, t) * yudag(t,s).imag())
					+ FOUR_THIRDS * g22 * yudag(r,s).real()*(
					WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
					; //RGE 3



				count++;
				f[c + count * DF] += lambda * (-4. * (WC2R(CHq1R, r, t) * yudag(t,s).imag()
					+ WC2I(CHq1I, r, t) * yudag(t,s).real())
					+ 12. * (WC2R(CHq3R, r, t) * yudag(t,s).imag() + WC2I(CHq3I, r, t) * yudag(t,s).real())
					+ 4. * (yudag(r,t).imag() * WC2R(CHuR, t, s) + yudag(r,t).real() * WC2I(CHuI, t, s))
					- 4. * (yddag(r,t).imag() * WC1(CHudR, s, t) - yddag(r,t).real() * WC1(CHudI, s, t))) //RGE1
					- 2. * (WC2R(CHq1R, r, t) * yudyuyud(t,s).imag() + WC2I(CHq1I, r, t) * yudyuyud(t,s).real())
					+ 6. * (WC2R(CHq3R, r, t) * yddydyud(t,s).imag()+(WC2I(CHq3I, r, t) * yddydyud(t,s).real()))
					+ 2. * (yudyuyud(r,t).imag() * WC2R(CHuR, t, s) + yudyuyud(r,t).real() * WC2I(CHuI, t, s))
					- 2. * (yddydydd(r,t).imag() * WC1(CHudR, s, t) - yddydydd(r,t).real() * WC1(CHudI, s, t))
					+ 4. * (WC1(CuHR, r, t) * gu(t,s).imag() + WC1(CuHI, r, t) * gu(t,s).real())
					+ 5. * (yudyu(r,t).imag() * WC1(CuHR, t, s) + yudyu(r,t).real() * WC1(CuHI, t, s))
					-(WC1(CdHR, r, t) * ydyud(t,s).imag() + WC1(CdHI, r, t) * ydyud(t,s).real())
					- 2. * (yddyd(r,t).imag() * WC1(CuHR, t, s) + yddyd(r,t).real() * WC1(CuHI, t, s))
					+ (gq(r,t).imag() * WC1(CuHR, t, s) + gq(r,t).real() * WC1(CuHI, t, s))
					+(WC1(CuHR, r, t) * gu(t,s).imag() + WC1(CuHI, r, t) * gu(t,s).real())
					//RGE 2
					- 12. * g2 * (yddyd(r,t).imag() * WC1(CuWR, t, s) + yddyd(r,t).real() * WC1(CuWI, t, s))
					- 6. * g2 * (WC1(CdWR, r, t) * ydyud(t,s).imag() + WC1(CdWI, r, t) * ydyud(t,s).real())
					- 3. * ((4. * g3 * cF3 * WC1(CuGR, r, t) + g2 * WC1(CuWR, r, t)
					+(3. * Yu + Yd) * WC1(CuBR, r, t)) * gu(t,s).imag()
					+(4. * g3 * cF3 * WC1(CuGI, r, t) + g2 * WC1(CuWI, r, t)
					+(3. * Yu + Yd) * g1 * WC1(CuBI, r, t)) * gu(t,s).real())
					- 3. * (yudyu(r,t).imag()*(4. * cF3 * g3 * WC1(CuGR, t, s)
					- g2 * WC1(CuWR, t, s) + 2. * (Yq + Yu) * g1 * WC1(CuBR, t, s))
					+ yudyu(r,t).real()*(4. * cF3 * g3 * WC1(CuGI, t, s)
					- g2 * WC1(CuWI, t, s) + 2. * (Yq + Yu) * g1 * WC1(CuBI, t, s)))
					-(3. * g22 - 12. * g12 * YhYq)*(
					yudag(r,t).imag() * WC2R(CHuR, t, s) + yudag(r,t).real() * WC2I(CHuI, t, s))
					+ 3. * g22 * (yddag(r,t).imag() * WC1(CHudR, s, t)
					- yddag(r,t).real() * WC1(CHudI, s, t))
					+ 12. * g12 * YhYu * (WC2R(CHq1R, r, t) * yudag(t,s).imag()
					+ WC2I(CHq1I, r, t) * yudag(t,s).real())
					- 12. * g12 * YhYu * (WC2R(CHq3R, r, t) * yudag(t,s).imag()
					+ WC2I(CHq3I, r, t) * yudag(t,s).real())
					+ FOUR_THIRDS * g22 * yudag(r,s).imag()*(
					WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
					; //RGE 3
				count++;

				//CdH : Re and Im
				f[c + count * DF] += lambda * (
					+4. * (WC2R(CHq1R, r, t) * yddag(t,s).real() - WC2I(CHq1I, r, t) * yddag(t,s).imag())
					+12.* (WC2R(CHq3R, r, t) * yddag(t,s).real() - WC2I(CHq3I, r, t) * yddag(t,s).imag())
					- 4. * (yddag(r,t).real() * WC2R(CHdR, t, s) - yddag(r,t).imag() * WC2I(CHdI, t, s))
					- 4. * (yudag(r,t).real() * WC1(CHudR, t, s) - yudag(r,t).imag() * WC1(CHudI, t, s)))
					//RGE 1
					+ 2. * (WC2R(CHq1R, r, t) * yddydydd(t,s).real() - WC2I(CHq1I, r, t) * yddydydd(t,s).imag())
					+ 6. * (WC2R(CHq3R, r, t) * yudyuydd(t,s).real() - WC2I(CHq3I, r, t) * yudyuydd(t,s).imag())
					- 2. * (yddydydd(r,t).real() * WC2R(CHdR, t, s) - yddydydd(r,t).imag() * WC2I(CHdI, t, s))
					- 2. * (yudyuyud(r,t).real() * WC1(CHudR, t, s) - yudyuyud(r,t).imag() * WC1(CHudI, t, s))
					+ 4. * (WC1(CdHR, r, t) * gd(t,s).real() - WC1(CdHI, r, t) * gd(t,s).imag())
					+ 5. * (yddyd(r,t).real() * WC1(CdHR, t, s) - yddyd(r,t).imag() * WC1(CdHI, t, s))
					-(WC1(CuHR, r, t) * yuydd(t,s).real() - WC1(CuHI, r, t) * yuydd(t,s).imag())
					- 2. * (yudyu(r,t).real() * WC1(CdHR, t, s) - yudyu(r,t).imag() * WC1(CdHI, t, s))
					+(gq(r,t).real() * WC1(CdHR, t, s) - gq(r,t).imag() * WC1(CdHI, t, s))
					+(WC1(CdHR, r, t) * gd(t,s).real() - WC1(CdHI, r, t) * gd(t,s).imag())
					//RGE 2
					- 12. * g2 * (yudyu(r,t).real() * WC1(CdWR, t, s) - yudyu(r,t).imag() * WC1(CdWI, t, s))
					- 6. * g2 * (WC1(CuWR, r, t) * yuydd(t,s).real() - WC1(CuWI, r, t) * yuydd(t,s).imag())
					- 3. * (
					(4. * cF3 * g3 * WC1(CdGR, r, t)
					+ g2 * WC1(CdWR, r, t)+(3. * Yd + Yu) * g1 * WC1(CdBR, r, t)) * gd(t,s).real()
					-(4. * cF3 * g3 * WC1(CdGI, r, t)
					+ g2 * WC1(CdWI, r, t)+(3. * Yd + Yu) * g1 * WC1(CdBI, r, t)) * gd(t,s).imag()
					)
					- 3. * (
					yddyd(r,t).real()*(4. * cF3 * g3 * WC1(CdGR, t, s)
					- g2 * WC1(CdWR, t, s) + 2. * (Yq + Yd) * g1 * WC1(CdBR, t, s))
					- yddyd(r,t).imag()*(4. * cF3 * g3 * WC1(CdGI, t, s)
					- g2 * WC1(CdWI, t, s) + 2. * (Yq + Yd) * g1 * WC1(CdBI, t, s))
					)
					+(3. * g22 + 12. * g12 * YhYq)*
					(yddag(r,t).real() * WC2R(CHdR, t, s) - yddag(r,t).imag() * WC2I(CHdI, t, s))
					+ 3. * g22 * (yudag(r,t).real() * WC1(CHudR, t, s)
					- yudag(r,t).imag() * WC1(CHudI, t, s))
					+ 12. * g12 * YhYd * (WC2R(CHq1R, r, t) * yddag(t,s).real()
					- WC2I(CHq1I, r, t) * yddag(t,s).imag())
					+ 12. * g12 * YhYd * (WC2R(CHq3R, r, t) * yddag(t,s).real()
					- WC2I(CHq3I, r, t) * yddag(t,s).imag())
					+ FOUR_THIRDS * g22 * yddag(r,s).real()*
					(WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
					; //RGE 3


				count++;
				f[c + count * DF] += lambda * (
					+4. * (WC2R(CHq1R, r, t) * yddag(t,s).imag() + WC2I(CHq1I, r, t) * yddag(t,s).real())
					+12.* (WC2R(CHq3R, r, t) * yddag(t,s).imag() + WC2I(CHq3I, r, t) * yddag(t,s).real())
					- 4. * (yddag(r,t).imag() * WC2R(CHdR, t, s) + yddag(r,t).real() * WC2I(CHdI, t, s))
					- 4. * (yudag(r,t).imag() * WC1(CHudR, t, s) + yudag(r,t).real() * WC1(CHudI, t, s)))
					//RGE 1
					+ 2. * (WC2R(CHq1R, r, t) * yddydydd(t,s).imag() + WC2I(CHq1I, r, t) * yddydydd(t,s).real())
					+ 6. * (WC2R(CHq3R, r, t) * yudyuydd(t,s).imag() + WC2I(CHq3I, r, t) * yudyuydd(t,s).real())
					- 2. * (yddydydd(r,t).imag() * WC2R(CHdR, t, s) + yddydydd(r,t).real() * WC2I(CHdI, t, s))
					- 2. * (yudyuyud(r,t).imag() * WC1(CHudR, t, s) + yudyuyud(r,t).real() * WC1(CHudI, t, s))
					+ 4. * (WC1(CdHR, r, t) * gd(t,s).imag() + WC1(CdHI, r, t) * gd(t,s).real())
					+ 5. * (yddyd(r,t).imag() * WC1(CdHR, t, s) + yddyd(r,t).real() * WC1(CdHI, t, s))
					-(WC1(CuHR, r, t) * yuydd(t,s).imag() + WC1(CuHI, r, t) * yuydd(t,s).real())
					- 2. * (yudyu(r,t).imag() * WC1(CdHR, t, s) + yudyu(r,t).real() * WC1(CdHI, t, s))
					+(gq(r,t).imag() * WC1(CdHR, t, s) + gq(r,t).real() * WC1(CdHI, t, s))
					+(WC1(CdHR, r, t) * gd(t,s).imag() + WC1(CdHI, r, t) * gd(t,s).real())
					//RGE 2
					- 12. * g2 * (yudyu(r,t).imag() * WC1(CdWR, t, s) + yudyu(r,t).real() * WC1(CdWI, t, s))
					- 6. * g2 * (WC1(CuWR, r, t) * yuydd(t,s).imag() + WC1(CuWI, r, t) * yuydd(t,s).real())
					- 3. * (
					(4. * cF3 * g3 * WC1(CdGR, r, t)
					+ g2 * WC1(CdWR, r, t)+(3. * Yd + Yu) * g1 * WC1(CdBR, r, t)) * gd(t,s).imag()
					+(4. * cF3 * g3 * WC1(CdGI, r, t)
					+ g2 * WC1(CdWI, r, t)+(3. * Yd + Yu) * g1 * WC1(CdBI, r, t)) * gd(t,s).real()
					)
					- 3. * (
					yddyd(r,t).imag()*(4. * cF3 * g3 * WC1(CdGR, t, s)
					- g2 * WC1(CdWR, t, s) + 2. * (Yq + Yd) * g1 * WC1(CdBR, t, s))
					+ yddyd(r,t).real()*(4. * cF3 * g3 * WC1(CdGI, t, s)
					- g2 * WC1(CdWI, t, s) + 2. * (Yq + Yd) * g1 * WC1(CdBI, t, s))
					)
					+(3. * g22 + 12. * g12 * YhYq)*
					(yddag(r,t).imag() * WC2R(CHdR, t, s) + yddag(r,t).real() * WC2I(CHdI, t, s))
					+ 3. * g22 * (yudag(r,t).imag() * WC1(CHudR, t, s)
					+ yudag(r,t).real() * WC1(CHudI, t, s))
					+ 12. * g12 * YhYd * (WC2R(CHq1R, r, t) * yddag(t,s).imag()
					+ WC2I(CHq1I, r, t) * yddag(t,s).real())
					+ 12. * g12 * YhYd * (WC2R(CHq3R, r, t) * yddag(t,s).imag()
					+ WC2I(CHq3I, r, t) * yddag(t,s).real())
					+ FOUR_THIRDS * g22 * yddag(r,s).imag()*
					(WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
					; //RGE 3

				count++;
				//CeH : Re and Im 
				f[c + count * DF] += lambda * (4. * (WC2R(CHl1R, r, t) * yedag(t,s).real() - WC2I(CHl1I, r, t) * yedag(t,s).imag())
					+ 12. * (WC2R(CHl3R, r, t) * yedag(t,s).real() - WC2I(CHl3I, r, t) * yedag(t,s).imag())
					- 4. * (yedag(r,t).real() * WC2R(CHeR, t, s) - yedag(r,t).imag() * WC2I(CHeI, t, s))) //RGE1
					+ 2. * (WC2R(CHl1R, r, t) * yedyeyed(t,s).real() - WC2I(CHl1I, r, t) * yedyeyed(t,s).imag())
					- 2. * (yedyeyed(r,t).real() * WC2R(CHeR, t, s) - yedyeyed(r,t).imag() * WC2I(CHeI, t, s))
					+ 5. * (yedye(r,t).real() * WC1(CeHR, t, s) - yedye(r,t).imag() * WC1(CeHI, t, s))
					+ 4. * (WC1(CeHR, r, t) * ge(t,s).real() - WC1(CeHI, r, t) * ge(t,s).imag())
					+(gl(r,t).real() * WC1(CeHR, t, s) - gl(r,t).imag() * WC1(CeHI, t, s))
					+(WC1(CeHR, r, t) * ge(t,s).real() - WC1(CeHI, r, t) * ge(t,s).imag())
					//RGE 2
					- 3. * ((3. * g1 * Ye * WC1(CeBR, r, t) + g2 * WC1(CeWR, r, t)) * ge(t,s).real()
					-(3. * g1 * Ye * WC1(CeBI, r, t) + g2 * WC1(CeWI, r, t)) * ge(t,s).imag())
					- 3. * (yedye(r,t).real()*(2. * g1 * (Yl + Ye) * WC1(CeBR, t, s) - g2 * WC1(CeWR, t, s))
					- yedye(r,t).imag()*(2. * g1 * (Yl + Ye) * WC1(CeBI, t, s) - g2 * WC1(CeWI, t, s)))
					+(3. * g22 + 12. * g12 * YhYl)*
					(yedag(r,t).real() * WC2R(CHeR, t, s) - yedag(r,t).imag() * WC2I(CHeI, t, s))
					+ 12. * g12 * YhYe * (WC2R(CHl1R, r, t) * yedag(t,s).real()
					- WC2I(CHl1I, r, t) * yedag(t,s).imag())
					+ 12. * g12 * YhYe * (WC2R(CHl3R, r, t) * yedag(t,s).real()
					- WC2I(CHl3I, r, t) * yedag(t,s).imag())
					+ FOUR_THIRDS * g22 * yedag(r,s).real()*
					(WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
					; //RGE 3

				count++;
				f[c + count * DF] += lambda * (4. * (WC2R(CHl1R, r, t) * yedag(t,s).imag() + WC2I(CHl1I, r, t) * yedag(t,s).real())
					+ 12. * (WC2R(CHl3R, r, t) * yedag(t,s).imag() + WC2I(CHl3I, r, t) * yedag(t,s).real())
					- 4. * (yedag(r,t).imag() * WC2R(CHeR, t, s) + yedag(r,t).real() * WC2I(CHeI, t, s))) //RGE1
					+ 2. * (WC2R(CHl1R, r, t) * yedyeyed(t,s).imag() + WC2I(CHl1I, r, t) * yedyeyed(t,s).real())
					- 2. * (yedyeyed(r,t).imag() * WC2R(CHeR, t, s) + yedyeyed(r,t).real() * WC2I(CHeI, t, s))
					+ 5. * (yedye(r,t).imag() * WC1(CeHR, t, s) + yedye(r,t).real() * WC1(CeHI, t, s))
					+ 4. * (WC1(CeHR, r, t) * ge(t,s).imag() + WC1(CeHI, r, t) * ge(t,s).real())
					+(gl(r,t).imag() * WC1(CeHR, t, s) + gl(r,t).real() * WC1(CeHI, t, s))
					+(WC1(CeHR, r, t) * ge(t,s).imag() + WC1(CeHI, r, t) * ge(t,s).real())
					//RGE 2
					- 3. * ((3. * g1 * Ye * WC1(CeBR, r, t) + g2 * WC1(CeWR, r, t)) * ge(t,s).imag()
					+(3. * g1 * Ye * WC1(CeBI, r, t) + g2 * WC1(CeWI, r, t)) * ge(t,s).real())
					- 3. * (yedye(r,t).imag()*(2. * g1 * (Yl + Ye) * WC1(CeBR, t, s) - g2 * WC1(CeWR, t, s))
					+ yedye(r,t).real()*(2. * g1 * (Yl + Ye) * WC1(CeBI, t, s) - g2 * WC1(CeWI, t, s)))
					+(3. * g22 + 12. * g12 * YhYl)*
					(yedag(r,t).imag() * WC2R(CHeR, t, s) + yedag(r,t).real() * WC2I(CHeI, t, s))
					+ 12. * g12 * YhYe * (WC2R(CHl1R, r, t) * yedag(t,s).imag()
					+ WC2I(CHl1I, r, t) * yedag(t,s).real())
					+ 12. * g12 * YhYe * (WC2R(CHl3R, r, t) * yedag(t,s).imag()
					+ WC2I(CHl3I, r, t) * yedag(t,s).real())
					+ FOUR_THIRDS * g22 * yedag(r,s).imag()*
					(WC2R(CHl3R, t, t) + NC * WC2R(CHq3R, t, t))
					; //RGE 3



				count++;

				//entries with 2 summed indices 
				for (p = 0; p < NG; p++) {
					count = 0;
					//CuH : Re and Im 
					f[c + count * DF] += lambda * (-8. * (WC7R(Cqu1R, r, p, t, s) * yudag(p,t).real() - WC7I(Cqu1I, r, p, t, s) * yudag(p,t).imag())
						- 8. * cF3 * (WC7R(Cqu8R, r, p, t, s) * yudag(p,t).real() - WC7I(Cqu8I, r, p, t, s) * yudag(p,t).imag())
						- 4. * (WC5(Clequ1R, p, t, r, s) * ye(t,p).real() - WC5(Clequ1I, p, t, r, s) * ye(t,p).imag())
						+ 4. * NC * (WC5(Cquqd1R, r, s, p, t) * yd(t,p).real() - WC5(Cquqd1I, r, s, p, t) * yd(t,p).imag())
						+ 2. * (WC5(Cquqd1R, p, s, r, t) * yd(t,p).real() - WC5(Cquqd1I, p, s, r, t) * yd(t,p).imag())
						+ 2. * cF3 * (WC5(Cquqd8R, p, s, r, t) * yd(t,p).real() - WC5(Cquqd8I, p, s, r, t) * yd(t,p).imag())
						) //RGE 1
						+ 8. * ((WC7R(Cqu1R, r, p, t, s) + cF3 * WC7R(Cqu8R, r, p, t, s)) * yudyuyud(p,t).real()
						-(WC7I(Cqu1I, r, p, t, s) + cF3 * WC7I(Cqu8I, r, p, t, s)) * yudyuyud(p,t).imag())
						- 2. * ((2. * NC * WC5(Cquqd1R, r, s, t, p) + WC5(Cquqd1R, t, s, r, p)
						+ cF3 * WC5(Cquqd8R, t, s, r, p)) * ydyddyd(p,t).real()
						-(2. * NC * WC5(Cquqd1I, r, s, t, p) + WC5(Cquqd1I, t, s, r, p)
						+ cF3 * WC5(Cquqd8I, t, s, r, p)) * ydyddyd(p,t).imag())
						+ 4. * (WC5(Clequ1R, t, p, r, s) * yeyedye(p,t).real() - WC5(Clequ1I, t, p, r, s) * yeyedye(p,t).imag())
						- 2. * ((yddag(r,t).real() * WC1(CdHR, p, t) + yddag(r,t).imag() * WC1(CdHI, p, t)) * yudag(p,s).real()
						-(yddag(r,t).imag() * WC1(CdHR, p, t) - yddag(r,t).real() * WC1(CdHI, p, t)) * yudag(p,s).imag());
					; //RGE 2
					count++;
					f[c + count * DF] += lambda * (-8. * (WC7R(Cqu1R, r, p, t, s) * yudag(p,t).imag() + WC7I(Cqu1I, r, p, t, s) * yudag(p,t).real())
						- 8. * cF3 * (WC7R(Cqu8R, r, p, t, s) * yudag(p,t).imag() + WC7I(Cqu8I, r, p, t, s) * yudag(p,t).real())
						- 4. * (WC5(Clequ1R, p, t, r, s) * ye(t,p).imag() + WC5(Clequ1I, p, t, r, s) * ye(t,p).real())
						+ 4. * NC * (WC5(Cquqd1R, r, s, p, t) * yd(t,p).imag() + WC5(Cquqd1I, r, s, p, t) * yd(t,p).real())
						+ 2. * (WC5(Cquqd1R, p, s, r, t) * yd(t,p).imag() + WC5(Cquqd1I, p, s, r, t) * yd(t,p).real())
						+ 2. * cF3 * (WC5(Cquqd8R, p, s, r, t) * yd(t,p).imag() + WC5(Cquqd8I, p, s, r, t) * yd(t,p).real())
						) //RGE 1
						+ 8. * ((WC7R(Cqu1R, r, p, t, s) + cF3 * WC7R(Cqu8R, r, p, t, s)) * yudyuyud(p,t).imag()
						+(WC7I(Cqu1I, r, p, t, s) + cF3 * WC7I(Cqu8I, r, p, t, s)) * yudyuyud(p,t).real())
						- 2. * ((2. * NC * WC5(Cquqd1R, r, s, t, p) + WC5(Cquqd1R, t, s, r, p)
						+ cF3 * WC5(Cquqd8R, t, s, r, p)) * ydyddyd(p,t).imag()
						+(2. * NC * WC5(Cquqd1I, r, s, t, p) + WC5(Cquqd1I, t, s, r, p)
						+ cF3 * WC5(Cquqd8I, t, s, r, p)) * ydyddyd(p,t).real())
						+ 4. * (WC5(Clequ1R, t, p, r, s) * yeyedye(p,t).imag() + WC5(Clequ1I, t, p, r, s) * yeyedye(p,t).real())
						- 2. * ((yddag(r,t).real() * WC1(CdHR, p, t) + yddag(r,t).imag() * WC1(CdHI, p, t)) * yudag(p,s).imag()
						+(yddag(r,t).imag() * WC1(CdHR, p, t) - yddag(r,t).real() * WC1(CdHI, p, t)) * yudag(p,s).real());
					; //RGE 2
					count++;

					//CdH : Re and Im
					f[c + count * DF] += lambda * (
						-8. * (WC7R(Cqd1R, r, p, t, s) * yddag(p,t).real() - WC7I(Cqd1I, r, p, t, s) * yddag(p,t).imag())
						- 8. * cF3 * (WC7R(Cqd8R, r, p, t, s) * yddag(p,t).real() - WC7I(Cqd8I, r, p, t, s) * yddag(p,t).imag())
						+ 4. * (WC5(CledqR, p, t, s, r) * yedag(p,t).real() + WC5(CledqI, p, t, s, r) * yedag(p,t).imag())
						+ 4. * NC * (WC5(Cquqd1R, p, t, r, s) * yu(t,p).real() - WC5(Cquqd1I, p, t, r, s) * yu(t,p).imag())
						+ 2. * (WC5(Cquqd1R, r, t, p, s) * yu(t,p).real() - WC5(Cquqd1I, r, t, p, s) * yu(t,p).imag())
						+ 2. * cF3 * (WC5(Cquqd8R, r, t, p, s) * yu(t,p).real() - WC5(Cquqd8I, r, t, p, s) * yu(t,p).imag())
						) //RGE 1
						+ 8. * ((WC7R(Cqd1R, r, p, t, s) + cF3 * WC7R(Cqd8R, r, p, t, s)) * yddydydd(p,t).real()
						-(WC7I(Cqd1I, r, p, t, s) + cF3 * WC7I(Cqd8I, r, p, t, s)) * yddydydd(p,t).imag())
						- 4. * (WC5(CledqR, p, t, s, r) * yedyeyed(p,t).real()
						+ WC5(CledqI, p, t, s, r) * yedyeyed(p,t).imag())
						- 2. * (
						(2. * NC * WC5(Cquqd1R, t, p, r, s) + WC5(Cquqd1R, r, p, t, s)
						+ cF3 * WC5(Cquqd8R, r, p, t, s)) * yuyudyu(p,t).real()
						-(2. * NC * WC5(Cquqd1I, t, p, r, s) + WC5(Cquqd1I, r, p, t, s)
						+ cF3 * WC5(Cquqd8I, r, p, t, s)) * yuyudyu(p,t).imag()
						)
						- 2. * ((yudag(r,t).real() * WC1(CuHR, p, t) + yudag(r,t).imag() * WC1(CuHI, p, t)) * yddag(p,s).real()
						-(yudag(r,t).imag() * WC1(CuHR, p, t) - yudag(r,t).real() * WC1(CuHI, p, t)) * yddag(p,s).imag())
						; //RGE 2
					count++;

					f[c + count * DF] += lambda * (
						-8. * (WC7R(Cqd1R, r, p, t, s) * yddag(p,t).imag() + WC7I(Cqd1I, r, p, t, s) * yddag(p,t).real())
						- 8. * cF3 * (WC7R(Cqd8R, r, p, t, s) * yddag(p,t).imag() + WC7I(Cqd8I, r, p, t, s) * yddag(p,t).real())
						+ 4. * (WC5(CledqR, p, t, s, r) * yedag(p,t).imag() - WC5(CledqI, p, t, s, r) * yedag(p,t).real())
						+ 4. * NC * (WC5(Cquqd1R, p, t, r, s) * yu(t,p).imag() + WC5(Cquqd1I, p, t, r, s) * yu(t,p).real())
						+ 2. * (WC5(Cquqd1R, r, t, p, s) * yu(t,p).imag() + WC5(Cquqd1I, r, t, p, s) * yu(t,p).real())
						+ 2. * cF3 * (WC5(Cquqd8R, r, t, p, s) * yu(t,p).imag() + WC5(Cquqd8I, r, t, p, s) * yu(t,p).real())
						) //RGE 1
						+ 8. * ((WC7R(Cqd1R, r, p, t, s) + cF3 * WC7R(Cqd8R, r, p, t, s)) * yddydydd(p,t).imag()
						+(WC7I(Cqd1I, r, p, t, s) + cF3 * WC7I(Cqd8I, r, p, t, s)) * yddydydd(p,t).real())
						- 4. * (WC5(CledqR, p, t, s, r) * yedyeyed(p,t).imag()
						- WC5(CledqI, p, t, s, r) * yedyeyed(p,t).real())
						- 2. * (
						(2. * NC * WC5(Cquqd1R, t, p, r, s) + WC5(Cquqd1R, r, p, t, s)
						+ cF3 * WC5(Cquqd8R, r, p, t, s)) * yuyudyu(p,t).imag()
						+(2. * NC * WC5(Cquqd1I, t, p, r, s) + WC5(Cquqd1I, r, p, t, s)
						+ cF3 * WC5(Cquqd8I, r, p, t, s)) * yuyudyu(p,t).real()
						)
						- 2. * ((yudag(r,t).real() * WC1(CuHR, p, t) + yudag(r,t).imag() * WC1(CuHI, p, t)) * yddag(p,s).imag()
						+(yudag(r,t).imag() * WC1(CuHR, p, t) - yudag(r,t).real() * WC1(CuHI, p, t)) * yddag(p,s).real())
						; //RGE 2
					count++;




					//CeH : Re and Im
					f[c + count * DF] += lambda * (-8. * (WC7R(CleR, r, p, t, s) * yedag(p,t).real() -
						WC7I(CleI, r, p, t, s) * yedag(p,t).imag())
						+ 4. * NC * (WC5(CledqR, r, s, p, t) * yddag(t,p).real() - WC5(CledqI, r, s, p, t) * yddag(t,p).imag())
						- 4. * NC * (WC5(Clequ1R, r, s, p, t) * yu(t,p).real() - WC5(Clequ1I, r, s, p, t) * yu(t,p).imag())
						) //RGE 1
						+ 8. * (WC7R(CleR, r, p, t, s) * yedyeyed(p,t).real() - WC7I(CleI, r, p, t, s) * yedyeyed(p,t).imag())
						- 4. * NC * (WC5(CledqR, r, s, p, t) * yddydydd(t,p).real() - WC5(CledqI, r, s, p, t) * yddydydd(t,p).imag())
						+ 4. * NC * (WC5(Clequ1R, r, s, t, p) * yuyudyu(p,t).real() - WC5(Clequ1I, r, s, t, p) * yuyudyu(p,t).imag())
						; //RGE 2

					count++;
					f[c + count * DF] += lambda * (-8. * (WC7R(CleR, r, p, t, s) * yedag(p,t).imag() +
						WC7I(CleI, r, p, t, s) * yedag(p,t).real())
						+ 4. * NC * (WC5(CledqR, r, s, p, t) * yddag(t,p).imag() + WC5(CledqI, r, s, p, t) * yddag(t,p).real())
						- 4. * NC * (WC5(Clequ1R, r, s, p, t) * yu(t,p).imag() + WC5(Clequ1I, r, s, p, t) * yu(t,p).real())
						) //RGE 1
						+ 8. * (WC7R(CleR, r, p, t, s) * yedyeyed(p,t).imag() + WC7I(CleI, r, p, t, s) * yedyeyed(p,t).real())
						- 4. * NC * (WC5(CledqR, r, s, p, t) * yddydydd(t,p).imag() + WC5(CledqI, r, s, p, t) * yddydydd(t,p).real())
						+ 4. * NC * (WC5(Clequ1R, r, s, t, p) * yuyudyu(p,t).imag() + WC5(Clequ1I, r, s, t, p) * yuyudyu(p,t).real())
						; //RGE 2
					count++;

				}


			}


			for (l = 0; l < 2 * N5; l++) {
				f[c + l * DF] *= loop_factor;
			}
			c++;
		}
	}
	c += (N5 * 2 - 1) * DF;



	//----------------------RGE SMEFT class 6--------------------------
	for (i = 0; i < NG; i++) {
		for (j = 0; j < NG; j++) {
			//Entries without  matrix products 
			l = 0;
			//CeW : Re and Im;
			f[c + l * DF] = ((3. * cF2 - b02) * g22 +
				(-3. * Ye2 + 8. * YeYl - 3. * Yl2) * g12) * WC1(CeWR, i, j)
				+ g1g2 * (3. * Yl - Ye) * WC1(CeBR, i, j)
				- g2 * (ye(j,i).real() * CHW + ye(j,i).imag() * CHWT)
				- g1 * (Yl + Ye)*(ye(j,i).real() * CHWB + ye(j,i).imag() * CHWBT)
				+ gammaH * WC1(CeWR, i, j);
			l++;
			f[c + l * DF] = ((3. * cF2 - b02) * g22 +
				(-3. * Ye2 + 8. * YeYl - 3. * Yl2) * g12) * WC1(CeWI, i, j)
				+ g1g2 * (3. * Yl - Ye) * WC1(CeBI, i, j)
				- g2 * (ye(j,i).real() * CHWT - ye(j,i).imag() * CHW)
				- g1 * (Yl + Ye)*(ye(j,i).real() * CHWBT - ye(j,i).imag() * CHWB)
				+ gammaH * WC1(CeWI, i, j);
			l++;
			//CeB : Re and Im;
			f[c + l * DF] = (-3. * cF2 * g22 +
				(3. * Ye2 + 4. * YeYl + 3. * Yl2 - b01) * g12) * WC1(CeBR, i, j)
				+ 4. * cF2 * g1g2 * (3. * Yl - Ye) * WC1(CeWR, i, j)
				- 2. * g1 * (Yl + Ye)*(ye(j,i).real() * CHB + ye(j,i).imag() * CHBT)
				- 1.5 * g2 * (ye(j,i).real() * CHWB + ye(j,i).imag() * CHWBT)
				+ gammaH * WC1(CeBR, i, j);
			l++;
			f[c + l * DF] = (-3. * cF2 * g22 +
				(3. * Ye2 + 4. * YeYl + 3. * Yl2 - b01) * g12) * WC1(CeBI, i, j)
				+ 4. * cF2 * g1g2 * (3. * Yl - Ye) * WC1(CeWI, i, j)
				- 2. * g1 * (Yl + Ye)*(ye(j,i).real() * CHBT - ye(j,i).imag() * CHB)
				- 1.5 * g2 * (ye(j,i).real() * CHWBT - ye(j,i).imag() * CHWB)
				+ gammaH * WC1(CeBI, i, j);
			l++;
			//CuG: Re and Im
			f[c + l * DF] = ((10. * cF3 - 4. * cA3 - b03) * g32 - 3. * cF2 * g22
				+ (-3. * Yu2 + 8. * YuYq - 3. * Yq2) * g12) * WC1(CuGR, i, j)
				+ 8. * cF2 * g2g3 * WC1(CuWR, i, j) + 4. * g1g3 * (Yu + Yq) * WC1(CuBR, i, j)
				- 4. * g3 * (yu(j,i).real() * CHG + yu(j,i).imag() * CHGT)
				+ 3. * g32 * cA3 * (yu(j,i).real() * CG + yu(j,i).imag() * CGT)
				+ gammaH * WC1(CuGR, i, j);
			l++;
			f[c + l * DF ] = ((10. * cF3 - 4. * cA3 - b03) * g32 - 3 * cF2 * g22
				+ (-3. * Yu2 + 8. * YuYq - 3. * Yq2) * g12) * WC1(CuGI, i, j)
				+ 8. * cF2 * g2g3 * WC1(CuWI, i, j) + 4. * g1g3 * (Yu + Yq) * WC1(CuBI, i, j)
				- 4. * g3 * (yu(j,i).real() * CHGT - yu(j,i).imag() * CHG)
				+ 3. * g32 * cA3 * (yu(j,i).real() * CGT - yu(j,i).imag() * CG)
				+ gammaH * WC1(CuGI, i, j);
			l++;
			//CuW: Re and Im
			f[c + l * DF ] = (2. * cF3 * g32 + (3. * cF2 - b02) * g22
				+ (-3. * Yu2 + 8. * YuYq - 3. * Yq2) * g12) * WC1(CuWR, i, j)
				+ 2. * cF3 * g2g3 * WC1(CuGR, i, j) + g1g2 * (3. * Yq - Yu) * WC1(CuBR, i, j)
				- g2 * (yu(j,i).real() * CHW + yu(j,i).imag() * CHWT)
				+ g1 * (Yq + Yu) * (yu(j,i).real() * CHWB + yu(j,i).imag() * CHWBT)
				+ gammaH * WC1(CuWR, i, j);
			l++;
			f[c + l * DF ] = (2. * cF3 * g32 + (3. * cF2 - b02) * g22
				+ (-3. * Yu2 + 8. * YuYq - 3. * Yq2) * g12) * WC1(CuWI, i, j)
				+ 2. * cF3 * g2g3 * WC1(CuGI, i, j) + g1g2 * (3. * Yq - Yu) * WC1(CuBI, i, j)
				- g2 * (yu(j,i).real() * CHWT - yu(j,i).imag() * CHW)
				+ g1 * (Yq + Yu) * (yu(j,i).real() * CHWBT - yu(j,i).imag() * CHWB)
				+ gammaH * WC1(CuWI, i, j);
			l++;
			//CuB: Re and Im
			f[c + l * DF ] = (2. * cF3 * g32 - 3. * cF2 * g22 +
				(3. * Yu2 + 4. * YuYq + 3. * Yq2 - b01) * g12) * WC1(CuBR, i, j)
				+ 4. * cF2 * g1g2 * (3. * Yq - Yu) * WC1(CuWR, i, j)
				+ 4. * cF3 * g1g3 * (Yq + Yu) * WC1(CuGR, i, j)
				-(yu(j,i).real()*(2. * g1 * (Yq + Yu) * CHB - 1.5 * g2 * CHWB)
				+ yu(j,i).imag()*(2. * g1 * (Yq + Yu) * CHBT - 1.5 * g2 * CHWBT))
				+ gammaH * WC1(CuBR, i, j);
			l++;
			f[c + l * DF ] = (2. * cF3 * g32 - 3. * cF2 * g22 +
				(3. * Yu2 + 4. * YuYq + 3. * Yq2 - b01) * g12) * WC1(CuBI, i, j)
				+ 4. * cF2 * g1g2 * (3. * Yq - Yu) * WC1(CuWI, i, j)
				+ 4. * cF3 * g1g3 * (Yq + Yu) * WC1(CuGI, i, j)
				- (yu(j,i).real()*(2. * g1 * (Yq + Yu) * CHBT - 1.5 * g2 * CHWBT)
				- yu(j,i).imag()*(2. * g1 * (Yq + Yu) * CHB - 1.5 * g2 * CHWB))
				+ gammaH * WC1(CuBI, i, j);
			l++;


			//CdG : Re and Im;
			f[c + l * DF] = ((10. * cF3 - 4. * cA3 - b03) * g32 - 3. * cF2 * g22
				+ (-3. * Yd2 + 8. * YdYq - 3. * Yq2) * g12) * WC1(CdGR, i, j)
				+ 8. * cF2 * g2g3 * WC1(CdWR, i, j) + 4. * g1g3 * (Yd + Yq) * WC1(CdBR, i, j)
				- 4. * g3 * (yd(j,i).real() * CHG + yd(j,i).imag() * CHGT)
				+ 3. * g32 * cA3 * (yd(j,i).real() * CG + yd(j,i).imag() * CGT)
				+ gammaH * WC1(CdGR, i, j);
			l++;
			f[c + l * DF] = ((10. * cF3 - 4. * cA3 - b03) * g32 - 3. * cF2 * g22
				+ (-3. * Yd2 + 8. * YdYq - 3. * Yq2) * g12) * WC1(CdGI, i, j)
				+ 8. * cF2 * g2g3 * WC1(CdWI, i, j) + 4. * g1g3 * (Yd + Yq) * WC1(CdBI, i, j)
				- 4. * g3 * (yd(j,i).real() * CHGT - yd(j,i).imag() * CHG)
				+ 3. * g32 * cA3 * (yd(j,i).real() * CGT - yd(j,i).imag() * CG)
				+ gammaH * WC1(CdGI, i, j);

			l++;
			//CdW : Re and Im;
			f[c + l * DF] = (2. * cF3 * g32 + (3. * cF2 - b02) * g22
				+ (-3. * Yd2 + 8. * YdYq - 3. * Yq2) * g12) * WC1(CdWR, i, j)
				+ 2. * cF3 * g2g3 * WC1(CdGR, i, j) + g1g2 * (3. * Yq - Yd) * WC1(CdBR, i, j)
				- g2 * (yd(j,i).real() * CHW + yd(j,i).imag() * CHWT)
				- g1 * (Yq + Yd)*(yd(j,i).real() * CHWB + yd(j,i).imag() * CHWBT)
				+ gammaH * WC1(CdWR, i, j);
			l++;
			f[c + l * DF] = (2. * cF3 * g32 + (3. * cF2 - b02) * g22
				+ (-3. * Yd2 + 8. * YdYq - 3. * Yq2) * g12) * WC1(CdWI, i, j)
				+ 2. * cF3 * g2g3 * WC1(CdGI, i, j) + g1g2 * (3. * Yq - Yd) * WC1(CdBI, i, j)
				- g2 * (yd(j,i).real() * CHWT - yd(j,i).imag() * CHW)
				- g1 * (Yq + Yd)*(yd(j,i).real() * CHWBT - yd(j,i).imag() * CHWB)
				+ gammaH * WC1(CdWI, i, j);

			l++;
			//CdB : Re and Im;
			f[c + l * DF] = (2. * cF3 * g32 - 3. * cF2 * g22 +
				(3. * Yd2 + 4. * YdYq + 3. * Yq2 - b01) * g12) * WC1(CdBR, i, j)
				+ 4. * cF2 * g1g2 * (3. * Yq - Yd) * WC1(CdWR, i, j)
				+ 4. * cF3 * g1g3 * (Yq + Yd) * WC1(CdGR, i, j)
				- 2. * g1 * (Yq + Yd)*(yd(j,i).real() * CHB + yd(j,i).imag() * CHBT)
				- 1.5 * g2 * (yd(j,i).real() * CHWB + yd(j,i).imag() * CHWBT)
				+ gammaH * WC1(CdBR, i, j);
			l++;
			f[c + l * DF] = (2. * cF3 * g32 - 3. * cF2 * g22 +
				(3. * Yd2 + 4. * YdYq + 3. * Yq2 - b01) * g12) * WC1(CdBI, i, j)
				+ 4. * cF2 * g1g2 * (3. * Yq - Yd) * WC1(CdWI, i, j)
				+ 4. * cF3 * g1g3 * (Yq + Yd) * WC1(CdGI, i, j)
				- 2. * g1 * (Yq + Yd)*(yd(j,i).real() * CHBT - yd(j,i).imag() * CHB)
				- 1.5 * g2 * (yd(j,i).real() * CHWBT - yd(j,i).imag() * CHWB)
				+ gammaH * WC1(CdBI, i, j);

			l++;
			//Entries with 1 matrix product (1 summed index)

			for (a = 0; a < NG; a++) {
				l = 0;
				//CeW : Re and Im;
				f[c + l * DF] += 2. * (WC1(CeWR, i, a) * ge(a,j).real() - WC1(CeWI, i, a) * ge(a,j).imag())
					+(gl(i,a).real() * WC1(CeWR, a, j) - gl(i,a).imag() * WC1(CeWI, a, j));
				l++;
				f[c + l * DF] += 2. * (WC1(CeWR, i, a) * ge(a,j).imag() + WC1(CeWI, i, a) * ge(a,j).real())
					+(gl(i,a).real() * WC1(CeWI, a, j) + gl(i,a).imag() * WC1(CeWR, a, j));
				l++;

				//CeB : Re and Im;
				f[c + l * DF] += 2. * (WC1(CeBR, i, a) * ge(a,j).real() - WC1(CeBI, i, a) * ge(a,j).imag())
					+ 2. * (yedye(i,a).real() * WC1(CeBR, a, j) - yedye(i,a).imag() * WC1(CeBI, a, j))
					+((gl(i,a).real() * WC1(CeBR, a, j) - gl(i,a).imag() * WC1(CeBI, a, j)));
				l++;
				f[c + l * DF] += 2. * (WC1(CeBR, i, a) * ge(a,j).imag() + WC1(CeBI, i, a) * ge(a,j).real())
					+ 2. * (yedye(i,a).real() * WC1(CeBI, a, j) + yedye(i,a).imag() * WC1(CeBR, a, j))
					+((gl(i,a).real() * WC1(CeBI, a, j) + gl(i,a).imag() * WC1(CeBR, a, j)));
				l++;

				//CuG: Re and Im
				f[c + l * DF] += 2. * ((yudyu(i,a).real() - yddyd(i,a).real()) * WC1(CuGR, a, j)- (yudyu(i,a).imag() - yddyd(i,a).imag()) * WC1(CuGI, a, j))
					+ 2. * (WC1(CuGR, i, a) * gu(a,j).real() - WC1(CuGI, i, a) * gu(a,j).imag())
					+ (gq(i,a).real() * WC1(CuGR, a, j) - gq(i,a).imag() * WC1(CuGI, a, j))
					-(WC1(CdGR, i, a) * ydyud(a,j).real() - WC1(CdGI, i, a) * ydyud(a,j).imag());
				l++;
				f[c + l * DF] += 2. * ((yudyu(i,a).real() - yddyd(i,a).real()) * WC1(CuGI, a, j)+ (yudyu(i,a).imag() - yddyd(i,a).imag()) * WC1(CuGR, a, j))
					+ 2. * (WC1(CuGR, i, a) * gu(a,j).imag() + WC1(CuGI, i, a) * gu(a,j).real())
					+ (gq(i,a).real() * WC1(CuGI, a, j) + gq(i,a).imag() * WC1(CuGR, a, j))
					-(WC1(CdGI, i, a) * ydyud(a,j).real() + WC1(CdGR, i, a) * ydyud(a,j).imag());
				l++;

				//CuW: Re and Im
				f[c + l * DF] += 2. * (yddyd(i,a).real() * WC1(CuWR, a, j) - yddyd(i,a).imag() * WC1(CuWI, a, j))
					+ 2. * (WC1(CuWR, i, a) * gu(a,j).real() - WC1(CuWI, i, a) * gu(a,j).imag())
					+ (gq(i,a).real() * WC1(CuWR, a, j) - gq(i,a).imag() * WC1(CuWI, a, j))
					-(WC1(CdWR, i, a) * ydyud(a,j).real() - WC1(CdWI, i, a) * ydyud(a,j).imag());
				l++;

				f[c + l * DF] += 2. * (yddyd(i,a).imag() * WC1(CuWR, a, j) + yddyd(i,a).real() * WC1(CuWI, a, j))
					+ 2. * (WC1(CuWR, i, a) * gu(a,j).imag() + WC1(CuWI, i, a) * gu(a,j).real())
					+ (gq(i,a).imag() * WC1(CuWR, a, j) + gq(i,a).real() * WC1(CuWI, a, j))
					-(WC1(CdWR, i, a) * ydyud(a,j).imag() + WC1(CdWI, i, a) * ydyud(a,j).real());
				l++;

				//CuB: Re and Im
				f[c + l * DF] += 2. * ((yudyu(i,a).real() - yddyd(i,a).real()) * WC1(CuBR, a, j)- (yudyu(i,a).imag() - yddyd(i,a).imag()) * WC1(CuBI, a, j))
					+ 2. * (WC1(CuBR, i, a) * gu(a,j).real() - WC1(CuBI, i, a) * gu(a,j).imag())
					+ (gq(i,a).real() * WC1(CuBR, a, j) - gq(i,a).imag() * WC1(CuBI, a, j))
					-(WC1(CdBR, i, a) * ydyud(a,j).real() - WC1(CdBI, i, a) * ydyud(a,j).imag());
				l++;


				f[c + l * DF] += 2. * ((yudyu(i,a).real() - yddyd(i,a).real()) * WC1(CuBI, a, j) + (yudyu(i,a).imag() - yddyd(i,a).imag()) * WC1(CuBR, a, j))
					+ 2. * (WC1(CuBR, i, a) * gu(a,j).imag() + WC1(CuBI, i, a) * gu(a,j).real())
					+ (gq(i,a).real() * WC1(CuBI, a, j) + gq(i,a).imag() * WC1(CuBR, a, j))
					-(WC1(CdBR, i, a) * ydyud(a,j).imag() + WC1(CdBI, i, a) * ydyud(a,j).real());
				l++;

				//CdG: Re and Im 
				f[c + l * DF] += -2. * ((yudyu(i,a).real() - yddyd(i,a).real()) * WC1(CdGR, a, j)-(yudyu(i,a).imag() - yddyd(i,a).imag()) * WC1(CdGI, a, j))
					-(WC1(CuGR, i, a) * yuydd(a,j).real() - WC1(CuGI, i, a) * yuydd(a,j).imag()) +
					2. * (WC1(CdGR, i, a) * gd(a,j).real() - WC1(CdGI, i, a) * gd(a,j).imag())
					+(gq(i,a).real() * WC1(CdGR, a, j) - gq(i,a).imag() * WC1(CdGI, a, j));
				l++;
				f[c + l * DF] += -2. * ((yudyu(i,a).real() - yddyd(i,a).real()) * WC1(CdGI, a, j)+(yudyu(i,a).imag() - yddyd(i,a).imag()) * WC1(CdGR, a, j))
					-(WC1(CuGR, i, a) * yuydd(a,j).imag() + WC1(CuGI, i, a) * yuydd(a,j).real()) +
					2. * (WC1(CdGR, i, a) * gd(a,j).imag() + WC1(CdGI, i, a) * gd(a,j).real())
					+(gq(i,a).real() * WC1(CdGI, a, j) + gq(i,a).imag() * WC1(CdGR, a, j));
				l++;

				//CdW: Re and Im 
				f[c + l * DF] += 2. * (yudyu(i,a).real() * WC1(CdWR, a, j) - yudyu(i,a).imag() * WC1(CdWI, a, j))
					-(WC1(CuWR, i, a) * yuydd(a,j).real() - WC1(CuWI, i, a) * yuydd(a,j).imag())
					+ 2. * (WC1(CdWR, i, a) * gd(a,j).real() - WC1(CdWI, i, a) * gd(a,j).imag())
					+(gq(i,a).real() * WC1(CdWR, a, j) - gq(i,a).imag() * WC1(CdWI, a, j));

				l++;
				f[c + l * DF] += 2. * (yudyu(i,a).real() * WC1(CdWI, a, j) + yudyu(i,a).imag() * WC1(CdWR, a, j))
					-(WC1(CuWR, i, a) * yuydd(a,j).imag() + WC1(CuWI, i, a) * yuydd(a,j).real())
					+ 2. * (WC1(CdWR, i, a) * gd(a,j).imag() + WC1(CdWI, i, a) * gd(a,j).real())
					+(gq(i,a).real() * WC1(CdWI, a, j) + gq(i,a).imag() * WC1(CdWR, a, j));
				l++;

				// CdB: Re and Im 
				f[c + l * DF] += -2. * ((yudyu(i,a).real() - yddyd(i,a).real()) * WC1(CdBR, a, j)-(yudyu(i,a).imag() - yddyd(i,a).imag()) * WC1(CdBI, a, j))
					-(WC1(CuBR, i, a) * yuydd(a,j).real() - WC1(CuBI, i, a) * yuydd(a,j).imag())
					+ 2. * (WC1(CdBR, i, a) * gd(a,j).real() - WC1(CdBI, i, a) * gd(a,j).imag())
					+(gq(i,a).real() * WC1(CdBR, a, j) - gq(i,a).imag() * WC1(CdBI, a, j));
				l++;
				f[c + l * DF] += -2. * ((yudyu(i,a).imag() - yddyd(i,a).imag()) * WC1(CdBR, a, j)+(yudyu(i,a).real() - yddyd(i,a).real()) * WC1(CdBI, a, j))
					-(WC1(CuBR, i, a) * yuydd(a,j).imag() + WC1(CuBI, i, a) * yuydd(a,j).real())
					+ 2. * (WC1(CdBR, i, a) * gd(a,j).imag() + WC1(CdBI, i, a) * gd(a,j).real())
					+(gq(i,a).imag() * WC1(CdBR, a, j) + gq(i,a).real() * WC1(CdBI, a, j));
				l++;

				//Entries with 2 matrix products (2 summed indices) 
				for (b = 0; b < NG; b++) {
					l = 0;
					//CeW : Re and Im;
					f[c + l * DF] += -2. * g2 * NC * (WC5(Clequ3R, i, j, a, b) * yu(b,a).real()
						- WC5(Clequ3I, i, j, a, b) * yu(b,a).imag());
					l++;
					f[c + l * DF] += -2. * g2 * NC * (WC5(Clequ3R, i, j, a, b) * yu(b,a).imag()
						+ WC5(Clequ3I, i, j, a, b) * yu(b,a).real());
					l++;
					//CeB : Re and Im;
					f[c + l * DF] += 4. * g1 * NC * (Yu + Yq)*
						(WC5(Clequ3R, i, j, a, b) * yu(b,a).real()
						- WC5(Clequ3I, i, j, a, b) * yu(b,a).imag());
					l++;
					f[c + l * DF] += 4. * g1 * NC * (Yu + Yq)*
						(WC5(Clequ3R, i, j, a, b) * yu(b,a).imag()
						+ WC5(Clequ3I, i, j, a, b) * yu(b,a).real());
					l++;
					//CuG: Re and Im
					f[c + l * DF] += -g3 * (
						(WC5(Cquqd1R, a, j, i, b) -(0.5 / NC) * WC5(Cquqd8R, a, j, i, b)) * yd(b,a).real()
						- (WC5(Cquqd1I, a, j, i, b) -(0.5 / NC) * WC5(Cquqd8I, a, j, i, b)) * yd(b,a).imag());
					l++;
					f[c + l * DF] += -g3 * (
						(WC5(Cquqd1R, a, j, i, b) -(0.5 / NC) * WC5(Cquqd8R, a, j, i, b)) * yd(b,a).imag()
						+ (WC5(Cquqd1I, a, j, i, b) -(0.5 / NC) * WC5(Cquqd8I, a, j, i, b)) * yd(b,a).real());
					l++;
					//CuW: Re and Im
					f[c + l * DF] += -2. * g2 *
						(WC5(Clequ3R, a, b, i, j) * ye(b,a).real() - WC5(Clequ3I, a, b, i, j) * ye(b,a).imag())
						+ 0.25 * g2 * (
						(WC5(Cquqd1R, a, j, i, b) + cF3 * WC5(Cquqd8R, a, j, i, b)) * yd(b,a).real()
						- (WC5(Cquqd1I, a, j, i, b) + cF3 * WC5(Cquqd8I, a, j, i, b)) * yd(b,a).imag()
						);
					l++;
					f[c + l * DF] += -2. * g2 *
						(WC5(Clequ3R, a, b, i, j) * ye(b,a).imag() + WC5(Clequ3I, a, b, i, j) * ye(b,a).real())
						+ 0.25 * g2 * (
						(WC5(Cquqd1R, a, j, i, b) + cF3 * WC5(Cquqd8R, a, j, i, b)) * yd(b,a).imag()
						+ (WC5(Cquqd1I, a, j, i, b) + cF3 * WC5(Cquqd8I, a, j, i, b)) * yd(b,a).real()
						);
					l++;
					//CuB: Re and Im
					f[c + l * DF] += 4. * g1 * (Ye + Yl)*(
						WC5(Clequ3R, a, b, i, j) * ye(b,a).real() - WC5(Clequ3I, a, b, i, j) * ye(b,a).imag())
						- 0.5 * g1 * (Yd + Yq)*(
						(WC5(Cquqd1R, a, j, i, b) + cF3 * WC5(Cquqd8R, a, j, i, b)) * yd(b,a).real()
						- (WC5(Cquqd1I, a, j, i, b) + cF3 * WC5(Cquqd8I, a, j, i, b)) * yd(b,a).imag()
						);
					l++;
					f[c + l * DF] += 4. * g1 * (Ye + Yl)*(
						WC5(Clequ3R, a, b, i, j) * ye(b,a).imag() + WC5(Clequ3I, a, b, i, j) * ye(b,a).real())
						- 0.5 * g1 * (Yd + Yq)*(
						(WC5(Cquqd1R, a, j, i, b) + cF3 * WC5(Cquqd8R, a, j, i, b)) * yd(b,a).imag()
						+ (WC5(Cquqd1I, a, j, i, b) + cF3 * WC5(Cquqd8I, a, j, i, b)) * yd(b,a).real()
						);
					l++;
					//CdG: Re and Im 
					f[c + l * DF] += -g3 * (
						(WC5(Cquqd1R, i, b, a, j)-(0.5 / NC) * WC5(Cquqd8R, i, b, a, j)) * yu(b,a).real()
						-(WC5(Cquqd1I, i, b, a, j)-(0.5 / NC) * WC5(Cquqd8I, i, b, a, j)) * yu(b,a).imag());
					l++;
					f[c + l * DF] += -g3 * (
						(WC5(Cquqd1R, i, b, a, j)-(0.5 / NC) * WC5(Cquqd8R, i, b, a, j)) * yu(b,a).imag()
						+(WC5(Cquqd1I, i, b, a, j)-(0.5 / NC) * WC5(Cquqd8I, i, b, a, j)) * yu(b,a).real());
					l++;
					//CdW: Re and Im 
					f[c + l * DF] += 0.25 * g2 * (
						(WC5(Cquqd1R, i, b, a, j) + cF3 * WC5(Cquqd8R, i, b, a, j)) * yu(b,a).real()
						- (WC5(Cquqd1I, i, b, a, j) + cF3 * WC5(Cquqd8I, i, b, a, j)) * yu(b,a).imag()
						);
					l++;
					f[c + l * DF] += 0.25 * g2 * (
						(WC5(Cquqd1R, i, b, a, j) + cF3 * WC5(Cquqd8R, i, b, a, j)) * yu(b,a).imag()
						+ (WC5(Cquqd1I, i, b, a, j) + cF3 * WC5(Cquqd8I, i, b, a, j)) * yu(b,a).real()
						);
					l++;
					//CdB: Re and Im 
					f[c + l * DF] += -0.5 * g1 * (Yu + Yq)*(
						(WC5(Cquqd1R, i, b, a, j) + cF3 * WC5(Cquqd8R, i, b, a, j)) * yu(b,a).real()
						-(WC5(Cquqd1I, i, b, a, j) + cF3 * WC5(Cquqd8I, i, b, a, j)) * yu(b,a).imag()
						);
					l++;
					f[c + l * DF] += -0.5 * g1 * (Yu + Yq)*(
						(WC5(Cquqd1R, i, b, a, j) + cF3 * WC5(Cquqd8R, i, b, a, j)) * yu(b,a).imag()
						+(WC5(Cquqd1I, i, b, a, j) + cF3 * WC5(Cquqd8I, i, b, a, j)) * yu(b,a).real()
						);
					l++;
				}
			}

			for (l = 0; l < 2 * N6; l++) {
				f[c + l * DF] *= loop_factor;
			}

			c++;
		}
	}
	c += (N6 * 2 - 1) * DF;


	//----------------------RGE SMEFT class 7--------------------------
	{

		//CHl1R
		for (i = 0; i < DWC2R; i++) {
			p = WC2R_indices[i][0];
			r = WC2R_indices[i][1];
			//RGE 2
			f[c] = 0.;
			f[c] = -0.5 * yedye(p,r).real() *(CHBOX + CHD) + 2. * gammaH * WC2R(CHl1R, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 1.5 * (yedye(p,t).real() *(WC2R(CHl1R, t, r) + 3. * WC2R(CHl3R, t, r))
					- yedye(p,t).imag() *(WC2I(CHl1I, t, r) + 3. * WC2I(CHl3I, t, r)))
					+ 1.5 * ((WC2R(CHl1R, p, t) + 3. * WC2R(CHl3R, p, t)) * yedye(t,r).real()
					-(WC2I(CHl1I, p, t) + 3. * WC2I(CHl3I, p, t)) * yedye(t,r).imag())
					+(gl(p,t).real() * WC2R(CHl1R, t, r) - gl(p,t).imag() * WC2I(CHl1I, t, r))
					+(WC2R(CHl1R, p, t) * gl(t,r).real() - WC2I(CHl1I, p, t) * gl(t,r).imag());

				for (s = 0; s < NG; s++) {
					f[c] += -((yedag(p,s).real() * WC2R(CHeR, s, t) - yedag(p,s).imag() * WC2I(CHeI, s, t)) * ye(t,r).real()
						-(yedag(p,s).imag() * WC2R(CHeR, s, t) + yedag(p,s).real() * WC2I(CHeI, s, t)) * ye(t,r).imag())
						+ 2. * (WC7R(CleR, p, r, s, t) * ge(t,s).real() - WC7I(CleI, p, r, s, t) * ge(t,s).imag())
						+((-2. * WC6R(CllR, p, r, s, t) - 2. * WC6R(CllR, s, t, p, r)
						- WC6R(CllR, p, t, s, r) - WC6R(CllR, s, r, p, t)) * yedye(t,s).real()
						-(-2. * WC6I(CllI, p, r, s, t) - 2. * WC6I(CllI, s, t, p, r)
						- WC6I(CllI, p, t, s, r) - WC6R(CllI, s, r, p, t)) * yedye(t,s).imag())
						- 2. * NC * (WC7R(Clq1R, p, r, s, t)*(yddyd(t,s).real() - yudyu(t,s).real())
						- WC7I(Clq1I, p, r, s, t)*(yddyd(t,s).imag() - yudyu(t,s).imag()))
						- 2. * NC * (WC7R(CluR, p, r, s, t) * gu(t,s).real() - WC7I(CluI, p, r, s, t) * gu(t,s).imag())
						+ 2. * NC * (WC7R(CldR, p, r, s, t) * gd(t,s).real() - WC7I(CldI, p, r, s, t) * gd(t,s).imag());
				}
			}
			//Necessary since [RGE2] and [RGE3] use different indices.
			r = WC2R_indices[i][0];
			s = WC2R_indices[i][1];
			//RGE 3
			f[c] += 0.5 * xiB * g12 * delta[r][s] * Yl
				+ FOUR_THIRDS * g12 * Yh2 * WC2R(CHl1R, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += FOUR_THIRDS * g12 * Yh * NC * Yd * WC7R(CldR, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYe * WC7R(CleR, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYl * (
					2. * WC6R(CllR, r, s, w, w) + WC6R(CllR, r, w, w, s)
					+ WC6R(CllR, w, s, r, w) + 2. * WC6R(CllR, w, w, r, s))
					+ EIGHT_THIRDS * g12 * NC * YhYq * WC7R(Clq1R, r, s, w, w)
					+ FOUR_THIRDS * g12 * NC * YhYu * WC7R(CluR, r, s, w, w);
			}

			f[c] *= loop_factor;
			c++;
		}
		//CHl1I
		for (i = 0; i < DWC2I; i++) {
			p = WC2I_indices[i][0];
			r = WC2I_indices[i][1];
			//RGE 2
			f[c] = -0.5 * yedye(p,r).imag() *(CHBOX + CHD) + 2. * gammaH * WC2I(CHl1I, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 1.5 * (yedye(p,t).imag() *(WC2R(CHl1R, t, r) + 3. * WC2R(CHl3R, t, r))
					+ yedye(p,t).real() *(WC2I(CHl1I, t, r) + 3. * WC2I(CHl3I, t, r)))
					+ 1.5 * ((WC2R(CHl1R, p, t) + 3. * WC2R(CHl3R, p, t)) * yedye(t,r).imag()
					+(WC2I(CHl1I, p, t) + 3. * WC2I(CHl3I, p, t)) * yedye(t,r).real())
					+(gl(p,t).imag() * WC2R(CHl1R, t, r) + gl(p,t).real() * WC2I(CHl1I, t, r))
					+(WC2R(CHl1R, p, t) * gl(t,r).imag() + WC2I(CHl1I, p, t) * gl(t,r).real());

				for (s = 0; s < NG; s++) {
					f[c] += -((yedag(p,s).real() * WC2R(CHeR, s, t) - yedag(p,s).imag() * WC2I(CHeI, s, t)) * ye(t,r).imag()
						+(yedag(p,s).imag() * WC2R(CHeR, s, t) + yedag(p,s).real() * WC2I(CHeI, s, t)) * ye(t,r).real())
						+ 2. * (WC7R(CleR, p, r, s, t) * ge(t,s).imag() + WC7I(CleI, p, r, s, t) * ge(t,s).real())
						+((-2. * WC6R(CllR, p, r, s, t) - 2. * WC6R(CllR, s, t, p, r)
						- WC6R(CllR, p, t, s, r) - WC6R(CllR, s, r, p, t)) * yedye(t,s).imag()
						+(-2. * WC6I(CllI, p, r, s, t) - 2. * WC6I(CllI, s, t, p, r)
						- WC6I(CllI, p, t, s, r) - WC6R(CllI, s, r, p, t)) * yedye(t,s).real())
						- 2. * NC * (WC7R(Clq1R, p, r, s, t)*(yddyd(t,s).imag() - yudyu(t,s).imag())
						+ WC7I(Clq1I, p, r, s, t)*(yddyd(t,s).real() - yudyu(t,s).real()))
						- 2. * NC * (WC7R(CluR, p, r, s, t) * gu(t,s).imag() + WC7I(CluI, p, r, s, t) * gu(t,s).real())
						+ 2. * NC * (WC7R(CldR, p, r, s, t) * gd(t,s).imag() + WC7I(CldI, p, r, s, t) * gd(t,s).real());
				}
			}

			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2I_indices[i][0];
			s = WC2I_indices[i][1];
			//RGE 3
			f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHl1I, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += FOUR_THIRDS * g12 * Yh * NC * Yd * WC7I(CldI, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYe * WC7I(CleI, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYl * (
					2. * WC6I(CllI, r, s, w, w) + WC6I(CllI, r, w, w, s)
					+ WC6I(CllI, w, s, r, w) + 2. * WC6I(CllI, w, w, r, s))
					+ EIGHT_THIRDS * g12 * NC * YhYq * WC7I(Clq1I, r, s, w, w)
					+ FOUR_THIRDS * g12 * NC * YhYu * WC7I(CluI, r, s, w, w);
			}


			f[c] *= loop_factor;
			c++;
		}



		//CHl3R
		for (i = 0; i < DWC2R; i++) {
			p = WC2R_indices[i][0];
			r = WC2R_indices[i][1];
			//RGE 2
			f[c] = -0.5 * yedye(p,r).real() * CHBOX + 2. * gammaH * WC2R(CHl3R, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 0.5 * (yedye(p,t).real()*(3. * WC2R(CHl1R, t, r) + WC2R(CHl3R, t, r))
					- yedye(p,t).imag()*(3. * WC2I(CHl1I, t, r) + WC2I(CHl3I, t, r)))
					+ 0.5 * ((3. * WC2R(CHl1R, p, t) + WC2R(CHl3R, p, t)) * yedye(t,r).real()
					-(3. * WC2I(CHl1I, p, t) + WC2I(CHl3I, p, t)) * yedye(t,r).imag())
					+(gl(p,t).real() * WC2R(CHl3R, t, r) - gl(p,t).imag() * WC2I(CHl3I, t, r))
					+(WC2R(CHl3R, p, t) * gl(t,r).real() - WC2I(CHl3I, p, t) * gl(t,r).imag());
				for (s = 0; s < NG; s++) {
					f[c] += -((WC6R(CllR, p, t, s, r) + WC6R(CllR, s, r, p, t)) * yedye(t,s).real()
						-(WC6I(CllI, p, t, s, r) + WC6I(CllI, s, r, p, t)) * yedye(t,s).imag())
						- 2. * NC * (WC7R(Clq3R, p, r, s, t)*(yddyd(t,s).real() + yudyu(t,s).real())
						- WC7I(Clq3I, p, r, s, t)*(yddyd(t,s).imag() + yudyu(t,s).imag()));
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2R_indices[i][0];
			s = WC2R_indices[i][1];
			//RGE 3
			f[c] += ONE_SIXTH * g22 * CHBOX * delta[r][s] + ONE_THIRD * g22 * WC2R(CHl3R, r, s)
				- 6. * g22 * WC2R(CHl3R, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += TWO_THIRDS * g22 * delta[r][s]*
					(WC2R(CHl3R, w, w) + NC * WC2R(CHq3R, w, w))
					+ ONE_THIRD * g22 * (WC6R(CllR, r, w, w, s) + WC6R(CllR, w, s, r, w))
					+ TWO_THIRDS * g22 * NC * WC7R(Clq3R, r, s, w, w);
			}
			f[c] *= loop_factor;
			c++;
		}


		//CHl3I
		for (i = 0; i < DWC2I; i++) {
			p = WC2I_indices[i][0];
			r = WC2I_indices[i][1];
			//RGE 2
			f[c] = -0.5 * yedye(p,r).imag() * CHBOX + 2. * gammaH * WC2I(CHl3I, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 0.5 * (yedye(p,t).imag()*(3. * WC2R(CHl1R, t, r) + WC2R(CHl3R, t, r))
					+ yedye(p,t).real()*(3. * WC2I(CHl1I, t, r) + WC2I(CHl3I, t, r)))
					+ 0.5 * ((3. * WC2R(CHl1R, p, t) + WC2R(CHl3R, p, t)) * yedye(t,r).imag()
					+(3. * WC2I(CHl1I, p, t) + WC2I(CHl3I, p, t)) * yedye(t,r).real())
					+(gl(p,t).imag() * WC2R(CHl3R, t, r) + gl(p,t).real() * WC2I(CHl3I, t, r))
					+(WC2R(CHl3R, p, t) * gl(t,r).imag() + WC2I(CHl3I, p, t) * gl(t,r).real());
				for (s = 0; s < NG; s++) {
					f[c] += -((WC6R(CllR, p, t, s, r) + WC6R(CllR, s, r, p, t)) * yedye(t,s).imag()
						+(WC6I(CllI, p, t, s, r) + WC6I(CllI, s, r, p, t)) * yedye(t,s).real())
						- 2. * NC * (WC7R(Clq3R, p, r, s, t)*(yddyd(t,s).imag() + yudyu(t,s).imag())
						+ WC7I(Clq3I, p, r, s, t)*(yddyd(t,s).real() + yudyu(t,s).real()));
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2I_indices[i][0];
			s = WC2I_indices[i][1];
			//RGE 3  
			f[c] += ONE_THIRD * g22 * WC2I(CHl3I, r, s)
				- 6. * g22 * WC2I(CHl3I, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += ONE_THIRD * g22 * (WC6I(CllI, r, w, w, s) + WC6I(CllI, w, s, r, w))
					+ TWO_THIRDS * g22 * NC * WC7I(Clq3I, r, s, w, w);
			}
			f[c] *= loop_factor;
			c++;
		}


		//CHeR
		for (i = 0; i < DWC2R; i++) {
			p = WC2R_indices[i][0];
			r = WC2R_indices[i][1];
			//RGE 2
			f[c] = ge(p,r).real()*(CHBOX + CHD) + 2. * gammaH * WC2R(CHeR, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 4. * (ge(p,t).real() * WC2R(CHeR, t, r) - ge(p,t).imag() * WC2I(CHeI, t, r))
					+ 4. * (WC2R(CHeR, p, t) * ge(t,r).real() - WC2I(CHeI, p, t) * ge(t,r).imag());

				for (s = 0; s < NG; s++) {
					f[c] += -2. * ((ye(p,s).real() * WC2R(CHl1R, s, t) - ye(p,s).imag() * WC2I(CHl1I, s, t)) * yedag(t,r).real()
						- (ye(p,s).imag() * WC2R(CHl1R, s, t) + ye(p,s).real() * WC2I(CHl1I, s, t)) * yedag(t,r).imag())
						- 2. * (WC7R(CleR, s, t, p, r) * yedye(t,s).real() - WC7I(CleI, s, t, p, r) * yedye(t,s).imag())
						+ 2. * ((WC8R(CeeR, p, r, s, t) + WC8R(CeeR, s, t, p, r)
						+ WC8R(CeeR, p, t, s, r) + WC8R(CeeR, s, r, p, t)) * ge(t,s).real()
						-(WC8I(CeeI, p, r, s, t) + WC8I(CeeI, s, t, p, r)
						+ WC8I(CeeI, p, t, s, r) + WC8I(CeeI, s, r, p, t)) * ge(t,s).imag())
						- 2. * NC * (WC7R(CeuR, p, r, s, t) * gu(t,s).real() - WC7I(CeuI, p, r, s, t) * gu(t,s).imag())
						+ 2. * NC * (WC7R(CedR, p, r, s, t) * gd(t,s).real() - WC7I(CedI, p, r, s, t) * gd(t,s).imag())
						- 2. * NC * (WC7R(CqeR, s, t, p, r)*(yddyd(t,s).real() - yudyu(t,s).real())
						- WC7I(CqeI, s, t, p, r)*(yddyd(t,s).imag() - yudyu(t,s).imag()));
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2R_indices[i][0];
			s = WC2R_indices[i][1];
			//RGE 3
			f[c] += 0.5 * xiB * g12 * delta[r][s] * Ye
				+ FOUR_THIRDS * g12 * Yh2 * WC2R(CHeR, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += FOUR_THIRDS * NC * g12 * YhYd * WC7R(CedR, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYe *
					(WC8R(CeeR, r, s, w, w) + WC8R(CeeR, r, w, w, s)
					+ WC8R(CeeR, w, s, r, w) + WC8R(CeeR, w, w, r, s))
					+ FOUR_THIRDS * g12 * NC * YhYu * WC7R(CeuR, r, s, w, w)
					+ EIGHT_THIRDS * g12 * YhYl * WC7R(CleR, w, w, r, s)
					+ EIGHT_THIRDS * g12 * NC * YhYq * WC7R(CqeR, w, w, r, s);
			}
			f[c] *= loop_factor;
			c++;
		}


		//CHeI
		for (i = 0; i < DWC2I; i++) {
			p = WC2I_indices[i][0];
			r = WC2I_indices[i][1];
			//RGE 2
			f[c] = ge(p,r).imag()*(CHBOX + CHD) + 2. * gammaH * WC2I(CHeI, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 4. * (ge(p,t).imag() * WC2R(CHeR, t, r) + ge(p,t).real() * WC2I(CHeI, t, r))
					+ 4. * (WC2R(CHeR, p, t) * ge(t,r).imag() + WC2I(CHeI, p, t) * ge(t,r).real());

				for (s = 0; s < NG; s++) {
					f[c] += -2. * ((ye(p,s).real() * WC2R(CHl1R, s, t) - ye(p,s).imag() * WC2I(CHl1I, s, t)) * yedag(t,r).imag()
						+(ye(p,s).imag() * WC2R(CHl1R, s, t) + ye(p,s).real() * WC2I(CHl1I, s, t)) * yedag(t,r).real())
						- 2. * (WC7R(CleR, s, t, p, r) * yedye(t,s).imag() + WC7I(CleI, s, t, p, r) * yedye(t,s).real())
						+ 2. * ((WC8R(CeeR, p, r, s, t) + WC8R(CeeR, s, t, p, r)
						+ WC8R(CeeR, p, t, s, r) + WC8R(CeeR, s, r, p, t)) * ge(t,s).imag()
						+(WC8I(CeeI, p, r, s, t) + WC8I(CeeI, s, t, p, r)
						+ WC8I(CeeI, p, t, s, r) + WC8I(CeeI, s, r, p, t)) * ge(t,s).real())
						- 2. * NC * (WC7R(CeuR, p, r, s, t) * gu(t,s).imag() + WC7I(CeuI, p, r, s, t) * gu(t,s).real())
						+ 2. * NC * (WC7R(CedR, p, r, s, t) * gd(t,s).imag() + WC7I(CedI, p, r, s, t) * gd(t,s).real())
						- 2. * NC * (WC7R(CqeR, s, t, p, r)*(yddyd(t,s).imag() - yudyu(t,s).imag())
						+ WC7I(CqeI, s, t, p, r)*(yddyd(t,s).real() - yudyu(t,s).real()));
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2I_indices[i][0];
			s = WC2I_indices[i][1];
			//RGE 3
			f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHeI, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += FOUR_THIRDS * NC * g12 * YhYd * WC7I(CedI, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYe *
					(WC8I(CeeI, r, s, w, w) + WC8I(CeeI, r, w, w, s)
					+ WC8I(CeeI, w, s, r, w) + WC8I(CeeI, w, w, r, s))
					+ FOUR_THIRDS * g12 * NC * YhYu * WC7I(CeuI, r, s, w, w)
					+ EIGHT_THIRDS * g12 * YhYl * WC7I(CleI, w, w, r, s)
					+ EIGHT_THIRDS * g12 * NC * YhYq * WC7I(CqeI, w, w, r, s);
			}
			f[c] *= loop_factor;
			c++;
		}

		//CHq1R
		for (i = 0; i < DWC2R; i++) {
			p = WC2R_indices[i][0];
			r = WC2R_indices[i][1];
			//RGE 2
			f[c] = 0.5 * (yudyu(p,r).real() - yddyd(p,r).real())*(CHBOX + CHD) + 2. * gammaH * WC2R(CHq1R, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 1.5 * ((yddyd(p,t).real() + yudyu(p,t).real()) * WC2R(CHq1R, t, r)
					-(yddyd(p,t).imag() + yudyu(p,t).imag()) * WC2I(CHq1I, t, r))
					+ 1.5 * (WC2R(CHq1R, p, t)*(yddyd(t,r).real() + yudyu(t,r).real())
					- WC2I(CHq1I, p, t)*(yddyd(t,r).imag() + yudyu(t,r).imag()))
					+ 4.5 * ((yddyd(p,t).real() - yudyu(p,t).real()) * WC2R(CHq3R, t, r)
					-(yddyd(p,t).imag() - yudyu(p,t).imag()) * WC2I(CHq3I, t, r))
					+ 4.5 * (WC2R(CHq3R, p, t)*(yddyd(t,r).real() - yudyu(t,r).real())
					- WC2I(CHq3I, p, t)*(yddyd(t,r).imag() - yudyu(t,r).imag()))
					+(gq(p,t).real() * WC2R(CHq1R, t, r) - gq(p,t).imag() * WC2I(CHq1I, t, r))
					+(WC2R(CHq1R, p, t) * gq(t,r).real() - WC2I(CHq1I, p, t) * gq(t,r).imag())
					;
				for (s = 0; s < NG; s++) {
					f[c] += -((yudag(p,s).real() * WC2R(CHuR, s, t) - yudag(p,s).imag() * WC2I(CHuI, s, t)) * yu(t,r).real()
						-(yudag(p,s).imag() * WC2R(CHuR, s, t) + yudag(p,s).real() * WC2I(CHuI, s, t)) * yu(t,r).imag())
						-((yddag(p,s).real() * WC2R(CHdR, s, t) - yddag(p,s).imag() * WC2I(CHdI, s, t)) * yd(t,r).real()
						-(yddag(p,s).imag() * WC2R(CHdR, s, t) + yddag(p,s).real() * WC2I(CHdI, s, t)) * yd(t,r).imag())
						+ 2. * (WC7R(CqeR, p, r, s, t) * ge(t,s).real() - WC7I(CqeI, p, r, s, t) * ge(t,s).imag())
						- 2. * (WC7R(Clq1R, s, t, p, r) * yedye(t,s).real() - WC7I(Clq1I, s, t, p, r) * yedye(t,s).imag())
						-(
						(2. * NC * WC6R(Cqq1R, p, r, s, t) + 2. * NC * WC6R(Cqq1R, s, t, p, r)
						+ WC6R(Cqq1R, p, t, s, r) + WC6R(Cqq1R, s, r, p, t)
						+ 3. * WC6R(Cqq3R, p, t, s, r) + 3. * WC6R(Cqq3R, s, r, p, t))*(yddyd(t,s).real() - yudyu(t,s).real())
						-(2. * NC * WC6I(Cqq1I, p, r, s, t) + 2. * NC * WC6I(Cqq1I, s, t, p, r)
						+ WC6I(Cqq1I, p, t, s, r) + WC6I(Cqq1I, s, r, p, t)
						+ 3. * WC6I(Cqq3I, p, t, s, r) + 3. * WC6I(Cqq3I, s, r, p, t))*(yddyd(t,s).imag() - yudyu(t,s).imag())
						)
						- 2. * NC * (WC7R(Cqu1R, p, r, s, t) * gu(t,s).real() - WC7I(Cqu1I, p, r, s, t) * gu(t,s).imag())
						+ 2. * NC * (WC7R(Cqd1R, p, r, s, t) * gd(t,s).real() - WC7I(Cqd1I, p, r, s, t) * gd(t,s).imag());

				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2R_indices[i][0];
			s = WC2R_indices[i][1];
			f[c] += 0.5 * xiB * g12 * delta[r][s] * Yq
				+ FOUR_THIRDS * g12 * Yh2 * WC2R(CHq1R, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += EIGHT_THIRDS * g12 * YhYl * WC7R(Clq1R, w, w, r, s)
					+ FOUR_THIRDS * g12 * NC * YhYd * WC7R(Cqd1R, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYe * WC7R(CqeR, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYq * (
					2. * NC * WC6R(Cqq1R, r, s, w, w) + WC6R(Cqq1R, r, w, w, s)
					+ WC6R(Cqq1R, w, s, r, w) + 2. * NC * WC6R(Cqq1R, w, w, r, s))
					+ 4. * g12 * YhYq *
					(WC6R(Cqq3R, r, w, w, s) + WC6R(Cqq3R, w, s, r, w))
					+ FOUR_THIRDS * g12 * NC * YhYu * WC7R(Cqu1R, r, s, w, w);
			}

			f[c] *= loop_factor;
			c++;
		}



		//CHq1I
		for (i = 0; i < DWC2I; i++) {
			p = WC2I_indices[i][0];
			r = WC2I_indices[i][1];
			//RGE 2
			f[c] = 0.5 * (yudyu(p,r).imag() - yddyd(p,r).imag())*(CHBOX + CHD) + 2. * gammaH * WC2I(CHq1I, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 1.5 * ((yddyd(p,t).imag() + yudyu(p,t).imag()) * WC2R(CHq1R, t, r)
					+(yddyd(p,t).real() + yudyu(p,t).real()) * WC2I(CHq1I, t, r))
					+ 1.5 * (WC2R(CHq1R, p, t)*(yddyd(t,r).imag() + yudyu(t,r).imag())
					+ WC2I(CHq1I, p, t)*(yddyd(t,r).real() + yudyu(t,r).real()))
					+ 4.5 * ((yddyd(p,t).imag() - yudyu(p,t).imag()) * WC2R(CHq3R, t, r)
					+(yddyd(p,t).real() - yudyu(p,t).real()) * WC2I(CHq3I, t, r))
					+ 4.5 * (WC2R(CHq3R, p, t)*(yddyd(t,r).imag() - yudyu(t,r).imag())
					+ WC2I(CHq3I, p, t)*(yddyd(t,r).real() - yudyu(t,r).real()))
					+(gq(p,t).imag() * WC2R(CHq1R, t, r) + gq(p,t).real() * WC2I(CHq1I, t, r))
					+(WC2R(CHq1R, p, t) * gq(t,r).imag() + WC2I(CHq1I, p, t) * gq(t,r).real())
					;
				for (s = 0; s < NG; s++) {
					f[c] += -((yudag(p,s).real() * WC2R(CHuR, s, t) - yudag(p,s).imag() * WC2I(CHuI, s, t)) * yu(t,r).imag()
						+(yudag(p,s).imag() * WC2R(CHuR, s, t) + yudag(p,s).real() * WC2I(CHuI, s, t)) * yu(t,r).real())
						-((yddag(p,s).real() * WC2R(CHdR, s, t) - yddag(p,s).imag() * WC2I(CHdI, s, t)) * yd(t,r).imag()
						+(yddag(p,s).imag() * WC2R(CHdR, s, t) + yddag(p,s).real() * WC2I(CHdI, s, t)) * yd(t,r).real())
						+ 2. * (WC7R(CqeR, p, r, s, t) * ge(t,s).imag() + WC7I(CqeI, p, r, s, t) * ge(t,s).real())
						- 2. * (WC7R(Clq1R, s, t, p, r) * yedye(t,s).imag() + WC7I(Clq1I, s, t, p, r) * yedye(t,s).real())
						-(
						(2. * NC * WC6R(Cqq1R, p, r, s, t) + 2. * NC * WC6R(Cqq1R, s, t, p, r)
						+ WC6R(Cqq1R, p, t, s, r) + WC6R(Cqq1R, s, r, p, t)
						+ 3. * WC6R(Cqq3R, p, t, s, r) + 3. * WC6R(Cqq3R, s, r, p, t))*(yddyd(t,s).imag() - yudyu(t,s).imag())
						+(2. * NC * WC6I(Cqq1I, p, r, s, t) + 2. * NC * WC6I(Cqq1I, s, t, p, r)
						+ WC6I(Cqq1I, p, t, s, r) + WC6I(Cqq1I, s, r, p, t)
						+ 3. * WC6I(Cqq3I, p, t, s, r) + 3. * WC6I(Cqq3I, s, r, p, t))*(yddyd(t,s).real() - yudyu(t,s).real())
						)
						- 2. * NC * (WC7R(Cqu1R, p, r, s, t) * gu(t,s).imag() + WC7I(Cqu1I, p, r, s, t) * gu(t,s).real())
						+ 2. * NC * (WC7R(Cqd1R, p, r, s, t) * gd(t,s).imag() + WC7I(Cqd1I, p, r, s, t) * gd(t,s).real());

				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2I_indices[i][0];
			s = WC2I_indices[i][1];
			f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHq1I, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += EIGHT_THIRDS * g12 * YhYl * WC7I(Clq1I, w, w, r, s)
					+ FOUR_THIRDS * g12 * NC * YhYd * WC7I(Cqd1I, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYe * WC7I(CqeI, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYq * (
					2. * NC * WC6I(Cqq1I, r, s, w, w) + WC6I(Cqq1I, r, w, w, s)
					+ WC6I(Cqq1I, w, s, r, w) + 2. * NC * WC6I(Cqq1I, w, w, r, s))
					+ 4. * g12 * YhYq *
					(WC6I(Cqq3I, r, w, w, s) + WC6I(Cqq3I, w, s, r, w))
					+ FOUR_THIRDS * g12 * NC * YhYu * WC7I(Cqu1I, r, s, w, w);
			}
			f[c] *= loop_factor;
			c++;
		}


		//CHq3R 
		for (i = 0; i < DWC2R; i++) {
			p = WC2R_indices[i][0];
			r = WC2R_indices[i][1];
			//RGE 2
			f[c] = -0.5 * (yudyu(p,r).real() + yddyd(p,r).real()) * CHBOX + 2. * gammaH * WC2R(CHq3R, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 1.5 * ((yddyd(p,t).real() - yudyu(p,t).real()) * WC2R(CHq1R, t, r)
					-(yddyd(p,t).imag() - yudyu(p,t).imag()) * WC2I(CHq1I, t, r))
					+ 1.5 * (WC2R(CHq1R, p, t)*(yddyd(t,r).real() - yudyu(t,r).real())
					- WC2I(CHq1I, p, t)*(yddyd(t,r).imag() - yudyu(t,r).imag()))

					+ 0.5 * ((yddyd(p,t).real() + yudyu(p,t).real()) * WC2R(CHq3R, t, r)
					-(yddyd(p,t).imag() + yudyu(p,t).imag()) * WC2I(CHq3I, t, r))
					+ 0.5 * (WC2R(CHq3R, p, t)*(yddyd(t,r).real() + yudyu(t,r).real())
					- WC2I(CHq3I, p, t)*(yddyd(t,r).imag() + yudyu(t,r).imag()))

					+(gq(p,t).real() * WC2R(CHq3R, t, r) - gq(p,t).imag() * WC2I(CHq3I, t, r))
					+(WC2R(CHq3R, p, t) * gq(t,r).real() - WC2I(CHq3I, p, t) * gq(t,r).imag());


				for (s = 0; s < NG; s++) {
					f[c] += -(
						(2. * NC * WC6R(Cqq3R, p, r, s, t) + 2. * NC * WC6R(Cqq3R, s, t, p, r)
						+ WC6R(Cqq1R, p, t, s, r) + WC6R(Cqq1R, s, r, p, t)
						- WC6R(Cqq3R, p, t, s, r) - WC6R(Cqq3R, s, r, p, t))*(yddyd(t,s).real() + yudyu(t,s).real())
						-(2. * NC * WC6I(Cqq3I, p, r, s, t) + 2. * NC * WC6I(Cqq3I, s, t, p, r)
						+ WC6I(Cqq1I, p, t, s, r) + WC6I(Cqq1I, s, r, p, t)
						- WC6I(Cqq3I, p, t, s, r) - WC6I(Cqq3I, s, r, p, t))* (yddyd(t,s).imag() + yudyu(t,s).imag())
						)
						- 2. * (WC7R(Clq3R, s, t, p, r) * yedye(t,s).real() - WC7I(Clq3I, s, t, p, r) * yedye(t,s).imag());
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2R_indices[i][0];
			s = WC2R_indices[i][1];
			//RGE 3
			f[c] += ONE_SIXTH * g22 * CHBOX * delta[r][s]
				+ ONE_THIRD * g22 * WC2R(CHq3R, r, s) - 6. * g22 * WC2R(CHq3R, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += TWO_THIRDS * g22 * delta[r][s]*
					(WC2R(CHl3R, w, w) + NC * WC2R(CHq3R, w, w))
					+ TWO_THIRDS * g22 * WC7R(Clq3R, w, w, r, s)
					+ ONE_THIRD * g22 * (WC6R(Cqq1R, r, w, w, s) + WC6R(Cqq1R, w, s, r, w))
					+ ONE_THIRD * g22 *
					(2. * NC * WC6R(Cqq3R, r, s, w, w) - WC6R(Cqq3R, r, w, w, s)
					- WC6R(Cqq3R, w, s, r, w) + 2. * NC * WC6R(Cqq3R, w, w, r, s));
			}
			f[c] *= loop_factor;
			c++;
		}

		//CHq3I
		for (i = 0; i < DWC2I; i++) {
			p = WC2I_indices[i][0];
			r = WC2I_indices[i][1];
			//RGE 2
			f[c] = -0.5 * (yudyu(p,r).imag() + yddyd(p,r).imag()) * CHBOX + 2. * gammaH * WC2I(CHq3I, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 1.5 * ((yddyd(p,t).imag() - yudyu(p,t).imag()) * WC2R(CHq1R, t, r)
					+(yddyd(p,t).real() - yudyu(p,t).real()) * WC2I(CHq1I, t, r))
					+ 1.5 * (WC2R(CHq1R, p, t)*(yddyd(t,r).imag() - yudyu(t,r).imag())
					+ WC2I(CHq1I, p, t)*(yddyd(t,r).real() - yudyu(t,r).real()))

					+ 0.5 * ((yddyd(p,t).imag() + yudyu(p,t).imag()) * WC2R(CHq3R, t, r)
					+(yddyd(p,t).real() + yudyu(p,t).real()) * WC2I(CHq3I, t, r))
					+ 0.5 * (WC2R(CHq3R, p, t)*(yddyd(t,r).imag() + yudyu(t,r).imag())
					+ WC2I(CHq3I, p, t)*(yddyd(t,r).real() + yudyu(t,r).real()))

					+(gq(p,t).imag() * WC2R(CHq3R, t, r) + gq(p,t).real() * WC2I(CHq3I, t, r))
					+(WC2R(CHq3R, p, t) * gq(t,r).imag() + WC2I(CHq3I, p, t) * gq(t,r).real());


				for (s = 0; s < NG; s++) {
					f[c] += -(
						(2. * NC * WC6R(Cqq3R, p, r, s, t) + 2. * NC * WC6R(Cqq3R, s, t, p, r)
						+ WC6R(Cqq1R, p, t, s, r) + WC6R(Cqq1R, s, r, p, t)
						- WC6R(Cqq3R, p, t, s, r) - WC6R(Cqq3R, s, r, p, t))*(yddyd(t,s).imag() + yudyu(t,s).imag())
						+(2. * NC * WC6I(Cqq3I, p, r, s, t) + 2. * NC * WC6I(Cqq3I, s, t, p, r)
						+ WC6I(Cqq1I, p, t, s, r) + WC6I(Cqq1I, s, r, p, t)
						- WC6I(Cqq3I, p, t, s, r) - WC6I(Cqq3I, s, r, p, t))* (yddyd(t,s).real() + yudyu(t,s).real())
						)
						- 2. * (WC7R(Clq3R, s, t, p, r) * yedye(t,s).imag() + WC7I(Clq3I, s, t, p, r) * yedye(t,s).real());
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2I_indices[i][0];
			s = WC2I_indices[i][1];
			//RGE 3
			f[c] += ONE_THIRD * g22 * WC2I(CHq3I, r, s) - 6. * g22 * WC2I(CHq3I, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += TWO_THIRDS * g22 * WC7I(Clq3I, w, w, r, s)
					+ ONE_THIRD * g22 *
					(WC6I(Cqq1I, r, w, w, s) + WC6I(Cqq1I, w, s, r, w))
					+ ONE_THIRD * g22 *
					(2. * NC * WC6I(Cqq3I, r, s, w, w) - WC6I(Cqq3I, r, w, w, s)
					- WC6I(Cqq3I, w, s, r, w) + 2. * NC * WC6I(Cqq3I, w, w, r, s));
			}
			f[c] *= loop_factor;
			c++;
		}


		//CHuR
		for (i = 0; i < DWC2R; i++) {
			p = WC2R_indices[i][0];
			r = WC2R_indices[i][1];
			//RGE 2
			f[c] = -gu(p,r).real()*(CHBOX + CHD) + 2. * gammaH * WC2R(CHuR, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 4. * (gu(p,t).real() * WC2R(CHuR, t, r) - gu(p,t).imag() * WC2I(CHuI, t, r))
					+ 4. * (WC2R(CHuR, p, t) * gu(t,r).real() - WC2I(CHuI, p, t) * gu(t,r).imag())
					+(yuydd(p,t).real() * WC1(CHudR, r, t) + yuydd(p,t).imag() * WC1(CHudI, r, t))
					+(WC1(CHudR, p, t) * ydyud(t,r).real() - WC1(CHudI, p, t) * ydyud(t,r).imag())
					;
				for (s = 0; s < NG; s++) {
					f[c] += -2. * ((yu(p,s).real() * WC2R(CHq1R, s, t) - yu(p,s).imag() * WC2I(CHq1I, s, t)) * yudag(t,r).real()
						-(yu(p,s).imag() * WC2R(CHq1R, s, t) + yu(p,s).real() * WC2I(CHq1I, s, t)) * yudag(t,r).imag())
						- 2. * (
						(NC * WC6R(CuuR, p, r, s, t) + NC * WC6R(CuuR, s, t, p, r)
						+ WC6R(CuuR, p, t, s, r) + WC6R(CuuR, s, r, p, t)) * gu(t,s).real()
						-(NC * WC6I(CuuI, p, r, s, t) + NC * WC6I(CuuI, s, t, p, r)
						+ WC6I(CuuI, p, t, s, r) + WC6I(CuuI, s, r, p, t)) * gu(t,s).imag())
						+ 2. * (WC7R(CeuR, s, t, p, r) * ge(t,s).real() - WC7I(CeuI, s, t, p, r) * ge(t,s).imag())
						- 2. * (WC7R(CluR, s, t, p, r) * yedye(t,s).real() - WC7I(CluI, s, t, p, r) * yedye(t,s).imag())
						+ 2. * NC * (WC7R(Cud1R, p, r, s, t) * gd(t,s).real() - WC7I(Cud1I, p, r, s, t) * gd(t,s).imag())
						- 2. * NC * (WC7R(Cqu1R, s, t, p, r)*(yddyd(t,s).real() - yudyu(t,s).real())
						- WC7I(Cqu1I, s, t, p, r)*(yddyd(t,s).imag() - yudyu(t,s).imag()));
				}
			}

			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2R_indices[i][0];
			s = WC2R_indices[i][1];
			//RGE 3
			f[c] += 0.5 * xiB * g12 * delta[r][s] * Yu
				+ FOUR_THIRDS * g12 * Yh2 * WC2R(CHuR, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += FOUR_THIRDS * g12 * YhYe * WC7R(CeuR, w, w, r, s)
					+ EIGHT_THIRDS * g12 * YhYl * WC7R(CluR, w, w, r, s)
					+ EIGHT_THIRDS * g12 * NC * YhYq * WC7R(Cqu1R, w, w, r, s)
					+ FOUR_THIRDS * g12 * NC * YhYd * WC7R(Cud1R, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYu * (
					NC * WC6R(CuuR, r, s, w, w) + WC6R(CuuR, r, w, w, s)
					+ WC6R(CuuR, w, s, r, w) + NC * WC6R(CuuR, w, w, r, s));

			}
			f[c] *= loop_factor;
			c++;
		}
		//CHuI
		for (i = 0; i < DWC2I; i++) {
			p = WC2I_indices[i][0];
			r = WC2I_indices[i][1];
			//RGE 2
			f[c] = -gu(p,r).imag()*(CHBOX + CHD) + 2. * gammaH * WC2I(CHuI, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 4. * (gu(p,t).imag() * WC2R(CHuR, t, r) + gu(p,t).real() * WC2I(CHuI, t, r))
					+ 4. * (WC2R(CHuR, p, t) * gu(t,r).imag() + WC2I(CHuI, p, t) * gu(t,r).real())
					+(yuydd(p,t).imag() * WC1(CHudR, r, t) - yuydd(p,t).real() * WC1(CHudI, r, t))
					+(WC1(CHudR, p, t) * ydyud(t,r).imag() + WC1(CHudI, p, t) * ydyud(t,r).real())
					;
				for (s = 0; s < NG; s++) {
					f[c] += -2. * ((yu(p,s).real() * WC2R(CHq1R, s, t) - yu(p,s).imag() * WC2I(CHq1I, s, t)) * yudag(t,r).imag()
						+(yu(p,s).imag() * WC2R(CHq1R, s, t) + yu(p,s).real() * WC2I(CHq1I, s, t)) * yudag(t,r).real())
						- 2. * (
						(NC * WC6R(CuuR, p, r, s, t) + NC * WC6R(CuuR, s, t, p, r)
						+ WC6R(CuuR, p, t, s, r) + WC6R(CuuR, s, r, p, t)) * gu(t,s).imag()
						+(NC * WC6I(CuuI, p, r, s, t) + NC * WC6I(CuuI, s, t, p, r)
						+ WC6I(CuuI, p, t, s, r) + WC6I(CuuI, s, r, p, t)) * gu(t,s).real())
						+ 2. * (WC7R(CeuR, s, t, p, r) * ge(t,s).imag() + WC7I(CeuI, s, t, p, r) * ge(t,s).real())
						- 2. * (WC7R(CluR, s, t, p, r) * yedye(t,s).imag() + WC7I(CluI, s, t, p, r) * yedye(t,s).real())
						+ 2. * NC * (WC7R(Cud1R, p, r, s, t) * gd(t,s).imag() + WC7I(Cud1I, p, r, s, t) * gd(t,s).real())
						- 2. * NC * (WC7R(Cqu1R, s, t, p, r)*(yddyd(t,s).imag() - yudyu(t,s).imag())
						+ WC7I(Cqu1I, s, t, p, r)*(yddyd(t,s).real() - yudyu(t,s).real()));
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2I_indices[i][0];
			s = WC2I_indices[i][1];
			//RGE 3
			f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHuI, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += FOUR_THIRDS * g12 * YhYe * WC7I(CeuI, w, w, r, s)
					+ EIGHT_THIRDS * g12 * YhYl * WC7I(CluI, w, w, r, s)
					+ EIGHT_THIRDS * g12 * NC * YhYq * WC7I(Cqu1I, w, w, r, s)
					+ FOUR_THIRDS * g12 * NC * YhYd * WC7I(Cud1I, r, s, w, w)
					+ FOUR_THIRDS * g12 * YhYu * (
					NC * WC6I(CuuI, r, s, w, w) + WC6I(CuuI, r, w, w, s)
					+ WC6I(CuuI, w, s, r, w) + NC * WC6I(CuuI, w, w, r, s));

			}
			f[c] *= loop_factor;
			c++;
		}




		//CHdR
		for (i = 0; i < DWC2R; i++) {
			p = WC2R_indices[i][0];
			r = WC2R_indices[i][1];
			//RGE 2
			f[c] = gd(p,r).real()*(CHBOX + CHD) + 2. * gammaH * WC2R(CHdR, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 4. * (gd(p,t).real() * WC2R(CHdR, t, r) - gd(p,t).imag() * WC2I(CHdI, t, r))
					+ 4. * (WC2R(CHdR, p, t) * gd(t,r).real() - WC2I(CHdI, p, t) * gd(t,r).imag())
					-(ydyud(p,t).real() * WC1(CHudR, t, r) - ydyud(p,t).imag() * WC1(CHudI, t, r))
					-(WC1(CHudR, t, p) * yuydd(t,r).real() + WC1(CHudI, t, p) * yuydd(t,r).imag());
				for (s = 0; s < NG; s++) {
					f[c] += -2. * ((yd(p,s).real() * WC2R(CHq1R, s, t) - yd(p,s).imag() * WC2I(CHq1I, s, t)) * yddag(t,r).real()
						-(yd(p,s).real() * WC2I(CHq1I, s, t) + yd(p,s).imag() * WC2R(CHq1R, s, t)) * yddag(t,r).imag())
						+ 2. * (
						(NC * WC6R(CddR, p, r, s, t) + NC * WC6R(CddR, s, t, p, r)
						+ WC6R(CddR, p, t, s, r) + WC6R(CddR, s, r, p, t)) * gd(t,s).real()
						-(NC * WC6I(CddI, p, r, s, t) + NC * WC6I(CddI, s, t, p, r)
						+ WC6I(CddI, p, t, s, r) + WC6I(CddI, s, r, p, t)) * gd(t,s).imag()
						)
						+ 2. * (WC7R(CedR, s, t, p, r) * ge(t,s).real() - WC7I(CedI, s, t, p, r) * ge(t,s).imag())
						- 2. * (WC7R(CldR, s, t, p, r) * yedye(t,s).real() - WC7I(CldI, s, t, p, r) * yedye(t,s).imag())
						- 2. * NC * (WC7R(Cud1R, s, t, p, r) * gu(t,s).real() - WC7I(Cud1I, s, t, p, r) * gu(t,s).imag())
						- 2. * NC * (WC7R(Cqd1R, s, t, p, r)*(yddyd(t,s).real() - yudyu(t,s).real())
						- WC7I(Cqd1I, s, t, p, r)*(yddyd(t,s).imag() - yudyu(t,s).imag()))
						;
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2R_indices[i][0];
			s = WC2R_indices[i][1];
			//RGE 3
			f[c] += 0.5 * xiB * g12 * delta[r][s] * Yd
				+ FOUR_THIRDS * g12 * Yh2 * WC2R(CHdR, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += FOUR_THIRDS * g12 * YhYd * (
					NC * WC6R(CddR, r, s, w, w) + WC6R(CddR, r, w, w, s)
					+ WC6R(CddR, w, s, r, w) + NC * WC6R(CddR, w, w, r, s))
					+ FOUR_THIRDS * g12 * YhYe * WC7R(CedR, w, w, r, s)
					+ EIGHT_THIRDS * g12 * YhYl * WC7R(CldR, w, w, r, s)
					+ EIGHT_THIRDS * g12 * NC * YhYq * WC7R(Cqd1R, w, w, r, s)
					+ FOUR_THIRDS * g12 * NC * YhYu * WC7R(Cud1R, w, w, r, s);
			}
			f[c] *= loop_factor;
			c++;
		}

		//CHdI  
		for (i = 0; i < DWC2I; i++) {
			p = WC2I_indices[i][0];
			r = WC2I_indices[i][1];
			//RGE 2
			f[c] = gd(p,r).imag()*(CHBOX + CHD) + 2. * gammaH * WC2I(CHdI, p, r);
			for (t = 0; t < NG; t++) {
				f[c] += 4. * (gd(p,t).imag() * WC2R(CHdR, t, r) + gd(p,t).real() * WC2I(CHdI, t, r))
					+ 4. * (WC2R(CHdR, p, t) * gd(t,r).imag() + WC2I(CHdI, p, t) * gd(t,r).real())
					-(ydyud(p,t).imag() * WC1(CHudR, t, r) + ydyud(p,t).real() * WC1(CHudI, t, r))
					-(WC1(CHudR, t, p) * yuydd(t,r).imag() - WC1(CHudI, t, p) * yuydd(t,r).real());
				for (s = 0; s < NG; s++) {
					f[c] += -2. * ((yd(p,s).real() * WC2R(CHq1R, s, t) - yd(p,s).imag() * WC2I(CHq1I, s, t)) * yddag(t,r).imag()
						+(yd(p,s).real() * WC2I(CHq1I, s, t) + yd(p,s).imag() * WC2R(CHq1R, s, t)) * yddag(t,r).real())
						+ 2. * (
						(NC * WC6R(CddR, p, r, s, t) + NC * WC6R(CddR, s, t, p, r)
						+ WC6R(CddR, p, t, s, r) + WC6R(CddR, s, r, p, t)) * gd(t,s).imag()
						+(NC * WC6I(CddI, p, r, s, t) + NC * WC6I(CddI, s, t, p, r)
						+ WC6I(CddI, p, t, s, r) + WC6I(CddI, s, r, p, t)) * gd(t,s).real()
						)
						+ 2. * (WC7R(CedR, s, t, p, r) * ge(t,s).imag() + WC7I(CedI, s, t, p, r) * ge(t,s).real())
						- 2. * (WC7R(CldR, s, t, p, r) * yedye(t,s).imag() + WC7I(CldI, s, t, p, r) * yedye(t,s).real())
						- 2. * NC * (WC7R(Cud1R, s, t, p, r) * gu(t,s).imag() + WC7I(Cud1I, s, t, p, r) * gu(t,s).real())
						- 2. * NC * (WC7R(Cqd1R, s, t, p, r)*(yddyd(t,s).imag() - yudyu(t,s).imag())
						+ WC7I(Cqd1I, s, t, p, r)*(yddyd(t,s).real() - yudyu(t,s).real()))
						;
				}
			}
			//Necessary since RGE 2 and RGE 3 use different indices.
			r = WC2I_indices[i][0];
			s = WC2I_indices[i][1];
			//RGE 3
			f[c] += FOUR_THIRDS * g12 * Yh2 * WC2I(CHdI, r, s);
			for (w = 0; w < NG; w++) {
				f[c] += FOUR_THIRDS * g12 * YhYd * (
					NC * WC6I(CddI, r, s, w, w) + WC6I(CddI, r, w, w, s)
					+ WC6I(CddI, w, s, r, w) + NC * WC6I(CddI, w, w, r, s))
					+ FOUR_THIRDS * g12 * YhYe * WC7I(CedI, w, w, r, s)
					+ EIGHT_THIRDS * g12 * YhYl * WC7I(CldI, w, w, r, s)
					+ EIGHT_THIRDS * g12 * NC * YhYq * WC7I(Cqd1I, w, w, r, s)
					+ FOUR_THIRDS * g12 * NC * YhYu * WC7I(Cud1I, w, w, r, s);
			}
			f[c] *= loop_factor;
			c++;
		}

		for (p = 0; p < NG; p++) {
			for (r = 0; r < NG; r++) {
				//CHudR
				f[c] = yuydd(p,r).real()*(2. * CHBOX - CHD) + 2. * gammaH * WC1(CHudR, p, r);
				for (t = 0; t < NG; t++) {
					f[c] += -2. * (yuydd(p,t).real() * WC2R(CHdR, t, r) - yuydd(p,t).imag() * WC2I(CHdI, t, r))
						+ 2. * (WC2R(CHuR, p, t) * yuydd(r,t).real() - WC2I(CHuI, p, t) * yuydd(r,t).imag())
						+ 3. * (gu(p,t).real() * WC1(CHudR, t, r) - gu(p,t).imag() * WC1(CHudI, t, r))
						+ 3. * (WC1(CHudR, p, t) * gd(t,r).real() - WC1(CHudI, p, t) * gd(t,r).imag())
						;
					for (s = 0; s < NG; s++) {
						f[c] += 4. * (
							(WC7R(Cud1R, p, t, s, r) + cF3 * WC7R(Cud8R, p, t, s, r)) * yuydd(t,s).real()
							-(WC7I(Cud1I, p, t, s, r) + cF3 * WC7I(Cud8I, p, t, s, r)) * yuydd(t,s).imag()
							);
					}
				}
				//RGE 3 
				f[c] += -3. * g12 * (Yu - Yd)*(Yu - Yd) * WC1(CHudR, p, r);

				f[c] *= loop_factor;

				//CHudI
				f[c + DF] = yuydd(p,r).imag()*(2. * CHBOX - CHD) + 2. * gammaH * WC1(CHudI, p, r);
				for (t = 0; t < NG; t++) {
					f[c + DF] += -2. * (yuydd(p,t).imag() * WC2R(CHdR, t, r) + yuydd(p,t).real() * WC2I(CHdI, t, r))
						+ 2. * (WC2R(CHuR, p, t) * yuydd(r,t).imag() + WC2I(CHuI, p, t) * yuydd(r,t).real())
						+ 3. * (gu(p,t).imag() * WC1(CHudR, t, r) + gu(p,t).real() * WC1(CHudI, t, r))
						+ 3. * (WC1(CHudR, p, t) * gd(t,r).imag() + WC1(CHudI, p, t) * gd(t,r).real())
						;
					for (s = 0; s < NG; s++) {
						f[c + DF] += 4. * (
							(WC7R(Cud1R, p, t, s, r) + cF3 * WC7R(Cud8R, p, t, s, r)) * yuydd(t,s).imag()
							+(WC7I(Cud1I, p, t, s, r) + cF3 * WC7I(Cud8I, p, t, s, r)) * yuydd(t,s).real()
							);
					}
				}
				//RGE 3 
				f[c + DF] += -3. * g12 * (Yu - Yd)*(Yu - Yd) * WC1(CHudI, p, r);


				f[c + DF] *= loop_factor;
				c++;
			}
		}
		c += DF;

	}
	//----------------------RGE SMEFT class 8_LLLL----------------------
	//CllR
	for (d = 0; d < DWC6R; d++) {
		p = WC6R_indices[d][0];
		r = WC6R_indices[d][1];
		s = WC6R_indices[d][2];
		t = WC6R_indices[d][3];

		f[c] = -0.5 * (yedye(p,r).real() * WC2R(CHl1R, s, t) - yedye(p,r).imag() * WC2I(CHl1I, s, t))
			- 0.5 * (yedye(s,t).real() * WC2R(CHl1R, p, r) - yedye(s,t).imag() * WC2I(CHl1I, p, r))
			+ 0.5 * (yedye(p,r).real() * WC2R(CHl3R, s, t) - yedye(p,r).imag() * WC2I(CHl3I, s, t))
			+ 0.5 * (yedye(s,t).real() * WC2R(CHl3R, p, r) - yedye(s,t).imag() * WC2I(CHl3I, p, r))
			-(yedye(s,r).real() * WC2R(CHl3R, p, t) - yedye(s,r).imag() * WC2I(CHl3I, p, t))
			-(yedye(p,t).real() * WC2R(CHl3R, s, r) - yedye(p,t).imag() * WC2I(CHl3I, s, r))//RGE 2
			+ TWO_THIRDS * g12 * YhYl * (WC2R(CHl1R, s, t) * delta[p][r] + WC2R(CHl1R, p, r) * delta[s][t])
			- ONE_SIXTH * g22 * (WC2R(CHl3R, s, t) * delta[p][r] + WC2R(CHl3R, p, r) * delta[s][t])
			+ ONE_THIRD * g22 * (WC2R(CHl3R, s, r) * delta[p][t] + WC2R(CHl3R, p, t) * delta[r][s])
			+ 6. * g22 * WC6R(CllR, p, t, s, r) - 3. * (g22 - 4. * Yl2 * g12) * WC6R(CllR, p, r, s, t)
			//RGE 3
			;
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += (gl(p,v).real() * WC6R(CllR, v, r, s, t) - gl(p,v).imag() * WC6I(CllI, v, r, s, t))
				+(gl(s,v).real() * WC6R(CllR, p, r, v, t) - gl(s,v).imag() * WC6I(CllI, p, r, v, t))
				+(WC6R(CllR, p, v, s, t) * gl(v,r).real() - WC6I(CllI, p, v, s, t) * gl(v,r).imag())
				+(WC6R(CllR, p, r, s, v) * gl(v,t).real() - WC6I(CllI, p, r, s, v) * gl(v,t).imag())
				//RGE 2
				+ FOUR_THIRDS * g12 * Yl2 * (
				WC6R(CllR, p, r, w, w) * delta[s][t] + WC6R(CllR, s, t, w, w) * delta[p][r]
				+ WC6R(CllR, w, w, s, t) * delta[p][r] + WC6R(CllR, w, w, p, r) * delta[s][t])
				+ TWO_THIRDS * g12 * Yl2 * (
				WC6R(CllR, p, w, w, r) * delta[s][t] + WC6R(CllR, s, w, w, t) * delta[p][r]
				+ WC6R(CllR, w, r, p, w) * delta[s][t] + WC6R(CllR, w, t, s, w) * delta[p][r])
				- ONE_SIXTH * g22 * (
				WC6R(CllR, p, w, w, r) * delta[s][t] + WC6R(CllR, s, w, w, t) * delta[p][r]
				+ WC6R(CllR, p, w, w, r) * delta[s][t] + WC6R(CllR, w, t, s, w) * delta[p][r])
				+ ONE_THIRD * g22 * (
				WC6R(CllR, s, w, w, r) * delta[p][t] + WC6R(CllR, p, w, w, t) * delta[r][s]
				+ WC6R(CllR, w, r, s, w) * delta[p][t] + WC6R(CllR, w, t, p, w) * delta[r][s])
				+ FOUR_THIRDS * g12 * NC * YlYq * (
				WC7R(Clq1R, p, r, w, w) * delta[s][t] + WC7R(Clq1R, s, t, w, w) * delta[p][r])
				+ ONE_THIRD * g22 * NC *
				(-WC7R(Clq3R, p, r, w, w) * delta[s][t] - WC7R(Clq3R, s, t, w, w) * delta[p][r]
				+ 2. * WC7R(Clq3R, s, r, w, w) * delta[p][t] + 2. * WC7R(Clq3R, p, t, w, w) * delta[r][s])
				+ TWO_THIRDS * g12 * NC * YuYl * (
				WC7R(CluR, p, r, w, w) * delta[s][t] + WC7R(CluR, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YdYl * (
				WC7R(CldR, p, r, w, w) * delta[s][t] + WC7R(CldR, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * YeYl * (
				WC7R(CleR, p, r, w, w) * delta[s][t] + WC7R(CleR, s, t, w, w) * delta[p][r]);
			//RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -0.5 * (
					+(yedag(s,v).real() * ye(w,t).real() - yedag(s,v).imag() * ye(w,t).imag()) * WC7R(CleR, p, r, v, w)
					-(yedag(s,v).real() * ye(w,t).imag() + yedag(s,v).imag() * ye(w,t).real()) * WC7I(CleI, p, r, v, w))
					- 0.5 * ((yedag(p,v).real() * ye(w,r).real() - yedag(p,v).imag() * ye(w,r).imag()) * WC7R(CleR, s, t, v, w)
					-(yedag(p,v).real() * ye(w,r).imag() + yedag(p,v).imag() * ye(w,r).real()) * WC7I(CleI, s, t, v, w))
					//RGE 2
					;
			}
		}

		f[c] *= loop_factor;
		c++;
	}
	//CllI
	for (d = 0; d < DWC6I; d++) {
		p = WC6I_indices[d][0];
		r = WC6I_indices[d][1];
		s = WC6I_indices[d][2];
		t = WC6I_indices[d][3];

		f[c] = -0.5 * (yedye(p,r).imag() * WC2R(CHl1R, s, t) + yedye(p,r).real() * WC2I(CHl1I, s, t))
			- 0.5 * (yedye(s,t).imag() * WC2R(CHl1R, p, r) + yedye(s,t).real() * WC2I(CHl1I, p, r))
			+ 0.5 * (yedye(p,r).imag() * WC2R(CHl3R, s, t) + yedye(p,r).real() * WC2I(CHl3I, s, t))
			+ 0.5 * (yedye(s,t).imag() * WC2R(CHl3R, p, r) + yedye(s,t).real() * WC2I(CHl3I, p, r))
			-(yedye(s,r).imag() * WC2R(CHl3R, p, t) + yedye(s,r).real() * WC2I(CHl3I, p, t))
			-(yedye(p,t).imag() * WC2R(CHl3R, s, r) + yedye(p,t).real() * WC2I(CHl3I, s, r))//RGE 2
			+ TWO_THIRDS * g12 * YhYl * (WC2I(CHl1I, s, t) * delta[p][r] + WC2I(CHl1I, p, r) * delta[s][t])
			- ONE_SIXTH * g22 * (WC2I(CHl3I, s, t) * delta[p][r] + WC2I(CHl3I, p, r) * delta[s][t])
			+ ONE_THIRD * g22 * (WC2I(CHl3I, s, r) * delta[p][t] + WC2I(CHl3I, p, t) * delta[r][s])
			+ 6. * g22 * WC6I(CllI, p, t, s, r) - 3. * (g22 - 4. * Yl2 * g12) * WC6I(CllI, p, r, s, t)
			//RGE 3
			;
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += (gl(p,v).imag() * WC6R(CllR, v, r, s, t) + gl(p,v).real() * WC6I(CllI, v, r, s, t))
				+(gl(s,v).imag() * WC6R(CllR, p, r, v, t) + gl(s,v).real() * WC6I(CllI, p, r, v, t))
				+(WC6R(CllR, p, v, s, t) * gl(v,r).imag() + WC6I(CllI, p, v, s, t) * gl(v,r).real())
				+(WC6R(CllR, p, r, s, v) * gl(v,t).imag() + WC6I(CllI, p, r, s, v) * gl(v,t).real())
				//RGE 2
				+ FOUR_THIRDS * g12 * Yl2 * (
				WC6I(CllI, p, r, w, w) * delta[s][t] + WC6I(CllI, s, t, w, w) * delta[p][r]
				+ WC6I(CllI, w, w, s, t) * delta[p][r] + WC6I(CllI, w, w, p, r) * delta[s][t])
  				+ TWO_THIRDS * g12 * Yl2 * (
				WC6I(CllI, p, w, w, r) * delta[s][t] + WC6I(CllI, s, w, w, t) * delta[p][r]
				+ WC6I(CllI, w, r, p, w) * delta[s][t] + WC6I(CllI, w, t, s, w) * delta[p][r])
				- ONE_SIXTH * g22 * (
				WC6I(CllI, p, w, w, r) * delta[s][t] + WC6I(CllI, s, w, w, t) * delta[p][r]
				+ WC6I(CllI, p, w, w, r) * delta[s][t] + WC6I(CllI, w, t, s, w) * delta[p][r])
				+ ONE_THIRD * g22 * (
				WC6I(CllI, s, w, w, r) * delta[p][t] + WC6I(CllI, p, w, w, t) * delta[r][s]
				+ WC6I(CllI, w, r, s, w) * delta[p][t] + WC6I(CllI, w, t, p, w) * delta[r][s])
				+ FOUR_THIRDS * g12 * NC * YlYq * (
				WC7I(Clq1I, p, r, w, w) * delta[s][t] + WC7I(Clq1I, s, t, w, w) * delta[p][r])
				+ ONE_THIRD * g22 * NC *
				(-WC7I(Clq3I, p, r, w, w) * delta[s][t] - WC7I(Clq3I, s, t, w, w) * delta[p][r]
				+ 2. * WC7I(Clq3I, s, r, w, w) * delta[p][t] + 2. * WC7I(Clq3I, p, t, w, w) * delta[r][s])
				+ TWO_THIRDS * g12 * NC * YuYl * (
				WC7I(CluI, p, r, w, w) * delta[s][t] + WC7I(CluI, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YdYl * (
				WC7I(CldI, p, r, w, w) * delta[s][t] + WC7I(CldI, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * YeYl * (
				WC7I(CleI, p, r, w, w) * delta[s][t] + WC7I(CleI, s, t, w, w) * delta[p][r]);
			//RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -0.5 * ((yedag(s,v).real() * ye(w,t).real() - yedag(s,v).imag() * ye(w,t).imag()) * WC7I(CleI, p, r, v, w)
					+(yedag(s,v).real() * ye(w,t).imag() + yedag(s,v).imag() * ye(w,t).real()) * WC7R(CleR, p, r, v, w))
					- 0.5 * ((yedag(p,v).real() * ye(w,r).real() - yedag(p,v).imag() * ye(w,r).imag()) * WC7I(CleI, s, t, v, w)
					+(yedag(p,v).real() * ye(w,r).imag() + yedag(p,v).imag() * ye(w,r).real()) * WC7R(CleR, s, t, v, w))
					//RGE 2
					;

			}
		}

		f[c] *= loop_factor;
		c++;
	}
	//Cqq1R
	for (d = 0; d < DWC6R; d++) {
		p = WC6R_indices[d][0];
		r = WC6R_indices[d][1];
		s = WC6R_indices[d][2];
		t = WC6R_indices[d][3];
		f[c] = 0.5 * ((yudyu(p,r).real() - yddyd(p,r).real()) * WC2R(CHq1R, s, t)
			-(yudyu(p,r).imag() - yddyd(p,r).imag()) * WC2I(CHq1I, s, t))
			+ 0.5 * ((yudyu(s,t).real() - yddyd(s,t).real()) * WC2R(CHq1R, p, r)
			-(yudyu(s,t).imag() - yddyd(s,t).imag()) * WC2I(CHq1I, p, r))//RGE 2
			+ TWO_THIRDS * g12 * YhYq * (WC2R(CHq1R, s, t) * delta[p][r] + WC2R(CHq1R, p, r) * delta[s][t])
			+ 3. * g32 * WC6R(Cqq1R, p, t, s, r) + 9. * g32 * WC6R(Cqq3R, p, t, s, r)
			+ 9. * g22 * WC6R(Cqq3R, p, r, s, t)
			-(6. / NC)*(g32 - 2. * NC * Yq2 * g12) * WC6R(Cqq1R, p, r, s, t)//RGE 3
			;
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += (gq(p,v).real() * WC6R(Cqq1R, v, r, s, t) - gq(p,v).imag() * WC6I(Cqq1I, v, r, s, t))
				+(gq(s,v).real() * WC6R(Cqq1R, p, r, v, t) - gq(s,v).imag() * WC6I(Cqq1I, p, r, v, t))
				+(WC6R(Cqq1R, p, v, s, t) * gq(v,r).real() - WC6I(Cqq1I, p, v, s, t) * gq(v,r).imag())
				+(WC6R(Cqq1R, p, r, s, v) * gq(v,t).real() - WC6I(Cqq1I, p, r, s, v) * gq(v,t).imag())
				//RGE 2
				+ FOUR_THIRDS * g12 * YlYq * (WC7R(Clq1R, w, w, s, t) * delta[p][r]
				+ WC7R(Clq1R, w, w, p, r) * delta[s][t])
				+ FOUR_THIRDS * g12 * Yq2 * NC * (
				WC6R(Cqq1R, p, r, w, w) * delta[s][t] + WC6R(Cqq1R, s, t, w, w) * delta[p][r]
				+ WC6R(Cqq1R, w, w, s, t) * delta[p][r] + WC6R(Cqq1R, w, w, p, r) * delta[s][t])
				+ TWO_THIRDS * g12 * Yq2 * (
				WC6R(Cqq1R, p, w, w, r) * delta[s][t] + WC6R(Cqq1R, s, w, w, t) * delta[p][r]
				+ WC6R(Cqq1R, w, r, p, w) * delta[s][t] + WC6R(Cqq1R, w, t, s, w) * delta[p][r])
				+ ONE_SIXTH * g32 * (
				WC6R(Cqq1R, s, w, w, r) * delta[p][t] + WC6R(Cqq1R, p, w, w, t) * delta[r][s]
				+ WC6R(Cqq1R, w, r, s, w) * delta[p][t] + WC6R(Cqq1R, w, t, p, w) * delta[r][s])
				- g32 * (ONE_THIRD / NC)*(
				WC6R(Cqq1R, p, w, w, r) * delta[s][t] + WC6R(Cqq1R, s, w, w, t) * delta[p][r]
				+ WC6R(Cqq1R, w, r, p, w) * delta[s][t] + WC6R(Cqq1R, w, t, s, w) * delta[p][r])
				+ 2. * g12 * Yq2 * (
				WC6R(Cqq3R, p, w, w, r) * delta[s][t] + WC6R(Cqq3R, s, w, w, t) * delta[p][r]
				+ WC6R(Cqq3R, w, r, p, w) * delta[s][t] + WC6R(Cqq3R, w, t, s, w) * delta[p][r])
				+ 0.5 * g32 * (
				WC6R(Cqq3R, s, w, w, r) * delta[p][t] + WC6R(Cqq3R, p, w, w, t) * delta[r][s]
				+ WC6R(Cqq3R, w, r, s, w) * delta[p][t] + WC6R(Cqq3R, w, t, p, w) * delta[r][s])
				-(1. / NC) * g32 * (
				WC6R(Cqq3R, p, w, w, r) * delta[s][t] + WC6R(Cqq3R, s, w, w, t) * delta[p][r]
				+ WC6R(Cqq3R, w, r, p, w) * delta[s][t] + WC6R(Cqq3R, w, t, s, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YuYq *
				(WC7R(Cqu1R, p, r, w, w) * delta[s][t] + WC7R(Cqu1R, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YdYq *
				(WC7R(Cqd1R, p, r, w, w) * delta[s][t] + WC7R(Cqd1R, s, t, w, w) * delta[p][r])
				+(g32 / 12.)*(
				WC7R(Cqu8R, s, r, w, w) * delta[p][t] + WC7R(Cqu8R, p, t, w, w) * delta[r][s]
				-(2. / NC) * WC7R(Cqu8R, p, r, w, w) * delta[s][t]-(2. / NC) * WC7R(Cqu8R, s, t, w, w) * delta[p][r])
				+(g32 / 12.)*(
				WC7R(Cqd8R, s, r, w, w) * delta[p][t] + WC7R(Cqd8R, p, t, w, w) * delta[r][s]
				-(2. / NC) * WC7R(Cqd8R, p, r, w, w) * delta[s][t]-(2. / NC) * WC7R(Cqd8R, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * YeYq *
				(WC7R(CqeR, p, r, w, w) * delta[s][t] + WC7R(CqeR, s, t, w, w) * delta[p][r])
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += +(0.25 / NC)*(
					(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC7R(Cqu8R, s, t, v, w)
					-(yudag(p,v).real() * yu(w,r).imag() + yudag(p,v).imag() * yu(w,r).real()) * WC7I(Cqu8I, s, t, v, w)
					+(yudag(s,v).real() * yu(w,t).real() - yudag(s,v).imag() * yu(w,t).imag()) * WC7R(Cqu8R, p, r, v, w)
					-(yudag(s,v).real() * yu(w,t).imag() + yudag(s,v).imag() * yu(w,t).real()) * WC7I(Cqu8I, p, r, v, w))
					+(0.25 / NC)*(
					(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC7R(Cqd8R, s, t, v, w)
					-(yddag(p,v).real() * yd(w,r).imag() + yddag(p,v).imag() * yd(w,r).real()) * WC7I(Cqd8I, s, t, v, w)
					+(yddag(s,v).real() * yd(w,t).real() - yddag(s,v).imag() * yd(w,t).imag()) * WC7R(Cqd8R, p, r, v, w)
					-(yddag(s,v).real() * yd(w,t).imag() + yddag(s,v).imag() * yd(w,t).real()) * WC7I(Cqd8I, p, r, v, w))

					- 0.125 * (
					(yudag(p,v).real() * yu(w,t).real() - yudag(p,v).imag() * yu(w,t).imag()) * WC7R(Cqu8R, s, r, v, w)
					-(yudag(p,v).real() * yu(w,t).imag() + yudag(p,v).imag() * yu(w,t).real()) * WC7I(Cqu8I, s, r, v, w)
					+(yudag(s,v).real() * yu(w,r).real() - yudag(s,v).imag() * yu(w,r).imag()) * WC7R(Cqu8R, p, t, v, w)
					-(yudag(s,v).real() * yu(w,r).imag() + yudag(s,v).imag() * yu(w,r).real()) * WC7I(Cqu8I, p, t, v, w))
					- 0.125 * (
					(yddag(p,v).real() * yd(w,t).real() - yddag(p,v).imag() * yd(w,t).imag()) * WC7R(Cqd8R, s, r, v, w)
					-(yddag(p,v).real() * yd(w,t).imag() + yddag(p,v).imag() * yd(w,t).real()) * WC7I(Cqd8I, s, r, v, w)
					+(yddag(s,v).real() * yd(w,r).real() - yddag(s,v).imag() * yd(w,r).imag()) * WC7R(Cqd8R, p, t, v, w)
					-(yddag(s,v).real() * yd(w,r).imag() + yddag(s,v).imag() * yd(w,r).real()) * WC7I(Cqd8I, p, t, v, w))

					+(0.0625 / NC) *(
					(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd8R, p, v, s, w)
					-(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd8I, p, v, s, w)
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd8R, s, v, p, w)
					-(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd8I, s, v, p, w))
					+(0.0625 / NC) *(
					(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd8R, r, v, t, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd8I, r, v, t, w)
					+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd8R, t, v, r, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd8I, t, v, r, w))

					+ 0.0625 * (
					(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd8R, s, v, p, w)
					-(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd8I, s, v, p, w)
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd8R, p, v, s, w)
					-(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd8I, p, v, s, w))
					+ 0.0625 * (
					(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd8R, t, v, r, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd8I, t, v, r, w)
					+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd8R, r, v, t, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd8I, r, v, t, w))

					- 0.5 * ((yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC7R(Cqu1R, s, t, v, w)
					-(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC7I(Cqu1I, s, t, v, w))
					- 0.5 * ((yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC7R(Cqd1R, s, t, v, w)
					-(yddag(s,v).imag() * yd(w,r).real() + yddag(s,v).real() * yd(w,r).imag()) * WC7I(Cqd1I, s, t, v, w))
					- 0.5 * ((yudag(s,v).real() * yu(w,t).real() - yudag(s,v).imag() * yu(w,t).imag()) * WC7R(Cqu1R, p, r, v, w)
					-(yudag(s,v).imag() * yu(w,t).real() + yudag(s,v).real() * yu(w,t).imag()) * WC7I(Cqu1I, p, r, v, w))
					- 0.5 * ((yddag(s,v).real() * yd(w,t).real() - yddag(s,v).imag() * yd(w,t).imag()) * WC7R(Cqd1R, p, r, v, w)
					-(yddag(s,v).imag() * yd(w,t).real() + yddag(s,v).real() * yd(w,t).imag()) * WC7I(Cqd1I, p, r, v, w))

					- 0.125 * ((yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd1R, p, v, s, w)
					-(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd1I, p, v, s, w))
					- 0.125 * ((yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd1R, r, v, t, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd1I, r, v, t, w))
					- 0.125 * ((yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd1R, s, v, p, w)
					-(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd1I, s, v, p, w))
					- 0.125 * ((yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd1R, t, v, r, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd1I, t, v, r, w))
					//RGE 2
					;
			}
		}

		f[c] *= loop_factor;
		c++;
	}


	//Cqq1I
	for (d = 0; d < DWC6I; d++) {
		p = WC6I_indices[d][0];
		r = WC6I_indices[d][1];
		s = WC6I_indices[d][2];
		t = WC6I_indices[d][3];
		f[c] = 0.5 * ((yudyu(p,r).imag() - yddyd(p,r).imag()) * WC2R(CHq1R, s, t)
			+(yudyu(p,r).real() - yddyd(p,r).real()) * WC2I(CHq1I, s, t))
			+ 0.5 * ((yudyu(s,t).imag() - yddyd(s,t).imag()) * WC2R(CHq1R, p, r)
			+(yudyu(s,t).real() - yddyd(s,t).real()) * WC2I(CHq1I, p, r))//RGE 2
			+ TWO_THIRDS * g12 * YhYq * (WC2I(CHq1I, s, t) * delta[p][r] + WC2I(CHq1I, p, r) * delta[s][t])
			+ 3. * g32 * WC6I(Cqq1I, p, t, s, r) + 9. * g32 * WC6I(Cqq3I, p, t, s, r)
			+ 9. * g22 * WC6I(Cqq3I, p, r, s, t)
			-(6. / NC)*(g32 - 2. * NC * Yq2 * g12) * WC6I(Cqq1I, p, r, s, t)//RGE 3
			;
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += (gq(p,v).imag() * WC6R(Cqq1R, v, r, s, t) + gq(p,v).real() * WC6I(Cqq1I, v, r, s, t))
				+(gq(s,v).imag() * WC6R(Cqq1R, p, r, v, t) + gq(s,v).real() * WC6I(Cqq1I, p, r, v, t))
				+(WC6R(Cqq1R, p, v, s, t) * gq(v,r).imag() + WC6I(Cqq1I, p, v, s, t) * gq(v,r).real())
				+(WC6R(Cqq1R, p, r, s, v) * gq(v,t).imag() + WC6I(Cqq1I, p, r, s, v) * gq(v,t).real()) //RGE 2
				//RGE 2
				+ FOUR_THIRDS * g12 * YlYq * (WC7I(Clq1I, w, w, s, t) * delta[p][r]
				+ WC7I(Clq1I, w, w, p, r) * delta[s][t])
				+ FOUR_THIRDS * g12 * Yq2 * NC * (
				WC6I(Cqq1I, p, r, w, w) * delta[s][t] + WC6I(Cqq1I, s, t, w, w) * delta[p][r]
				+ WC6I(Cqq1I, w, w, s, t) * delta[p][r] + WC6I(Cqq1I, w, w, p, r) * delta[s][t])
				+ TWO_THIRDS * g12 * Yq2 * (
				WC6I(Cqq1I, p, w, w, r) * delta[s][t] + WC6I(Cqq1I, s, w, w, t) * delta[p][r]
				+ WC6I(Cqq1I, w, r, p, w) * delta[s][t] + WC6I(Cqq1I, w, t, s, w) * delta[p][r])
				+ ONE_SIXTH * g32 * (
				WC6I(Cqq1I, s, w, w, r) * delta[p][t] + WC6I(Cqq1I, p, w, w, t) * delta[r][s]
				+ WC6I(Cqq1I, w, r, s, w) * delta[p][t] + WC6I(Cqq1I, w, t, p, w) * delta[r][s])
				- g32 * (ONE_THIRD / NC)*(
				WC6I(Cqq1I, p, w, w, r) * delta[s][t] + WC6I(Cqq1I, s, w, w, t) * delta[p][r]
				+ WC6I(Cqq1I, w, r, p, w) * delta[s][t] + WC6I(Cqq1I, w, t, s, w) * delta[p][r])
				+ 2. * g12 * Yq2 * (
				WC6I(Cqq3I, p, w, w, r) * delta[s][t] + WC6I(Cqq3I, s, w, w, t) * delta[p][r]
				+ WC6I(Cqq3I, w, r, p, w) * delta[s][t] + WC6I(Cqq3I, w, t, s, w) * delta[p][r])
				+ 0.5 * g32 * (
				WC6I(Cqq3I, s, w, w, r) * delta[p][t] + WC6I(Cqq3I, p, w, w, t) * delta[r][s]
				+ WC6I(Cqq3I, w, r, s, w) * delta[p][t] + WC6I(Cqq3I, w, t, p, w) * delta[r][s])
				-(1. / NC) * g32 * (
				WC6I(Cqq3I, p, w, w, r) * delta[s][t] + WC6I(Cqq3I, s, w, w, t) * delta[p][r]
				+ WC6I(Cqq3I, w, r, p, w) * delta[s][t] + WC6I(Cqq3I, w, t, s, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YuYq *
				(WC7I(Cqu1I, p, r, w, w) * delta[s][t] + WC7I(Cqu1I, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YdYq *
				(WC7I(Cqd1I, p, r, w, w) * delta[s][t] + WC7I(Cqd1I, s, t, w, w) * delta[p][r])
				+(g32 / 12.)*(
				WC7I(Cqu8I, s, r, w, w) * delta[p][t] + WC7I(Cqu8I, p, t, w, w) * delta[r][s]
				-(2. / NC) * WC7I(Cqu8I, p, r, w, w) * delta[s][t]-(2. / NC) * WC7I(Cqu8I, s, t, w, w) * delta[p][r])
				+(g32 / 12.)*(
				WC7I(Cqd8I, s, r, w, w) * delta[p][t] + WC7I(Cqd8I, p, t, w, w) * delta[r][s]
				-(2. / NC) * WC7I(Cqd8I, p, r, w, w) * delta[s][t]-(2. / NC) * WC7I(Cqd8I, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * YeYq *
				(WC7I(CqeI, p, r, w, w) * delta[s][t] + WC7I(CqeI, s, t, w, w) * delta[p][r])
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += +(0.25 / NC)*(
					(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC7I(Cqu8I, s, t, v, w)
					+(yudag(p,v).real() * yu(w,r).imag() + yudag(p,v).imag() * yu(w,r).real()) * WC7R(Cqu8R, s, t, v, w)
					+(yudag(s,v).real() * yu(w,t).real() - yudag(s,v).imag() * yu(w,t).imag()) * WC7I(Cqu8I, p, r, v, w)
					+(yudag(s,v).real() * yu(w,t).imag() + yudag(s,v).imag() * yu(w,t).real()) * WC7R(Cqu8R, p, r, v, w))
					+(0.25 / NC)*(
					(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC7I(Cqd8I, s, t, v, w)
					+(yddag(p,v).real() * yd(w,r).imag() + yddag(p,v).imag() * yd(w,r).real()) * WC7R(Cqd8R, s, t, v, w)
					+(yddag(s,v).real() * yd(w,t).real() - yddag(s,v).imag() * yd(w,t).imag()) * WC7I(Cqd8I, p, r, v, w)
					+(yddag(s,v).real() * yd(w,t).imag() + yddag(s,v).imag() * yd(w,t).real()) * WC7R(Cqd8R, p, r, v, w))

					- 0.125 * (
					(yudag(p,v).real() * yu(w,t).real() - yudag(p,v).imag() * yu(w,t).imag()) * WC7I(Cqu8I, s, r, v, w)
					+(yudag(p,v).real() * yu(w,t).imag() + yudag(p,v).imag() * yu(w,t).real()) * WC7R(Cqu8R, s, r, v, w)
					+(yudag(s,v).real() * yu(w,r).real() - yudag(s,v).imag() * yu(w,r).imag()) * WC7I(Cqu8I, p, t, v, w)
					+(yudag(s,v).real() * yu(w,r).imag() + yudag(s,v).imag() * yu(w,r).real()) * WC7R(Cqu8R, p, t, v, w))
					- 0.125 * (
					(yddag(p,v).real() * yd(w,t).real() - yddag(p,v).imag() * yd(w,t).imag()) * WC7I(Cqd8I, s, r, v, w)
					+(yddag(p,v).real() * yd(w,t).imag() + yddag(p,v).imag() * yd(w,t).real()) * WC7R(Cqd8R, s, r, v, w)
					+(yddag(s,v).real() * yd(w,r).real() - yddag(s,v).imag() * yd(w,r).imag()) * WC7I(Cqd8I, p, t, v, w)
					+(yddag(s,v).real() * yd(w,r).imag() + yddag(s,v).imag() * yd(w,r).real()) * WC7R(Cqd8R, p, t, v, w))

					+(0.0625 / NC) *(
					(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd8I, p, v, s, w)
					+(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd8R, p, v, s, w)
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd8I, s, v, p, w)
					+(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd8R, s, v, p, w))
					+(0.0625 / NC) *(
					-(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd8I, r, v, t, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd8R, r, v, t, w)
					-(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd8I, t, v, r, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd8R, t, v, r, w))

					+ 0.0625 * (
					(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd8I, s, v, p, w)
					+(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd8R, s, v, p, w)
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd8I, p, v, s, w)
					+(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd8R, p, v, s, w))
					+ 0.0625 * (
					-(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd8I, t, v, r, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd8R, t, v, r, w)
					-(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd8I, r, v, t, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd8R, r, v, t, w))

					- 0.5 * ((yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC7I(Cqu1I, s, t, v, w)
					+(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC7R(Cqu1R, s, t, v, w))
					- 0.5 * ((yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC7I(Cqd1I, s, t, v, w)
					+(yddag(s,v).imag() * yd(w,r).real() + yddag(s,v).real() * yd(w,r).imag()) * WC7R(Cqd1R, s, t, v, w))
					- 0.5 * ((yudag(s,v).real() * yu(w,t).real() - yudag(s,v).imag() * yu(w,t).imag()) * WC7I(Cqu1I, p, r, v, w)
					+(yudag(s,v).imag() * yu(w,t).real() + yudag(s,v).real() * yu(w,t).imag()) * WC7R(Cqu1R, p, r, v, w))
					- 0.5 * ((yddag(s,v).real() * yd(w,t).real() - yddag(s,v).imag() * yd(w,t).imag()) * WC7I(Cqd1I, p, r, v, w)
					+(yddag(s,v).imag() * yd(w,t).real() + yddag(s,v).real() * yd(w,t).imag()) * WC7R(Cqd1R, p, r, v, w))

					- 0.125 * ((yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd1I, p, v, s, w)
					+(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd1R, p, v, s, w))
					- 0.125 * (-(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd1I, r, v, t, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd1R, r, v, t, w))
					- 0.125 * ((yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd1I, s, v, p, w)
					+(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd1R, s, v, p, w))
					- 0.125 * (-(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd1I, t, v, r, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd1R, t, v, r, w))
					;

			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//Cqq3R
	for (d = 0; d < DWC6R; d++) {
		p = WC6R_indices[d][0];
		r = WC6R_indices[d][1];
		s = WC6R_indices[d][2];
		t = WC6R_indices[d][3];
		f[c] = -0.5 * ((yudyu(p,r).real() + yddyd(p,r).real()) * WC2R(CHq3R, s, t)
			- (yudyu(p,r).imag() + yddyd(p,r).imag()) * WC2I(CHq3I, s, t))
			- 0.5 * ((yudyu(s,t).real() + yddyd(s,t).real()) * WC2R(CHq3R, p, r)
			- (yudyu(s,t).imag() + yddyd(s,t).imag()) * WC2I(CHq3I, p, r)) //RGE 2
			+ ONE_SIXTH * g22 * (WC2R(CHq3R, s, t) * delta[p][r] + WC2R(CHq3R, p, r) * delta[s][t])
			- 3. * g32 * WC6R(Cqq3R, p, t, s, r)-(6. / NC) * g32 * WC6R(Cqq3R, p, r, s, t)
			- 6. * g22 * WC6R(Cqq3R, p, r, s, t) + 12. * Yq2 * g12 * WC6R(Cqq3R, p, r, s, t)
			+ 3. * g32 * WC6R(Cqq1R, p, t, s, r) + 3. * g22 * WC6R(Cqq1R, p, r, s, t) //RGE3
			;


		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index

			f[c] += (gq(p,v).real() * WC6R(Cqq3R, v, r, s, t) - gq(p,v).imag() * WC6I(Cqq3I, v, r, s, t))
				+ (gq(s,v).real() * WC6R(Cqq3R, p, r, v, t) - gq(s,v).imag() * WC6I(Cqq3I, p, r, v, t))
				+(WC6R(Cqq3R, p, v, s, t) * gq(v,r).real() - WC6I(Cqq3I, p, v, s, t) * gq(v,r).imag())
				+(WC6R(Cqq3R, p, r, s, v) * gq(v,t).real() - WC6I(Cqq3I, p, r, s, v) * gq(v,t).imag())
				//RGE 2
				+ ONE_THIRD * g22 * (
				WC7R(Clq3R, w, w, s, t) * delta[p][r] + WC7R(Clq3R, w, w, p, r) * delta[s][t])
				+ ONE_SIXTH * g22 * (
				+WC6R(Cqq1R, s, w, w, t) * delta[p][r] + WC6R(Cqq1R, w, t, s, w) * delta[p][r]
				+ WC6R(Cqq1R, p, w, w, r) * delta[s][t] + WC6R(Cqq1R, w, r, p, w) * delta[s][t])
				+ ONE_SIXTH * g32 * (
				WC6R(Cqq1R, p, w, w, t) * delta[r][s] + WC6R(Cqq1R, s, w, w, r) * delta[p][t]
				+ WC6R(Cqq1R, w, t, p, w) * delta[r][s] + WC6R(Cqq1R, w, r, s, w) * delta[p][t])
				+ ONE_THIRD * g22 * NC * (
				WC6R(Cqq3R, p, r, w, w) * delta[s][t] + WC6R(Cqq3R, s, t, w, w) * delta[p][r]
				+ WC6R(Cqq3R, w, w, s, t) * delta[p][r] + WC6R(Cqq3R, w, w, p, r) * delta[s][t])
				- ONE_SIXTH * g22 * (
				WC6R(Cqq3R, p, w, w, r) * delta[s][t] + WC6R(Cqq3R, s, w, w, t) * delta[p][r]
				+ WC6R(Cqq3R, w, r, p, w) * delta[s][t] + WC6R(Cqq3R, w, t, s, w) * delta[p][r])
				+ 0.5 * g32 * (
				WC6R(Cqq3R, p, w, w, t) * delta[r][s] + WC6R(Cqq3R, s, w, w, r) * delta[p][t]
				+ WC6R(Cqq3R, w, t, p, w) * delta[r][s] + WC6R(Cqq3R, w, r, s, w) * delta[p][t])
				+ (1./12.) * g32 * (
				+WC7R(Cqu8R, p, t, w, w) * delta[r][s] + WC7R(Cqu8R, s, r, w, w) * delta[p][t]
				+ WC7R(Cqd8R, p, t, w, w) * delta[r][s] + WC7R(Cqd8R, s, r, w, w) * delta[p][t])
				//RGE 3
				;
			for (w = 0; w < NG; w++) {
				f[c] += -0.125 * (
					+(yudag(p,v).real() * yu(w,t).real() - yudag(p,v).imag() * yu(w,t).imag()) * WC7R(Cqu8R, s, r, v, w)
					-(yudag(p,v).imag() * yu(w,t).real() + yudag(p,v).real() * yu(w,t).imag()) * WC7I(Cqu8I, s, r, v, w)
					+(yudag(s,v).real() * yu(w,r).real() - yudag(s,v).imag() * yu(w,r).imag()) * WC7R(Cqu8R, p, t, v, w)
					-(yudag(s,v).imag() * yu(w,r).real() + yudag(s,v).real() * yu(w,r).imag()) * WC7I(Cqu8I, p, t, v, w)
					)
					- 0.125 * (
					+(yddag(p,v).real() * yd(w,t).real() - yddag(p,v).imag() * yd(w,t).imag()) * WC7R(Cqd8R, s, r, v, w)
					-(yddag(p,v).imag() * yd(w,t).real() + yddag(p,v).real() * yd(w,t).imag()) * WC7I(Cqd8I, s, r, v, w)
					+(yddag(s,v).real() * yd(w,r).real() - yddag(s,v).imag() * yd(w,r).imag()) * WC7R(Cqd8R, p, t, v, w)
					-(yddag(s,v).imag() * yd(w,r).real() + yddag(s,v).real() * yd(w,r).imag()) * WC7I(Cqd8I, p, t, v, w)
					)

					-(0.0625 / NC)*(
					+(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd8R, p, v, s, w)
					-(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd8I, p, v, s, w)
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd8R, s, v, p, w)
					-(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd8I, s, v, p, w)
					)
					-(0.0625 / NC)*(
					+(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd8R, r, v, t, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd8I, r, v, t, w)
					+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd8R, t, v, r, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd8I, t, v, r, w)
					)

					- 0.0625 * (
					+(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd8R, s, v, p, w)
					-(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd8I, s, v, p, w)
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd8R, p, v, s, w)
					-(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd8I, p, v, s, w)
					)
					- 0.0625 * (
					+(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd8R, t, v, r, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd8I, t, v, r, w)
					+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd8R, r, v, t, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd8I, r, v, t, w)
					)

					+ 0.125 * (
					+(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd1R, p, v, s, w)
					-(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd1I, p, v, s, w)
					)
					+ 0.125 * (
					+(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd1R, r, v, t, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd1I, r, v, t, w)
					)
					+ 0.125 * (
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd1R, s, v, p, w)
					-(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd1I, s, v, p, w)
					)
					+ 0.125 * (
					+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd1R, t, v, r, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd1I, t, v, r, w)
					)
					//RGE 2
					;
			}
		}
		f[c] *= loop_factor;
		c++;


	}
	//Cqq3I
	for (d = 0; d < DWC6I; d++) {
		p = WC6I_indices[d][0];
		r = WC6I_indices[d][1];
		s = WC6I_indices[d][2];
		t = WC6I_indices[d][3];
		f[c] = -0.5 * ((yudyu(p,r).imag() + yddyd(p,r).imag()) * WC2R(CHq3R, s, t)
			+ (yudyu(p,r).real() + yddyd(p,r).real()) * WC2I(CHq3I, s, t))
			- 0.5 * ((yudyu(s,t).imag() + yddyd(s,t).imag()) * WC2R(CHq3R, p, r)
			+ (yudyu(s,t).real() + yddyd(s,t).real()) * WC2I(CHq3I, p, r)) //RGE 2
			+ ONE_SIXTH * g22 * (WC2I(CHq3I, s, t) * delta[p][r] + WC2I(CHq3I, p, r) * delta[s][t])
			- 3. * g32 * WC6I(Cqq3I, p, t, s, r)-(6. / NC) * g32 * WC6I(Cqq3I, p, r, s, t)
			- 6. * g22 * WC6I(Cqq3I, p, r, s, t) + 12. * Yq2 * g12 * WC6I(Cqq3I, p, r, s, t)
			+ 3. * g32 * WC6I(Cqq1I, p, t, s, r) + 3. * g22 * WC6I(Cqq1I, p, r, s, t) //RGE3
			;
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += (gq(p,v).imag() * WC6R(Cqq3R, v, r, s, t) + gq(p,v).real() * WC6I(Cqq3I, v, r, s, t))
				+ (gq(s,v).imag() * WC6R(Cqq3R, p, r, v, t) + gq(s,v).real() * WC6I(Cqq3I, p, r, v, t))
				+(WC6R(Cqq3R, p, v, s, t) * gq(v,r).imag() + WC6I(Cqq3I, p, v, s, t) * gq(v,r).real())
				+(WC6R(Cqq3R, p, r, s, v) * gq(v,t).imag() + WC6I(Cqq3I, p, r, s, v) * gq(v,t).real())
				//RGE 2
				+ ONE_THIRD * g22 * (
				WC7I(Clq3I, w, w, s, t) * delta[p][r] + WC7I(Clq3I, w, w, p, r) * delta[s][t])
				+ ONE_SIXTH * g22 * (
				+WC6I(Cqq1I, s, w, w, t) * delta[p][r] + WC6I(Cqq1I, w, t, s, w) * delta[p][r]
				+ WC6I(Cqq1I, p, w, w, r) * delta[s][t] + WC6I(Cqq1I, w, r, p, w) * delta[s][t])
				+ ONE_SIXTH * g32 * (
				WC6I(Cqq1I, p, w, w, t) * delta[r][s] + WC6I(Cqq1I, s, w, w, r) * delta[p][t]
				+ WC6I(Cqq1I, w, t, p, w) * delta[r][s] + WC6I(Cqq1I, w, r, s, w) * delta[p][t])
				+ ONE_THIRD * g22 * NC * (
				WC6I(Cqq3I, p, r, w, w) * delta[s][t] + WC6I(Cqq3I, s, t, w, w) * delta[p][r]
				+ WC6I(Cqq3I, w, w, s, t) * delta[p][r] + WC6I(Cqq3I, w, w, p, r) * delta[s][t])
				- ONE_SIXTH * g22 * (
				WC6I(Cqq3I, p, w, w, r) * delta[s][t] + WC6I(Cqq3I, s, w, w, t) * delta[p][r]
				+ WC6I(Cqq3I, w, r, p, w) * delta[s][t] + WC6I(Cqq3I, w, t, s, w) * delta[p][r])
				+ 0.5 * g32 * (
				WC6I(Cqq3I, p, w, w, t) * delta[r][s] + WC6I(Cqq3I, s, w, w, r) * delta[p][t]
				+ WC6I(Cqq3I, w, t, p, w) * delta[r][s] + WC6I(Cqq3I, w, r, s, w) * delta[p][t])
				+ (1. / 12.) * g32 * (
				+WC7I(Cqu8I, p, t, w, w) * delta[r][s] + WC7I(Cqu8I, s, r, w, w) * delta[p][t]
				+ WC7I(Cqd8I, p, t, w, w) * delta[r][s] + WC7I(Cqd8I, s, r, w, w) * delta[p][t])
				//RGE 3
				;

			for (w = 0; w < NG; w++) {
				f[c] += -0.125 * (
					+(yudag(p,v).real() * yu(w,t).real() - yudag(p,v).imag() * yu(w,t).imag()) * WC7I(Cqu8I, s, r, v, w)
					+(yudag(p,v).imag() * yu(w,t).real() + yudag(p,v).real() * yu(w,t).imag()) * WC7R(Cqu8R, s, r, v, w)
					+(yudag(s,v).real() * yu(w,r).real() - yudag(s,v).imag() * yu(w,r).imag()) * WC7I(Cqu8I, p, t, v, w)
					+(yudag(s,v).imag() * yu(w,r).real() + yudag(s,v).real() * yu(w,r).imag()) * WC7R(Cqu8R, p, t, v, w)
					)
					- 0.125 * (
					+(yddag(p,v).real() * yd(w,t).real() - yddag(p,v).imag() * yd(w,t).imag()) * WC7I(Cqd8I, s, r, v, w)
					+(yddag(p,v).imag() * yd(w,t).real() + yddag(p,v).real() * yd(w,t).imag()) * WC7R(Cqd8R, s, r, v, w)
					+(yddag(s,v).real() * yd(w,r).real() - yddag(s,v).imag() * yd(w,r).imag()) * WC7I(Cqd8I, p, t, v, w)
					+ (yddag(s,v).imag() * yd(w,r).real() + yddag(s,v).real() * yd(w,r).imag()) * WC7R(Cqd8R, p, t, v, w)
					)

					-(0.0625 / NC)*(
					+(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd8I, p, v, s, w)
					+(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd8R, p, v, s, w)
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd8I, s, v, p, w)
					+(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd8R, s, v, p, w)
					)
					-(0.0625 / NC)*(
					-(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd8I, r, v, t, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd8R, r, v, t, w)
					-(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd8I, t, v, r, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd8R, t, v, r, w)
					)

					- 0.0625 * (
					+(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd8I, s, v, p, w)
					+(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd8R, s, v, p, w)
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd8I, p, v, s, w)
					+(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd8R, p, v, s, w)
					)
					- 0.0625 * (
					-(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd8I, t, v, r, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd8R, t, v, r, w)
					-(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd8I, r, v, t, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd8R, r, v, t, w)
					)

					+ 0.125 * (
					+(yd(w,t).real() * yu(v,r).real() - yd(w,t).imag() * yu(v,r).imag()) * WC5(Cquqd1I, p, v, s, w)
					+(yd(w,t).imag() * yu(v,r).real() + yd(w,t).real() * yu(v,r).imag()) * WC5(Cquqd1R, p, v, s, w)
					)
					+ 0.125 * (
					-(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC5(Cquqd1I, r, v, t, w)
					+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC5(Cquqd1R, r, v, t, w)
					)
					+ 0.125 * (
					+(yd(w,r).real() * yu(v,t).real() - yd(w,r).imag() * yu(v,t).imag()) * WC5(Cquqd1I, s, v, p, w)
					+(yd(w,r).imag() * yu(v,t).real() + yd(w,r).real() * yu(v,t).imag()) * WC5(Cquqd1R, s, v, p, w)
					)
					+ 0.125 * (
					-(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC5(Cquqd1I, t, v, r, w)
					+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC5(Cquqd1R, t, v, r, w)
					)
					//RGE 2
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}

	//Clq1R
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			-(yedye(p,r).real() * WC2R(CHq1R, s, t) - yedye(p,r).imag() * WC2I(CHq1I, s, t))
			+(
			+(yudyu(s,t).real() - yddyd(s,t).real()) * WC2R(CHl1R, p, r)
			-(yudyu(s,t).imag() - yddyd(s,t).imag()) * WC2I(CHl1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYl * WC2R(CHq1R, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYq * WC2R(CHl1R, p, r) * delta[s][t]
			+ 12. * YlYq * g12 * WC7R(Clq1R, p, r, s, t)
			+ 9. * g22 * WC7R(Clq3R, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).real() * WC7R(Clq1R, v, r, s, t) - gl(p,v).imag() * WC7I(Clq1I, v, r, s, t))
				+(gq(s,v).real() * WC7R(Clq1R, p, r, v, t) - gq(s,v).imag() * WC7I(Clq1I, p, r, v, t))
				+(WC7R(Clq1R, p, v, s, t) * gl(v,r).real() - WC7I(Clq1I, p, v, s, t) * gl(v,r).imag())
				+(WC7R(Clq1R, p, r, s, v) * gq(v,t).real() - WC7I(Clq1I, p, r, s, v) * gq(v,t).imag())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YlYq * delta[s][t]*(
				+WC6R(CllR, p, r, w, w) + WC6R(CllR, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YlYq * delta[s][t]*(
				+WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * NC * Yq2 * WC7R(Clq1R, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * Yl2 * WC7R(Clq1R, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * NC * YlYq * delta[p][r]*(
				+WC6R(Cqq1R, s, t, w, w) + WC6R(Cqq1R, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YlYq * delta[p][r]*(
				+WC6R(Cqq1R, s, w, w, t) + WC6R(Cqq1R, w, t, s, w)
				)
				+ 4. * g12 * YlYq * delta[p][r]* (
				+WC6R(Cqq3R, s, w, w, t) + WC6R(Cqq3R, w, t, s, w)
				)
				+ FOUR_THIRDS * g12 * NC * YuYl * WC7R(Cqu1R, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYl * WC7R(Cqd1R, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYl * WC7R(CqeR, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYq * WC7R(CluR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YdYq * WC7R(CldR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * YeYq * WC7R(CleR, p, r, w, w) * delta[s][t]

				//RGE 3
				;
			for (w = 0; w < NG; w++) {
				f[c] += 0.25 * (
					+(yu(w,t).real() * ye(v,r).real() - yu(w,t).imag() * ye(v,r).imag()) * WC5(Clequ1R, p, v, s, w)
					-(yu(w,t).imag() * ye(v,r).real() + yu(w,t).real() * ye(v,r).imag()) * WC5(Clequ1I, p, v, s, w)
					)
					+ 0.25 * (
					+(yudag(s,w).real() * yedag(p,v).real() - yudag(s,w).imag() * yedag(p,v).imag()) * WC5(Clequ1R, r, v, t, w)
					+(yudag(s,w).imag() * yedag(p,v).real() + yudag(s,w).real() * yedag(p,v).imag()) * WC5(Clequ1I, r, v, t, w)
					)
					-(
					+(yudag(s,v).real() * yu(w,t).real() - yudag(s,v).imag() * yu(w,t).imag()) * WC7R(CluR, p, r, v, w)
					-(yudag(s,v).imag() * yu(w,t).real() + yudag(s,v).real() * yu(w,t).imag()) * WC7I(CluI, p, r, v, w)
					)
					-(
					+(yddag(s,v).real() * yd(w,t).real() - yddag(s,v).imag() * yd(w,t).imag()) * WC7R(CldR, p, r, v, w)
					-(yddag(s,v).imag() * yd(w,t).real() + yddag(s,v).real() * yd(w,t).imag()) * WC7I(CldI, p, r, v, w)
					)
					-(
					+(ye(p,v).real() * ye(w,r).real() - ye(p,v).imag() * ye(w,r).imag()) * WC7R(CqeR, s, t, v, w)
					-(ye(p,v).imag() * ye(w,r).real() + ye(p,v).real() * ye(w,r).imag()) * WC7I(CqeI, s, t, v, w)
					)
					+ 0.25 * (
					+(yddag(s,w).real() * ye(v,r).real() - yddag(s,w).imag() * ye(v,r).imag()) * WC5(CledqR, p, v, w, t)
					-(yddag(s,w).imag() * ye(v,r).real() + yddag(s,w).real() * ye(v,r).imag()) * WC5(CledqI, p, v, w, t)
					+(yedag(p,v).real() * yd(w,t).real() - yedag(p,v).imag() * yd(w,t).imag()) * WC5(CledqR, r, v, w, s)
					+(yedag(p,v).imag() * yd(w,t).real() + yedag(p,v).real() * yd(w,t).imag()) * WC5(CledqI, r, v, w, s)
					)
					- 3. * (
					+(ye(v,r).real() * yu(w,t).real() - ye(v,r).imag() * yu(w,t).imag()) * WC5(Clequ3R, p, v, s, w)
					-(ye(v,r).imag() * yu(w,t).real() + ye(v,r).real() * yu(w,t).imag()) * WC5(Clequ3I, p, v, s, w)
					+(yedag(p,v).real() * yudag(s,w).real() - yedag(p,v).imag() * yudag(s,w).imag()) * WC5(Clequ3R, r, v, t, w)
					+(yedag(p,v).imag() * yudag(s,w).real() + yedag(p,v).real() * yudag(s,w).imag()) * WC5(Clequ3I, r, v, t, w)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}


	//Clq1I
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			-(yedye(p,r).imag() * WC2R(CHq1R, s, t) + yedye(p,r).real() * WC2I(CHq1I, s, t))
			+(
			+(yudyu(s,t).imag() - yddyd(s,t).imag()) * WC2R(CHl1R, p, r)
			+(yudyu(s,t).real() - yddyd(s,t).real()) * WC2I(CHl1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYl * WC2I(CHq1I, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYq * WC2I(CHl1I, p, r) * delta[s][t]
			+ 12. * YlYq * g12 * WC7I(Clq1I, p, r, s, t)
			+ 9. * g22 * WC7I(Clq3I, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).imag() * WC7R(Clq1R, v, r, s, t) + gl(p,v).real() * WC7I(Clq1I, v, r, s, t))
				+(gq(s,v).imag() * WC7R(Clq1R, p, r, v, t) + gq(s,v).real() * WC7I(Clq1I, p, r, v, t))
				+(WC7R(Clq1R, p, v, s, t) * gl(v,r).imag() + WC7I(Clq1I, p, v, s, t) * gl(v,r).real())
				+(WC7R(Clq1R, p, r, s, v) * gq(v,t).imag() + WC7I(Clq1I, p, r, s, v) * gq(v,t).real())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YlYq * delta[s][t]*(
				+WC6I(CllI, p, r, w, w) + WC6I(CllI, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YlYq * delta[s][t]*(
				+WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * NC * Yq2 * WC7I(Clq1I, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * Yl2 * WC7I(Clq1I, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * NC * YlYq * delta[p][r]*(
				+WC6I(Cqq1I, s, t, w, w) + WC6I(Cqq1I, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YlYq * delta[p][r]*(
				+WC6I(Cqq1I, s, w, w, t) + WC6I(Cqq1I, w, t, s, w)
				)
				+ 4. * g12 * YlYq * delta[p][r]* (
				+WC6I(Cqq3I, s, w, w, t) + WC6I(Cqq3I, w, t, s, w)
				)
				+ FOUR_THIRDS * g12 * NC * YuYl * WC7I(Cqu1I, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYl * WC7I(Cqd1I, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYl * WC7I(CqeI, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYq * WC7I(CluI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YdYq * WC7I(CldI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * YeYq * WC7I(CleI, p, r, w, w) * delta[s][t]

				//RGE 3
				;
			for (w = 0; w < NG; w++) {
				f[c] += 0.25 * (
					+(yu(w,t).real() * ye(v,r).real() - yu(w,t).imag() * ye(v,r).imag()) * WC5(Clequ1I, p, v, s, w)
					+(yu(w,t).imag() * ye(v,r).real() + yu(w,t).real() * ye(v,r).imag()) * WC5(Clequ1R, p, v, s, w)
					)
					+ 0.25 * (
					-(yudag(s,w).real() * yedag(p,v).real() - yudag(s,w).imag() * yedag(p,v).imag()) * WC5(Clequ1I, r, v, t, w)
					+(yudag(s,w).imag() * yedag(p,v).real() + yudag(s,w).real() * yedag(p,v).imag()) * WC5(Clequ1R, r, v, t, w)
					)
					-(
					+(yudag(s,v).real() * yu(w,t).real() - yudag(s,v).imag() * yu(w,t).imag()) * WC7I(CluI, p, r, v, w)
					+(yudag(s,v).imag() * yu(w,t).real() + yudag(s,v).real() * yu(w,t).imag()) * WC7R(CluR, p, r, v, w)
					)
					-(
					+(yddag(s,v).real() * yd(w,t).real() - yddag(s,v).imag() * yd(w,t).imag()) * WC7I(CldI, p, r, v, w)
					+(yddag(s,v).imag() * yd(w,t).real() + yddag(s,v).real() * yd(w,t).imag()) * WC7R(CldR, p, r, v, w)
					)
					-(
					+(ye(p,v).real() * ye(w,r).real() - ye(p,v).imag() * ye(w,r).imag()) * WC7I(CqeI, s, t, v, w)
					+(ye(p,v).imag() * ye(w,r).real() + ye(p,v).real() * ye(w,r).imag()) * WC7R(CqeR, s, t, v, w)
					)
					+ 0.25 * (
					+(yddag(s,w).real() * ye(v,r).real() - yddag(s,w).imag() * ye(v,r).imag()) * WC5(CledqI, p, v, w, t)
					+(yddag(s,w).imag() * ye(v,r).real() + yddag(s,w).real() * ye(v,r).imag()) * WC5(CledqR, p, v, w, t)
					-(yedag(p,v).real() * yd(w,t).real() - yedag(p,v).imag() * yd(w,t).imag()) * WC5(CledqI, r, v, w, s)
					+(yedag(p,v).imag() * yd(w,t).real() + yedag(p,v).real() * yd(w,t).imag()) * WC5(CledqR, r, v, w, s)
					)
					- 3. * (
					+(ye(v,r).real() * yu(w,t).real() - ye(v,r).imag() * yu(w,t).imag()) * WC5(Clequ3I, p, v, s, w)
					+(ye(v,r).imag() * yu(w,t).real() + ye(v,r).real() * yu(w,t).imag()) * WC5(Clequ3R, p, v, s, w)
					-(yedag(p,v).real() * yudag(s,w).real() - yedag(p,v).imag() * yudag(s,w).imag()) * WC5(Clequ3I, r, v, t, w)
					+(yedag(p,v).imag() * yudag(s,w).real() + yedag(p,v).real() * yudag(s,w).imag()) * WC5(Clequ3R, r, v, t, w)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}



	//Clq3R
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			-(yedye(p,r).real() * WC2R(CHq3R, s, t) - yedye(p,r).imag() * WC2I(CHq3I, s, t))
			-(
			+(yudyu(s,t).real() + yddyd(s,t).real()) * WC2R(CHl3R, p, r)
			-(yudyu(s,t).imag() + yddyd(s,t).imag()) * WC2I(CHl3I, p, r)
			)
			//RGE 2
			+ ONE_THIRD * g22 * WC2R(CHq3R, s, t) * delta[p][r]
			+ ONE_THIRD * g22 * WC2R(CHl3R, p, r) * delta[s][t]
			+ 3. * g22 * WC7R(Clq1R, p, r, s, t)
			- 6. * (g22 - 2. * YlYq * g12) * WC7R(Clq3R, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).real() * WC7R(Clq3R, v, r, s, t) - gl(p,v).imag() * WC7I(Clq3I, v, r, s, t))
				+(gq(s,v).real() * WC7R(Clq3R, p, r, v, t) - gq(s,v).imag() * WC7I(Clq3I, p, r, v, t))
				+(WC7R(Clq3R, p, v, s, t) * gl(v,r).real() - WC7I(Clq3I, p, v, s, t) * gl(v,r).imag())
				+(WC7R(Clq3R, p, r, s, v) * gq(v,t).real() - WC7I(Clq3I, p, r, s, v) * gq(v,t).imag())
				//RGE 2
				+ TWO_THIRDS * g22 * NC * WC7R(Clq3R, p, r, w, w) * delta[s][t]
				+ TWO_THIRDS * g22 * WC7R(Clq3R, w, w, s, t) * delta[p][r]
				+ ONE_THIRD * g22 * delta[p][r] * (
				+WC6R(Cqq1R, s, w, w, t) + WC6R(Cqq1R, w, t, s, w)
				)
				+ TWO_THIRDS * g22 * NC * delta[p][r]*(
				+WC6R(Cqq3R, s, t, w, w) + WC6R(Cqq3R, w, w, s, t)
				)
				- ONE_THIRD * g22 * delta[p][r]*(
				+WC6R(Cqq3R, s, w, w, t) + WC6R(Cqq3R, w, t, s, w)
				)
				+ ONE_THIRD * g22 * delta[s][t]*(
				+WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					-0.25 * (
					+(yu(w,t).real() * ye(v,r).real() - yu(w,t).imag() * ye(v,r).imag()) * WC5(Clequ1R, p, v, s, w)
					-(yu(w,t).imag() * ye(v,r).real() + yu(w,t).real() * ye(v,r).imag()) * WC5(Clequ1I, p, v, s, w)
					)
					- 0.25 * (
					+(yudag(s,w).real() * yedag(p,v).real() - yudag(s,w).imag() * yedag(p,v).imag()) * WC5(Clequ1R, r, v, t, w)
					+(yudag(s,w).imag() * yedag(p,v).real() + yudag(s,w).real() * yedag(p,v).imag()) * WC5(Clequ1I, r, v, t, w)
					)
					+ 0.25 * (
					+(yddag(s,w).real() * ye(v,r).real() - yddag(s,w).imag() * ye(v,r).imag()) * WC5(CledqR, p, v, w, t)
					-(yddag(s,w).imag() * ye(v,r).real() + yddag(s,w).real() * ye(v,r).imag()) * WC5(CledqI, p, v, w, t)
					+(yedag(p,v).real() * yd(w,t).real() - yedag(p,v).imag() * yd(w,t).imag()) * WC5(CledqR, r, v, w, s)
					+(yedag(p,v).imag() * yd(w,t).real() + yedag(p,v).real() * yd(w,t).imag()) * WC5(CledqI, r, v, w, s)
					)
					+ 3. * (
					+(ye(v,r).real() * yu(w,t).real() - ye(v,r).imag() * yu(w,t).imag()) * WC5(Clequ3R, p, v, s, w)
					-(ye(v,r).imag() * yu(w,t).real() + ye(v,r).real() * yu(w,t).imag()) * WC5(Clequ3I, p, v, s, w)
					+(yedag(p,v).real() * yudag(s,w).real() - yedag(p,v).imag() * yudag(s,w).imag()) * WC5(Clequ3R, r, v, t, w)
					+(yedag(p,v).imag() * yudag(s,w).real() + yedag(p,v).real() * yudag(s,w).imag()) * WC5(Clequ3I, r, v, t, w)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}

	//Clq3I
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			-(yedye(p,r).imag() * WC2R(CHq3R, s, t) + yedye(p,r).real() * WC2I(CHq3I, s, t))
			-(
			+(yudyu(s,t).imag() + yddyd(s,t).imag()) * WC2R(CHl3R, p, r)
			+(yudyu(s,t).real() + yddyd(s,t).real()) * WC2I(CHl3I, p, r)
			)
			//RGE 2
			+ ONE_THIRD * g22 * WC2I(CHq3I, s, t) * delta[p][r]
			+ ONE_THIRD * g22 * WC2I(CHl3I, p, r) * delta[s][t]
			+ 3. * g22 * WC7I(Clq1I, p, r, s, t)
			- 6. * (g22 - 2. * YlYq * g12) * WC7I(Clq3I, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).imag() * WC7R(Clq3R, v, r, s, t) + gl(p,v).real() * WC7I(Clq3I, v, r, s, t))
				+(gq(s,v).imag() * WC7R(Clq3R, p, r, v, t) + gq(s,v).real() * WC7I(Clq3I, p, r, v, t))
				+(WC7R(Clq3R, p, v, s, t) * gl(v,r).imag() + WC7I(Clq3I, p, v, s, t) * gl(v,r).real())
				+(WC7R(Clq3R, p, r, s, v) * gq(v,t).imag() + WC7I(Clq3I, p, r, s, v) * gq(v,t).real())
				//RGE 2
				+ TWO_THIRDS * g22 * NC * WC7I(Clq3I, p, r, w, w) * delta[s][t]
				+ TWO_THIRDS * g22 * WC7I(Clq3I, w, w, s, t) * delta[p][r]
				+ ONE_THIRD * g22 * delta[p][r] * (
				+WC6I(Cqq1I, s, w, w, t) + WC6I(Cqq1I, w, t, s, w)
				)
				+ TWO_THIRDS * g22 * NC * delta[p][r]*(
				+WC6I(Cqq3I, s, t, w, w) + WC6I(Cqq3I, w, w, s, t)
				)
				- ONE_THIRD * g22 * delta[p][r]*(
				+WC6I(Cqq3I, s, w, w, t) + WC6I(Cqq3I, w, t, s, w)
				)
				+ ONE_THIRD * g22 * delta[s][t]*(
				+WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					-0.25 * (
					+(yu(w,t).real() * ye(v,r).real() - yu(w,t).imag() * ye(v,r).imag()) * WC5(Clequ1I, p, v, s, w)
					+(yu(w,t).imag() * ye(v,r).real() + yu(w,t).real() * ye(v,r).imag()) * WC5(Clequ1R, p, v, s, w)
					)
					- 0.25 * (
					-(yudag(s,w).real() * yedag(p,v).real() - yudag(s,w).imag() * yedag(p,v).imag()) * WC5(Clequ1I, r, v, t, w)
					+(yudag(s,w).imag() * yedag(p,v).real() + yudag(s,w).real() * yedag(p,v).imag()) * WC5(Clequ1R, r, v, t, w)
					)
					+ 0.25 * (
					+(yddag(s,w).real() * ye(v,r).real() - yddag(s,w).imag() * ye(v,r).imag()) * WC5(CledqI, p, v, w, t)
					+(yddag(s,w).imag() * ye(v,r).real() + yddag(s,w).real() * ye(v,r).imag()) * WC5(CledqR, p, v, w, t)
					-(yedag(p,v).real() * yd(w,t).real() - yedag(p,v).imag() * yd(w,t).imag()) * WC5(CledqI, r, v, w, s)
					+(yedag(p,v).imag() * yd(w,t).real() + yedag(p,v).real() * yd(w,t).imag()) * WC5(CledqR, r, v, w, s)
					)
					+ 3. * (
					+(ye(v,r).real() * yu(w,t).real() - ye(v,r).imag() * yu(w,t).imag()) * WC5(Clequ3I, p, v, s, w)
					+(ye(v,r).imag() * yu(w,t).real() + ye(v,r).real() * yu(w,t).imag()) * WC5(Clequ3R, p, v, s, w)
					-(yedag(p,v).real() * yudag(s,w).real() - yedag(p,v).imag() * yudag(s,w).imag()) * WC5(Clequ3I, r, v, t, w)
					+(yedag(p,v).imag() * yudag(s,w).real() + yedag(p,v).real() * yudag(s,w).imag()) * WC5(Clequ3R, r, v, t, w)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}




	//----------------------RGE SMEFT class 8_RRRR----------------------
	//CeeR
	for (d = 0; d < DWC8R; d++) {
		p = WC8R_indices[d][0];
		r = WC8R_indices[d][1];
		s = WC8R_indices[d][2];
		t = WC8R_indices[d][3];
		f[c] = (ge(p,r).real() * WC2R(CHeR, s, t) - ge(p,r).imag() * WC2I(CHeI, s, t))
			+(ge(s,t).real() * WC2R(CHeR, p, r) - ge(s,t).imag() * WC2I(CHeI, p, r)) //RGE 2
			+ TWO_THIRDS * g12 * YhYe * WC2R(CHeR, s, t) * delta[p][r]
			+ TWO_THIRDS * g12 * YhYe * WC2R(CHeR, p, r) * delta[s][t]
			+ 12. * Ye2 * g12 * WC8R(CeeR, p, r, s, t); //RGE 3

		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index

			f[c] += +(ge(p,v).real() * WC8R(CeeR, v, r, s, t) - ge(p,v).imag() * WC8I(CeeI, v, r, s, t))
				+(ge(s,v).real() * WC8R(CeeR, p, r, v, t) - ge(s,v).imag() * WC8I(CeeI, p, r, v, t))
				+(WC8R(CeeR, p, v, s, t) * ge(v,r).real() - WC8I(CeeI, p, v, s, t) * ge(v,r).imag())
				+(WC8R(CeeR, p, r, s, v) * ge(v,t).real() - WC8I(CeeI, p, r, s, v) * ge(v,t).imag())
				//RGE 2
				+ FOUR_THIRDS * g12 * YeYl * (
				WC7R(CleR, w, w, p, r) * delta[s][t] + WC7R(CleR, w, w, s, t) * delta[p][r])
				+ FOUR_THIRDS * g12 * NC * YeYq * (
				WC7R(CqeR, w, w, p, r) * delta[s][t] + WC7R(CqeR, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YuYe * (
				WC7R(CeuR, p, r, w, w) * delta[s][t] + WC7R(CeuR, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YdYe * (
				WC7R(CedR, p, r, w, w) * delta[s][t] + WC7R(CedR, s, t, w, w) * delta[p][r])

				+ TWO_THIRDS * g12 * Ye2 * (
				+WC8R(CeeR, p, r, w, w) * delta[s][t] + WC8R(CeeR, s, t, w, w) * delta[p][r]
				+ WC8R(CeeR, w, w, p, r) * delta[s][t] + WC8R(CeeR, w, w, s, t) * delta[p][r]
				+ WC8R(CeeR, p, w, w, r) * delta[s][t] + WC8R(CeeR, s, w, w, t) * delta[p][r]
				+ WC8R(CeeR, w, t, s, w) * delta[p][r] + WC8R(CeeR, w, r, p, w) * delta[s][t])
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -(
					(yedag(w,r).real() * ye(p,v).real() - yedag(w,r).imag() * ye(p,v).imag()) * WC7R(CleR, v, w, s, t)
					-(yedag(w,r).imag() * ye(p,v).real() + yedag(w,r).real() * ye(p,v).imag()) * WC7I(CleI, v, w, s, t)
					)
					-(
					(yedag(w,t).real() * ye(s,v).real() - yedag(w,t).imag() * ye(s,v).imag()) * WC7R(CleR, v, w, p, r)
					-(yedag(w,t).imag() * ye(s,v).real() + yedag(w,t).real() * ye(s,v).imag()) * WC7I(CleI, v, w, p, r)
					)
					//RGE 2 
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CeeI
	for (d = 0; d < DWC8I; d++) {
		p = WC8I_indices[d][0];
		r = WC8I_indices[d][1];
		s = WC8I_indices[d][2];
		t = WC8I_indices[d][3];

		f[c] = (ge(p,r).imag() * WC2R(CHeR, s, t) + ge(p,r).real() * WC2I(CHeI, s, t))
			+(ge(s,t).imag() * WC2R(CHeR, p, r) + ge(s,t).real() * WC2I(CHeI, p, r)) //RGE 2
			+ TWO_THIRDS * g12 * YhYe * WC2I(CHeI, s, t) * delta[p][r]
			+ TWO_THIRDS * g12 * YhYe * WC2I(CHeI, p, r) * delta[s][t]
			+ 12. * Ye2 * g12 * WC8I(CeeI, p, r, s, t); //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index

			f[c] += +(ge(p,v).imag() * WC8R(CeeR, v, r, s, t) + ge(p,v).real() * WC8I(CeeI, v, r, s, t))
				+(ge(s,v).imag() * WC8R(CeeR, p, r, v, t) + ge(s,v).real() * WC8I(CeeI, p, r, v, t))
				+(WC8R(CeeR, p, v, s, t) * ge(v,r).imag() + WC8I(CeeI, p, v, s, t) * ge(v,r).real())
				+(WC8R(CeeR, p, r, s, v) * ge(v,t).imag() + WC8I(CeeI, p, r, s, v) * ge(v,t).real())
				//RGE 2
				+ FOUR_THIRDS * g12 * YeYl * (
				WC7I(CleI, w, w, p, r) * delta[s][t] + WC7I(CleI, w, w, s, t) * delta[p][r])
				+ FOUR_THIRDS * g12 * NC * YeYq * (
				WC7I(CqeI, w, w, p, r) * delta[s][t] + WC7I(CqeI, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YuYe * (
				WC7I(CeuI, p, r, w, w) * delta[s][t] + WC7I(CeuI, s, t, w, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YdYe * (
				WC7I(CedI, p, r, w, w) * delta[s][t] + WC7I(CedI, s, t, w, w) * delta[p][r])

				+ TWO_THIRDS * g12 * Ye2 * (
				+WC8I(CeeI, p, r, w, w) * delta[s][t] + WC8I(CeeI, s, t, w, w) * delta[p][r]
				+ WC8I(CeeI, w, w, p, r) * delta[s][t] + WC8I(CeeI, w, w, s, t) * delta[p][r]
				+ WC8I(CeeI, p, w, w, r) * delta[s][t] + WC8I(CeeI, s, w, w, t) * delta[p][r]
				+ WC8I(CeeI, w, t, s, w) * delta[p][r] + WC8I(CeeI, w, r, p, w) * delta[s][t])
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -(
					(yedag(w,r).real() * ye(p,v).real() - yedag(w,r).imag() * ye(p,v).imag()) * WC7I(CleI, v, w, s, t)
					+(yedag(w,r).imag() * ye(p,v).real() + yedag(w,r).real() * ye(p,v).imag()) * WC7R(CleR, v, w, s, t)
					)
					-(
					(yedag(w,t).real() * ye(s,v).real() - yedag(w,t).imag() * ye(s,v).imag()) * WC7I(CleI, v, w, p, r)
					+(yedag(w,t).imag() * ye(s,v).real() + yedag(w,t).real() * ye(s,v).imag()) * WC7R(CleR, v, w, p, r)
					)
					//RGE 2 
					;
			}

		}
		f[c] *= loop_factor;
		c++;
	}
	//CuuR
	for (d = 0; d < DWC6R; d++) {
		p = WC6R_indices[d][0];
		r = WC6R_indices[d][1];
		s = WC6R_indices[d][2];
		t = WC6R_indices[d][3];
		f[c] = -(gu(p,r).real() * WC2R(CHuR, s, t) - gu(p,r).imag() * WC2I(CHuI, s, t))
			-(gu(s,t).real() * WC2R(CHuR, p, r) - gu(s,t).imag() * WC2I(CHuI, p, r))
			//RGE 2
			+ TWO_THIRDS * g12 * YhYu *
			(WC2R(CHuR, s, t) * delta[p][r] + WC2R(CHuR, p, r) * delta[s][t])
			+ 6. * g32 * WC6R(CuuR, p, t, s, r) - 6. * (g32 / NC - 2. * Yu2 * g12) * WC6R(CuuR, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += +(gu(p,v).real() * WC6R(CuuR, v, r, s, t) - gu(p,v).imag() * WC6I(CuuI, v, r, s, t))
				+(gu(s,v).real() * WC6R(CuuR, p, r, v, t) - gu(s,v).imag() * WC6I(CuuI, p, r, v, t))
				+(WC6R(CuuR, p, v, s, t) * gu(v,r).real() - WC6I(CuuI, p, v, s, t) * gu(v,r).imag())
				+(WC6R(CuuR, p, r, s, v) * gu(v,t).real() - WC6I(CuuI, p, r, s, v) * gu(v,t).imag())
				//RGE 2
				+ TWO_THIRDS * g12 * YuYe * (
				WC7R(CeuR, w, w, s, t) * delta[p][r] + WC7R(CeuR, w, w, p, r) * delta[s][t])
				+ FOUR_THIRDS * g12 * YuYl * (
				WC7R(CluR, w, w, p, r) * delta[s][t] + WC7R(CluR, w, w, s, t) * delta[p][r])
				+ FOUR_THIRDS * g12 * NC * YuYq * (
				WC7R(Cqu1R, w, w, s, t) * delta[p][r] + WC7R(Cqu1R, w, w, p, r) * delta[s][t])
				+ ONE_THIRD * g32 * (
				WC7R(Cqu8R, w, w, p, t) * delta[r][s] + WC7R(Cqu8R, w, w, s, r) * delta[p][t])
				-(ONE_THIRD / NC) * g32 * (
				WC7R(Cqu8R, w, w, s, t) * delta[p][r] + WC7R(Cqu8R, w, w, p, r) * delta[s][t])
				+ TWO_THIRDS * g12 * NC * Yu2 * (
				WC6R(CuuR, p, r, w, w) * delta[s][t] + WC6R(CuuR, s, t, w, w) * delta[p][r]
				+ WC6R(CuuR, w, w, p, r) * delta[s][t] + WC6R(CuuR, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * Yu2 * (
				WC6R(CuuR, p, w, w, r) * delta[s][t] + WC6R(CuuR, s, w, w, t) * delta[p][r]
				+ WC6R(CuuR, w, r, p, w) * delta[s][t] + WC6R(CuuR, w, t, s, w) * delta[p][r])
				+ ONE_THIRD * g32 * (
				WC6R(CuuR, p, w, w, t) * delta[r][s] + WC6R(CuuR, s, w, w, r) * delta[p][t]
				+ WC6R(CuuR, w, t, p, w) * delta[r][s] + WC6R(CuuR, w, r, s, w) * delta[p][t])
				-(ONE_THIRD / NC) * g32 * (
				WC6R(CuuR, p, w, w, r) * delta[s][t] + WC6R(CuuR, s, w, w, t) * delta[p][r]
				+ WC6R(CuuR, w, r, p, w) * delta[s][t] + WC6R(CuuR, w, t, s, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YuYd * (
				WC7R(Cud1R, p, r, w, w) * delta[s][t] + WC7R(Cud1R, s, t, w, w) * delta[p][r])
				+ ONE_SIXTH * g32 * (
				WC7R(Cud8R, p, t, w, w) * delta[r][s] + WC7R(Cud8R, s, r, w, w) * delta[p][t])
				-(ONE_SIXTH / NC) * g32 *
				(WC7R(Cud8R, p, r, w, w) * delta[s][t] + WC7R(Cud8R, s, t, w, w) * delta[p][r])
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -(
					+(yudag(w,r).real() * yu(p,v).real() - yudag(w,r).imag() * yu(p,v).imag()) * WC7R(Cqu1R, v, w, s, t)
					-(yudag(w,r).imag() * yu(p,v).real() + yudag(w,r).real() * yu(p,v).imag()) * WC7I(Cqu1I, v, w, s, t))
					-(
					+(yudag(w,t).real() * yu(s,v).real() - yudag(w,t).imag() * yu(s,v).imag()) * WC7R(Cqu1R, v, w, p, r)
					-(yudag(w,t).imag() * yu(s,v).real() + yudag(w,t).real() * yu(s,v).imag()) * WC7I(Cqu1I, v, w, p, r))

					+(0.5 / NC)*(
					+(yu(p,v).real() * yudag(w,r).real() - yu(p,v).imag() * yudag(w,r).imag()) * WC7R(Cqu8R, v, w, s, t)
					-(yu(p,v).imag() * yudag(w,r).real() + yu(p,v).real() * yudag(w,r).imag()) * WC7I(Cqu8I, v, w, s, t))
					+(0.5 / NC)*(
					+(yu(s,v).real() * yudag(w,t).real() - yu(s,v).imag() * yudag(w,t).imag()) * WC7R(Cqu8R, v, w, p, r)
					-(yu(s,v).imag() * yudag(w,t).real() + yu(s,v).real() * yudag(w,t).imag()) * WC7I(Cqu8I, v, w, p, r))

					- 0.5 * (
					+(yudag(w,r).real() * yu(s,v).real() - yudag(w,r).imag() * yu(s,v).imag()) * WC7R(Cqu8R, v, w, p, t)
					-(yudag(w,r).imag() * yu(s,v).real() + yudag(w,r).real() * yu(s,v).imag()) * WC7I(Cqu8I, v, w, p, t))
					- 0.5 * (
					+(yudag(w,t).real() * yu(p,v).real() - yudag(w,t).imag() * yu(p,v).imag()) * WC7R(Cqu8R, v, w, s, r)
					-(yudag(w,t).imag() * yu(p,v).real() + yudag(w,t).real() * yu(p,v).imag()) * WC7I(Cqu8I, v, w, s, r))

					//RGE 2
					;
			}
		}

		f[c] *= loop_factor;
		c++;
	}
	//CuuI
	for (d = 0; d < DWC6I; d++) {
		p = WC6I_indices[d][0];
		r = WC6I_indices[d][1];
		s = WC6I_indices[d][2];
		t = WC6I_indices[d][3];
		f[c] = -(gu(p,r).imag() * WC2R(CHuR, s, t) + gu(p,r).real() * WC2I(CHuI, s, t))
			-(gu(s,t).imag() * WC2R(CHuR, p, r) + gu(s,t).real() * WC2I(CHuI, p, r))
			//RGE 2
			+ TWO_THIRDS * g12 * YhYu *
			(WC2I(CHuI, s, t) * delta[p][r] + WC2I(CHuI, p, r) * delta[s][t])
			+ 6. * g32 * WC6I(CuuI, p, t, s, r) - 6. * (g32 / NC - 2. * Yu2 * g12) * WC6I(CuuI, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += +(gu(p,v).imag() * WC6R(CuuR, v, r, s, t) + gu(p,v).real() * WC6I(CuuI, v, r, s, t))
				+(gu(s,v).imag() * WC6R(CuuR, p, r, v, t) + gu(s,v).real() * WC6I(CuuI, p, r, v, t))
				+(WC6R(CuuR, p, v, s, t) * gu(v,r).imag() + WC6I(CuuI, p, v, s, t) * gu(v,r).real())
				+(WC6R(CuuR, p, r, s, v) * gu(v,t).imag() + WC6I(CuuI, p, r, s, v) * gu(v,t).real())
				//RGE 2
				+ TWO_THIRDS * g12 * YuYe * (
				WC7I(CeuI, w, w, s, t) * delta[p][r] + WC7I(CeuI, w, w, p, r) * delta[s][t])
				+ FOUR_THIRDS * g12 * YuYl * (
				WC7I(CluI, w, w, p, r) * delta[s][t] + WC7I(CluI, w, w, s, t) * delta[p][r])
				+ FOUR_THIRDS * g12 * NC * YuYq * (
				WC7I(Cqu1I, w, w, s, t) * delta[p][r] + WC7I(Cqu1I, w, w, p, r) * delta[s][t])
				+ ONE_THIRD * g32 * (
				WC7I(Cqu8I, w, w, p, t) * delta[r][s] + WC7I(Cqu8I, w, w, s, r) * delta[p][t])
				-(ONE_THIRD / NC) * g32 * (
				WC7I(Cqu8I, w, w, s, t) * delta[p][r] + WC7I(Cqu8I, w, w, p, r) * delta[s][t])
				+ TWO_THIRDS * g12 * NC * Yu2 * (
				WC6I(CuuI, p, r, w, w) * delta[s][t] + WC6I(CuuI, s, t, w, w) * delta[p][r]
				+ WC6I(CuuI, w, w, p, r) * delta[s][t] + WC6I(CuuI, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * Yu2 * (
				WC6I(CuuI, p, w, w, r) * delta[s][t] + WC6I(CuuI, s, w, w, t) * delta[p][r]
				+ WC6I(CuuI, w, r, p, w) * delta[s][t] + WC6I(CuuI, w, t, s, w) * delta[p][r])
				+ ONE_THIRD * g32 * (
				WC6I(CuuI, p, w, w, t) * delta[r][s] + WC6I(CuuI, s, w, w, r) * delta[p][t]
				+ WC6I(CuuI, w, t, p, w) * delta[r][s] + WC6I(CuuI, w, r, s, w) * delta[p][t])
				-(ONE_THIRD / NC) * g32 * (
				WC6I(CuuI, p, w, w, r) * delta[s][t] + WC6I(CuuI, s, w, w, t) * delta[p][r]
				+ WC6I(CuuI, w, r, p, w) * delta[s][t] + WC6I(CuuI, w, t, s, w) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YuYd * (
				WC7I(Cud1I, p, r, w, w) * delta[s][t] + WC7I(Cud1I, s, t, w, w) * delta[p][r])
				+ ONE_SIXTH * g32 * (
				WC7I(Cud8I, p, t, w, w) * delta[r][s] + WC7I(Cud8I, s, r, w, w) * delta[p][t])
				-(ONE_SIXTH / NC) * g32 *
				(WC7I(Cud8I, p, r, w, w) * delta[s][t] + WC7I(Cud8I, s, t, w, w) * delta[p][r])
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -(
					+(yudag(w,r).real() * yu(p,v).real() - yudag(w,r).imag() * yu(p,v).imag()) * WC7I(Cqu1I, v, w, s, t)
					+(yudag(w,r).imag() * yu(p,v).real() + yudag(w,r).real() * yu(p,v).imag()) * WC7R(Cqu1R, v, w, s, t))
					-(
					+(yudag(w,t).real() * yu(s,v).real() - yudag(w,t).imag() * yu(s,v).imag()) * WC7I(Cqu1I, v, w, p, r)
					+(yudag(w,t).imag() * yu(s,v).real() + yudag(w,t).real() * yu(s,v).imag()) * WC7R(Cqu1R, v, w, p, r))

					+(0.5 / NC)*(
					+(yu(p,v).real() * yudag(w,r).real() - yu(p,v).imag() * yudag(w,r).imag()) * WC7I(Cqu8I, v, w, s, t)
					+(yu(p,v).imag() * yudag(w,r).real() + yu(p,v).real() * yudag(w,r).imag()) * WC7R(Cqu8R, v, w, s, t))
					+(0.5 / NC)*(
					+(yu(s,v).real() * yudag(w,t).real() - yu(s,v).imag() * yudag(w,t).imag()) * WC7I(Cqu8I, v, w, p, r)
					+(yu(s,v).imag() * yudag(w,t).real() + yu(s,v).real() * yudag(w,t).imag()) * WC7R(Cqu8R, v, w, p, r))

					- 0.5 * (
					+(yudag(w,r).real() * yu(s,v).real() - yudag(w,r).imag() * yu(s,v).imag()) * WC7I(Cqu8I, v, w, p, t)
					+(yudag(w,r).imag() * yu(s,v).real() + yudag(w,r).real() * yu(s,v).imag()) * WC7R(Cqu8R, v, w, p, t))
					- 0.5 * (
					+(yudag(w,t).real() * yu(p,v).real() - yudag(w,t).imag() * yu(p,v).imag()) * WC7I(Cqu8I, v, w, s, r)
					+(yudag(w,t).imag() * yu(p,v).real() + yudag(w,t).real() * yu(p,v).imag()) * WC7R(Cqu8R, v, w, s, r))

					//RGE 2
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CddR
	for (d = 0; d < DWC6R; d++) {
		p = WC6R_indices[d][0];
		r = WC6R_indices[d][1];
		s = WC6R_indices[d][2];
		t = WC6R_indices[d][3];
		f[c] = +(gd(p,r).real() * WC2R(CHdR, s, t) - gd(p,r).imag() * WC2I(CHdI, s, t))
			+(gd(s,t).real() * WC2R(CHdR, p, r) - gd(s,t).imag() * WC2I(CHdI, p, r))
			//RGE 2
			+ TWO_THIRDS * g12 * YhYd *
			(WC2R(CHdR, s, t) * delta[p][r] + WC2R(CHdR, p, r) * delta[s][t])
			+ 6. * g32 * WC6R(CddR, p, t, s, r) - 6. * (g32 / NC - 2. * Yd2 * g12) * WC6R(CddR, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += +(gd(p,v).real() * WC6R(CddR, v, r, s, t) - gd(p,v).imag() * WC6I(CddI, v, r, s, t))
				+(gd(s,v).real() * WC6R(CddR, p, r, v, t) - gd(s,v).imag() * WC6I(CddI, p, r, v, t))
				+(WC6R(CddR, p, v, s, t) * gd(v,r).real() - WC6I(CddI, p, v, s, t) * gd(v,r).imag())
				+(WC6R(CddR, p, r, s, v) * gd(v,t).real() - WC6I(CddI, p, r, s, v) * gd(v,t).imag())
				//RGE 2
				+ TWO_THIRDS * g12 * NC * Yd2 * (
				WC6R(CddR, p, r, w, w) * delta[s][t] + WC6R(CddR, s, t, w, w) * delta[p][r]
				+ WC6R(CddR, w, w, p, r) * delta[s][t] + WC6R(CddR, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * Yd2 * (
				WC6R(CddR, p, w, w, r) * delta[s][t] + WC6R(CddR, s, w, w, t) * delta[p][r]
				+ WC6R(CddR, w, t, s, w) * delta[p][r] + WC6R(CddR, w, r, p, w) * delta[s][t])
				+ ONE_THIRD * g32 * (
				WC6R(CddR, p, w, w, t) * delta[r][s] + WC6R(CddR, s, w, w, r) * delta[p][t]
				+ WC6R(CddR, w, t, p, w) * delta[r][s] + WC6R(CddR, w, r, s, w) * delta[p][t])
				-(ONE_THIRD / NC) * g32 * (
				WC6R(CddR, p, w, w, r) * delta[s][t] + WC6R(CddR, s, w, w, t) * delta[p][r]
				+ WC6R(CddR, w, t, s, w) * delta[p][r] + WC6R(CddR, w, r, p, w) * delta[s][t])

				+ FOUR_THIRDS * g12 * YdYl * (
				WC7R(CldR, w, w, p, r) * delta[s][t] + WC7R(CldR, w, w, s, t) * delta[p][r])
				+ FOUR_THIRDS * g12 * NC * YdYq * (
				WC7R(Cqd1R, w, w, p, r) * delta[s][t] + WC7R(Cqd1R, w, w, s, t) * delta[p][r])
				+ ONE_THIRD * g32 * (
				WC7R(Cqd8R, w, w, s, r) * delta[p][t] + WC7R(Cqd8R, w, w, p, t) * delta[r][s])
				-(ONE_THIRD / NC) * g32 * (
				WC7R(Cqd8R, w, w, p, r) * delta[s][t] + WC7R(Cqd8R, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * YdYe * (
				WC7R(CedR, w, w, p, r) * delta[s][t] + WC7R(CedR, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YuYd * (
				WC7R(Cud1R, w, w, p, r) * delta[s][t] + WC7R(Cud1R, w, w, s, t) * delta[p][r])
				+ ONE_SIXTH * g32 * (
				WC7R(Cud8R, w, w, p, t) * delta[r][s] + WC7R(Cud8R, w, w, s, r) * delta[p][t])
				-(ONE_SIXTH / NC) * g32 * (
				WC7R(Cud8R, w, w, p, r) * delta[s][t] + WC7R(Cud8R, w, w, s, t) * delta[p][r])
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -(
					+(yddag(w,r).real() * yd(p,v).real() - yddag(w,r).imag() * yd(p,v).imag()) * WC7R(Cqd1R, v, w, s, t)
					-(yddag(w,r).imag() * yd(p,v).real() + yddag(w,r).real() * yd(p,v).imag()) * WC7I(Cqd1I, v, w, s, t))
					-(
					+(yddag(w,t).real() * yd(s,v).real() - yddag(w,t).imag() * yd(s,v).imag()) * WC7R(Cqd1R, v, w, p, r)
					-(yddag(w,t).imag() * yd(s,v).real() + yddag(w,t).real() * yd(s,v).imag()) * WC7I(Cqd1I, v, w, p, r))

					+(0.5 / NC)*(
					+(yd(p,v).real() * yddag(w,r).real() - yd(p,v).imag() * yddag(w,r).imag()) * WC7R(Cqd8R, v, w, s, t)
					-(yd(p,v).imag() * yddag(w,r).real() + yd(p,v).real() * yddag(w,r).imag()) * WC7I(Cqd8I, v, w, s, t))
					+(0.5 / NC)*(
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7R(Cqd8R, v, w, p, r)
					-(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7I(Cqd8I, v, w, p, r))

					- 0.5 * (
					+(yddag(w,r).real() * yd(s,v).real() - yddag(w,r).imag() * yd(s,v).imag()) * WC7R(Cqd8R, v, w, p, t)
					-(yddag(w,r).imag() * yd(s,v).real() + yddag(w,r).real() * yd(s,v).imag()) * WC7I(Cqd8I, v, w, p, t))
					- 0.5 * (
					+(yddag(w,t).real() * yd(p,v).real() - yddag(w,t).imag() * yd(p,v).imag()) * WC7R(Cqd8R, v, w, s, r)
					-(yddag(w,t).imag() * yd(p,v).real() + yddag(w,t).real() * yd(p,v).imag()) * WC7I(Cqd8I, v, w, s, r))

					//RGE 2
					;
			}
		}

		f[c] *= loop_factor;
		c++;
	}
	//CddI
	for (d = 0; d < DWC6I; d++) {
		p = WC6I_indices[d][0];
		r = WC6I_indices[d][1];
		s = WC6I_indices[d][2];
		t = WC6I_indices[d][3];
		f[c] = +(gd(p,r).imag() * WC2R(CHdR, s, t) + gd(p,r).real() * WC2I(CHdI, s, t))
			+(gd(s,t).imag() * WC2R(CHdR, p, r) + gd(s,t).real() * WC2I(CHdI, p, r))
			//RGE 2
			+ TWO_THIRDS * g12 * YhYd *
			(WC2I(CHdI, s, t) * delta[p][r] + WC2I(CHdI, p, r) * delta[s][t])
			+ 6. * g32 * WC6I(CddI, p, t, s, r) - 6. * (g32 / NC - 2. * Yd2 * g12) * WC6I(CddI, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += +(gd(p,v).imag() * WC6R(CddR, v, r, s, t) + gd(p,v).real() * WC6I(CddI, v, r, s, t))
				+(gd(s,v).imag() * WC6R(CddR, p, r, v, t) + gd(s,v).real() * WC6I(CddI, p, r, v, t))
				+(WC6R(CddR, p, v, s, t) * gd(v,r).imag() + WC6I(CddI, p, v, s, t) * gd(v,r).real())
				+(WC6R(CddR, p, r, s, v) * gd(v,t).imag() + WC6I(CddI, p, r, s, v) * gd(v,t).real())
				//RGE 2
				+ TWO_THIRDS * g12 * NC * Yd2 * (
				WC6I(CddI, p, r, w, w) * delta[s][t] + WC6I(CddI, s, t, w, w) * delta[p][r]
				+ WC6I(CddI, w, w, p, r) * delta[s][t] + WC6I(CddI, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * Yd2 * (
				WC6I(CddI, p, w, w, r) * delta[s][t] + WC6I(CddI, s, w, w, t) * delta[p][r]
				+ WC6I(CddI, w, t, s, w) * delta[p][r] + WC6I(CddI, w, r, p, w) * delta[s][t])
				+ ONE_THIRD * g32 * (
				WC6I(CddI, p, w, w, t) * delta[r][s] + WC6I(CddI, s, w, w, r) * delta[p][t]
				+ WC6I(CddI, w, t, p, w) * delta[r][s] + WC6I(CddI, w, r, s, w) * delta[p][t])
				-(ONE_THIRD / NC) * g32 * (
				WC6I(CddI, p, w, w, r) * delta[s][t] + WC6I(CddI, s, w, w, t) * delta[p][r]
				+ WC6I(CddI, w, t, s, w) * delta[p][r] + WC6I(CddI, w, r, p, w) * delta[s][t])

				+ FOUR_THIRDS * g12 * YdYl * (
				WC7I(CldI, w, w, p, r) * delta[s][t] + WC7I(CldI, w, w, s, t) * delta[p][r])
				+ FOUR_THIRDS * g12 * NC * YdYq * (
				WC7I(Cqd1I, w, w, p, r) * delta[s][t] + WC7I(Cqd1I, w, w, s, t) * delta[p][r])
				+ ONE_THIRD * g32 * (
				WC7I(Cqd8I, w, w, s, r) * delta[p][t] + WC7I(Cqd8I, w, w, p, t) * delta[r][s])
				-(ONE_THIRD / NC) * g32 * (
				WC7I(Cqd8I, w, w, p, r) * delta[s][t] + WC7I(Cqd8I, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * YdYe * (
				WC7I(CedI, w, w, p, r) * delta[s][t] + WC7I(CedI, w, w, s, t) * delta[p][r])
				+ TWO_THIRDS * g12 * NC * YuYd * (
				WC7I(Cud1I, w, w, p, r) * delta[s][t] + WC7I(Cud1I, w, w, s, t) * delta[p][r])
				+ ONE_SIXTH * g32 * (
				WC7I(Cud8I, w, w, p, t) * delta[r][s] + WC7I(Cud8I, w, w, s, r) * delta[p][t])
				-(ONE_SIXTH / NC) * g32 * (
				WC7I(Cud8I, w, w, p, r) * delta[s][t] + WC7I(Cud8I, w, w, s, t) * delta[p][r])
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -(
					+(yddag(w,r).real() * yd(p,v).real() - yddag(w,r).imag() * yd(p,v).imag()) * WC7I(Cqd1I, v, w, s, t)
					+(yddag(w,r).imag() * yd(p,v).real() + yddag(w,r).real() * yd(p,v).imag()) * WC7R(Cqd1R, v, w, s, t))
					-(
					+(yddag(w,t).real() * yd(s,v).real() - yddag(w,t).imag() * yd(s,v).imag()) * WC7I(Cqd1I, v, w, p, r)
					+(yddag(w,t).imag() * yd(s,v).real() + yddag(w,t).real() * yd(s,v).imag()) * WC7R(Cqd1R, v, w, p, r))

					+(0.5 / NC)*(
					+(yd(p,v).real() * yddag(w,r).real() - yd(p,v).imag() * yddag(w,r).imag()) * WC7I(Cqd8I, v, w, s, t)
					+(yd(p,v).imag() * yddag(w,r).real() + yd(p,v).real() * yddag(w,r).imag()) * WC7R(Cqd8R, v, w, s, t))
					+(0.5 / NC)*(
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7I(Cqd8I, v, w, p, r)
					+(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7R(Cqd8R, v, w, p, r))

					- 0.5 * (
					+(yddag(w,r).real() * yd(s,v).real() - yddag(w,r).imag() * yd(s,v).imag()) * WC7I(Cqd8I, v, w, p, t)
					+(yddag(w,r).imag() * yd(s,v).real() + yddag(w,r).real() * yd(s,v).imag()) * WC7R(Cqd8R, v, w, p, t))
					- 0.5 * (
					+(yddag(w,t).real() * yd(p,v).real() - yddag(w,t).imag() * yd(p,v).imag()) * WC7I(Cqd8I, v, w, s, r)
					+(yddag(w,t).imag() * yd(p,v).real() + yddag(w,t).real() * yd(p,v).imag()) * WC7R(Cqd8R, v, w, s, r))

					//RGE 2
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CeuR
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] = 2. * (ge(p,r).real() * WC2R(CHuR, s, t) - ge(p,r).imag() * WC2I(CHuI, s, t))
			- 2. * (gu(s,t).real() * WC2R(CHeR, p, r) - gu(s,t).imag() * WC2I(CHeI, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYe * WC2R(CHuR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYu * WC2R(CHeR, p, r) * delta[s][t]
			+ 12. * g12 * YuYe * WC7R(CeuR, p, r, s, t); //RGE 3

		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += (ge(p,v).real() * WC7R(CeuR, v, r, s, t) - ge(p,v).imag() * WC7I(CeuI, v, r, s, t))
				+(gu(s,v).real() * WC7R(CeuR, p, r, v, t) - gu(s,v).imag() * WC7I(CeuI, p, r, v, t))
				+(WC7R(CeuR, p, v, s, t) * ge(v,r).real() - WC7I(CeuI, p, v, s, t) * ge(v,r).imag())
				+(WC7R(CeuR, p, r, s, v) * gu(v,t).real() - WC7I(CeuI, p, r, s, v) * gu(v,t).imag())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YuYl * WC7R(CleR, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * YeYl * WC7R(CluR, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * NC * YuYq * WC7R(CqeR, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YeYq * WC7R(Cqu1R, w, w, s, t) * delta[p][r]

				+ FOUR_THIRDS * g12 * NC * YuYe * (
				+WC6R(CuuR, s, t, w, w) * delta[p][r] + WC6R(CuuR, w, w, s, t) * delta[p][r]
				)
				+ FOUR_THIRDS * g12 * YuYe * (
				+WC6R(CuuR, s, w, w, t) * delta[p][r] + WC6R(CuuR, w, t, s, w) * delta[p][r]
				)
				+ FOUR_THIRDS * g12 * YuYe * (
				+WC8R(CeeR, p, r, w, w) * delta[s][t] + WC8R(CeeR, w, w, p, r) * delta[s][t]
				+ WC8R(CeeR, p, w, w, r) * delta[s][t] + WC8R(CeeR, w, r, p, w) * delta[s][t]
				)
				+ FOUR_THIRDS * g12 * NC * Yu2 * WC7R(CeuR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * Ye2 * WC7R(CeuR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYe * WC7R(Cud1R, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7R(CedR, p, r, w, w) * delta[s][t]
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += (
					+(ye(p,v).real() * yu(s,w).real() - ye(p,v).imag() * yu(s,w).imag()) * WC5(Clequ1R, v, r, w, t)
					- (ye(p,v).imag() * yu(s,w).real() + ye(p,v).real() * yu(s,w).imag()) * WC5(Clequ1I, v, r, w, t)
					+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC5(Clequ1R, v, p, w, s)
					+(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC5(Clequ1I, v, p, w, s)
					)
					- 2. * (
					+(ye(p,v).real() * yedag(w,r).real() - ye(p,v).imag() * yedag(w,r).imag()) * WC7R(CluR, v, w, s, t)
					-(ye(p,v).imag() * yedag(w,r).real() + ye(p,v).real() * yedag(w,r).imag()) * WC7I(CluI, v, w, s, t)
					)
					- 12. * (
					+(ye(p,v).real() * yu(s,w).real() - ye(p,v).imag() * yu(s,w).imag()) * WC5(Clequ3R, v, r, w, t)
					- (ye(p,v).imag() * yu(s,w).real() + ye(p,v).real() * yu(s,w).imag()) * WC5(Clequ3I, v, r, w, t)
					+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC5(Clequ3R, v, p, w, s)
					+(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC5(Clequ3I, v, p, w, s)
					)
					- 2. * (
					+(yu(s,v).real() * yudag(w,t).real() - yu(s,v).imag() * yudag(w,t).imag()) * WC7R(CqeR, v, w, p, r)
					- -(yu(s,v).imag() * yudag(w,t).real() + yu(s,v).real() * yudag(w,t).imag()) * WC7I(CqeI, v, w, p, r)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CeuI
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] = 2. * (ge(p,r).imag() * WC2R(CHuR, s, t) + ge(p,r).real() * WC2I(CHuI, s, t))
			- 2. * (gu(s,t).imag() * WC2R(CHeR, p, r) + gu(s,t).real() * WC2I(CHeI, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYe * WC2I(CHuI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYu * WC2I(CHeI, p, r) * delta[s][t]
			+ 12. * g12 * YuYe * WC7I(CeuI, p, r, s, t); //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += (ge(p,v).imag() * WC7R(CeuR, v, r, s, t) + ge(p,v).real() * WC7I(CeuI, v, r, s, t))
				+(gu(s,v).imag() * WC7R(CeuR, p, r, v, t) + gu(s,v).real() * WC7I(CeuI, p, r, v, t))
				+(WC7R(CeuR, p, v, s, t) * ge(v,r).imag() + WC7I(CeuI, p, v, s, t) * ge(v,r).real())
				+(WC7R(CeuR, p, r, s, v) * gu(v,t).imag() + WC7I(CeuI, p, r, s, v) * gu(v,t).real())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YuYl * WC7I(CleI, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * YeYl * WC7I(CluI, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * NC * YuYq * WC7I(CqeI, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YeYq * WC7I(Cqu1I, w, w, s, t) * delta[p][r]

				+ FOUR_THIRDS * g12 * NC * YuYe * (
				+WC6I(CuuI, s, t, w, w) * delta[p][r] + WC6I(CuuI, w, w, s, t) * delta[p][r]
				)
				+ FOUR_THIRDS * g12 * YuYe * (
				+WC6I(CuuI, s, w, w, t) * delta[p][r] + WC6I(CuuI, w, t, s, w) * delta[p][r]
				)
				+ FOUR_THIRDS * g12 * YuYe * (
				+WC8I(CeeI, p, r, w, w) * delta[s][t] + WC8I(CeeI, w, w, p, r) * delta[s][t]
				+ WC8I(CeeI, p, w, w, r) * delta[s][t] + WC8I(CeeI, w, r, p, w) * delta[s][t]
				)
				+ FOUR_THIRDS * g12 * NC * Yu2 * WC7I(CeuI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * Ye2 * WC7I(CeuI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYe * WC7I(Cud1I, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7I(CedI, p, r, w, w) * delta[s][t]
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += (
					+(ye(p,v).real() * yu(s,w).real() - ye(p,v).imag() * yu(s,w).imag()) * WC5(Clequ1I, v, r, w, t)
					+(ye(p,v).imag() * yu(s,w).real() + ye(p,v).real() * yu(s,w).imag()) * WC5(Clequ1R, v, r, w, t)
					-(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC5(Clequ1I, v, p, w, s)
					+(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC5(Clequ1R, v, p, w, s)
					)
					- 2. * (
					+(ye(p,v).real() * yedag(w,r).real() - ye(p,v).imag() * yedag(w,r).imag()) * WC7I(CluI, v, w, s, t)
					+(ye(p,v).imag() * yedag(w,r).real() + ye(p,v).real() * yedag(w,r).imag()) * WC7R(CluR, v, w, s, t)
					)
					- 12. * (
					+(ye(p,v).real() * yu(s,w).real() - ye(p,v).imag() * yu(s,w).imag()) * WC5(Clequ3I, v, r, w, t)
					+ (ye(p,v).imag() * yu(s,w).real() + ye(p,v).real() * yu(s,w).imag()) * WC5(Clequ3R, v, r, w, t)
					- (yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC5(Clequ3I, v, p, w, s)
					+ (yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC5(Clequ3R, v, p, w, s)
					)
					- 2. * (
					+(yu(s,v).real() * yudag(w,t).real() - yu(s,v).imag() * yudag(w,t).imag()) * WC7I(CqeI, v, w, p, r)
					+(yu(s,v).imag() * yudag(w,t).real() + yu(s,v).real() * yudag(w,t).imag()) * WC7R(CqeR, v, w, p, r)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CedR
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];

		f[c] = 2. * (ge(p,r).real() * WC2R(CHdR, s, t) - ge(p,r).imag() * WC2I(CHdI, s, t))
			+ 2. * (gd(s,t).real() * WC2R(CHeR, p, r) - gd(s,t).imag() * WC2I(CHeI, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYe * WC2R(CHdR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYd * WC2R(CHeR, p, r) * delta[s][t]
			+ 12. * YdYe * g12 * WC7R(CedR, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] += (ge(p,v).real() * WC7R(CedR, v, r, s, t) - ge(p,v).imag() * WC7I(CedI, v, r, s, t))
				+(gd(s,v).real() * WC7R(CedR, p, r, v, t) - gd(s,v).imag() * WC7I(CedI, p, r, v, t))
				+(WC7R(CedR, p, v, s, t) * ge(v,r).real() - WC7I(CedI, p, v, s, t) * ge(v,r).imag())
				+(WC7R(CedR, p, r, s, v) * gd(v,t).real() - WC7I(CedI, p, r, s, v) * gd(v,t).imag())
				//RGE 2
				+ FOUR_THIRDS * g12 * NC * YdYe * (
				+WC6R(CddR, s, t, w, w) * delta[p][r] + WC6R(CddR, w, w, s, t) * delta[p][r]
				)
				+ FOUR_THIRDS * g12 * YdYe * (
				+WC6R(CddR, s, w, w, t) * delta[p][r] + WC6R(CddR, w, t, s, w) * delta[p][r]
				)
				+ FOUR_THIRDS * g12 * YdYe * (
				+WC8R(CeeR, p, r, w, w) * delta[s][t] + WC8R(CeeR, w, w, p, r) * delta[s][t]
				+ WC8R(CeeR, p, w, w, r) * delta[s][t] + WC8R(CeeR, w, r, p, w) * delta[s][t]
				)
				+ EIGHT_THIRDS * g12 * YeYl * WC7R(CldR, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * YdYl * WC7R(CleR, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YeYq * WC7R(Cqd1R, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * NC * YdYq * WC7R(CqeR, w, w, p, r) * delta[s][t]

				+ FOUR_THIRDS * g12 * NC * YuYd * WC7R(CeuR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYe * WC7R(Cud1R, w, w, s, t) * delta[p][r]

				+ FOUR_THIRDS * g12 * NC * Yd2 * WC7R(CedR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * Ye2 * WC7R(CedR, w, w, s, t) * delta[p][r]
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -2. * (
					+(ye(p,v).real() * yedag(w,r).real() - ye(p,v).imag() * yedag(w,r).imag()) * WC7R(CldR, v, w, s, t)
					-(ye(p,v).imag() * yedag(w,r).real() + ye(p,v).real() * yedag(w,r).imag()) * WC7I(CldI, v, w, s, t)
					)
					- 2. * (
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7R(CqeR, v, w, p, r)
					-(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7I(CqeI, v, w, p, r)
					)
					+(
					+(ye(p,v).real() * yddag(w,t).real() - ye(p,v).imag() * yddag(w,t).imag()) * WC5(CledqR, v, r, s, w)
					-(ye(p,v).imag() * yddag(w,t).real() + ye(p,v).real() * yddag(w,t).imag()) * WC5(CledqI, v, r, s, w)
					)
					+(
					+(yedag(v,r).real() * yd(s,w).real() - yedag(v,r).imag() * yd(s,w).imag()) * WC5(CledqR, v, p, t, w)
					+(yedag(v,r).imag() * yd(s,w).real() + yedag(v,r).real() * yd(s,w).imag()) * WC5(CledqI, v, p, t, w)
					)
					//RGE 2
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CedI
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];

		f[c] = 2. * (ge(p,r).imag() * WC2R(CHdR, s, t) + ge(p,r).real() * WC2I(CHdI, s, t))
			+ 2. * (gd(s,t).imag() * WC2R(CHeR, p, r) + gd(s,t).real() * WC2I(CHeI, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYe * WC2I(CHdI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYd * WC2I(CHeI, p, r) * delta[s][t]
			+ 12. * YdYe * g12 * WC7I(CedI, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index

			f[c] += (ge(p,v).imag() * WC7R(CedR, v, r, s, t) + ge(p,v).real() * WC7I(CedI, v, r, s, t))
				+(gd(s,v).imag() * WC7R(CedR, p, r, v, t) + gd(s,v).real() * WC7I(CedI, p, r, v, t))
				+(WC7R(CedR, p, v, s, t) * ge(v,r).imag() + WC7I(CedI, p, v, s, t) * ge(v,r).real())
				+(WC7R(CedR, p, r, s, v) * gd(v,t).imag() + WC7I(CedI, p, r, s, v) * gd(v,t).real())
				//RGE 2
				+ FOUR_THIRDS * g12 * NC * YdYe * (
				+WC6I(CddI, s, t, w, w) * delta[p][r] + WC6I(CddI, w, w, s, t) * delta[p][r]
				)
				+ FOUR_THIRDS * g12 * YdYe * (
				+WC6I(CddI, s, w, w, t) * delta[p][r] + WC6I(CddI, w, t, s, w) * delta[p][r]
				)
				+ FOUR_THIRDS * g12 * YdYe * (
				+WC8I(CeeI, p, r, w, w) * delta[s][t] + WC8I(CeeI, w, w, p, r) * delta[s][t]
				+ WC8I(CeeI, p, w, w, r) * delta[s][t] + WC8I(CeeI, w, r, p, w) * delta[s][t]
				)
				+ EIGHT_THIRDS * g12 * YeYl * WC7I(CldI, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * YdYl * WC7I(CleI, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YeYq * WC7I(Cqd1I, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * NC * YdYq * WC7I(CqeI, w, w, p, r) * delta[s][t]

				+ FOUR_THIRDS * g12 * NC * YuYd * WC7I(CeuI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYe * WC7I(Cud1I, w, w, s, t) * delta[p][r]

				+ FOUR_THIRDS * g12 * NC * Yd2 * WC7I(CedI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * Ye2 * WC7I(CedI, w, w, s, t) * delta[p][r]
				; //RGE 3
			;
			for (w = 0; w < NG; w++) {
				f[c] += -2. * (
					+(ye(p,v).real() * yedag(w,r).real() - ye(p,v).imag() * yedag(w,r).imag()) * WC7I(CldI, v, w, s, t)
					+(ye(p,v).imag() * yedag(w,r).real() + ye(p,v).real() * yedag(w,r).imag()) * WC7R(CldR, v, w, s, t)
					)
					- 2. * (
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7I(CqeI, v, w, p, r)
					+(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7R(CqeR, v, w, p, r)
					)
					+(
					+(ye(p,v).real() * yddag(w,t).real() - ye(p,v).imag() * yddag(w,t).imag()) * WC5(CledqI, v, r, s, w)
					+(ye(p,v).imag() * yddag(w,t).real() + ye(p,v).real() * yddag(w,t).imag()) * WC5(CledqR, v, r, s, w)
					)
					+(
					-(yedag(v,r).real() * yd(s,w).real() - yedag(v,r).imag() * yd(s,w).imag()) * WC5(CledqI, v, p, t, w)
					+(yedag(v,r).imag() * yd(s,w).real() + yedag(v,r).real() * yd(s,w).imag()) * WC5(CledqR, v, p, t, w)
					)
					//RGE 2
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//Cud1R
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			-2. * (gu(p,r).real() * WC2R(CHdR, s, t) - gu(p,r).imag() * WC2I(CHdI, s, t))
			+ 2. * (gd(s,t).real() * WC2R(CHuR, p, r) - gd(s,t).imag() * WC2I(CHuI, p, r))
			+(2. / NC)*
			(ydyud(s,r).real() * WC1(CHudR, p, t) - ydyud(s,r).imag() * WC1(CHudI, p, t))
			+(2. / NC)*
			(yuydd(p,t).real() * WC1(CHudR, r, s) + yuydd(p,t).imag() * WC1(CHudI, r, s))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYu * WC2R(CHdR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYd * WC2R(CHuR, p, r) * delta[s][t]
			+ 12. * YuYd * g12 * WC7R(Cud1R, p, r, s, t)
			+ 3. * ((NC2 - 1) / NC2) * g32 * WC7R(Cud8R, p, r, s, t)
			//RGE 3
			;
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gu(p,v).real() * WC7R(Cud1R, v, r, s, t) - gu(p,v).imag() * WC7I(Cud1I, v, r, s, t))
				+(gd(s,v).real() * WC7R(Cud1R, p, r, v, t) - gd(s,v).imag() * WC7I(Cud1I, p, r, v, t))
				+(WC7R(Cud1R, p, v, s, t) * gu(v,r).real() - WC7I(Cud1I, p, v, s, t) * gu(v,r).imag())
				+(WC7R(Cud1R, p, r, s, v) * gd(v,t).real() - WC7I(Cud1I, p, r, s, v) * gd(v,t).imag())
				//RGE 2
				+ FOUR_THIRDS * g12 * NC * YuYd * delta[s][t]*(
				WC6R(CuuR, p, r, w, w) + WC6R(CuuR, w, w, p, r))
				+ FOUR_THIRDS * g12 * YuYd * delta[s][t]*(
				WC6R(CuuR, p, w, w, r) + WC6R(CuuR, w, r, p, w))
				+ FOUR_THIRDS * g12 * NC * YuYd * delta[p][r]*(
				WC6R(CddR, s, t, w, w) + WC6R(CddR, w, w, s, t))
				+ FOUR_THIRDS * g12 * YuYd * delta[p][r]*(
				WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w))
				+ EIGHT_THIRDS * g12 * NC * YdYq * WC7R(Cqu1R, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YuYq * WC7R(Cqd1R, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * YdYl * WC7R(CluR, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * YuYl * WC7R(CldR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * Yd2 * WC7R(Cud1R, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * Yu2 * WC7R(Cud1R, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * YdYe * WC7R(CeuR, w, w, p, r) * delta[s][t]
				+ FOUR_THIRDS * g12 * YuYe * WC7R(CedR, w, w, s, t) * delta[p][r]
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += (1. / NC)*(
					+(yd(s,v).real() * yu(p,w).real() - yd(s,v).imag() * yu(p,w).imag()) * WC5(Cquqd1R, v, r, w, t)
					-(yd(s,v).imag() * yu(p,w).real() + yd(s,v).real() * yu(p,w).imag()) * WC5(Cquqd1I, v, r, w, t)
					+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC5(Cquqd1R, v, p, w, s)
					+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC5(Cquqd1I, v, p, w, s)
					)
					-(
					+(yd(s,w).real() * yu(p,v).real() - yd(s,w).imag() * yu(p,v).imag()) * WC5(Cquqd1R, v, r, w, t)
					-(yd(s,w).imag() * yu(p,v).real() + yd(s,w).real() * yu(p,v).imag()) * WC5(Cquqd1I, v, r, w, t)
					)
					-(
					+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC5(Cquqd1R, v, p, w, s)
					+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC5(Cquqd1I, v, p, w, s)
					)
					+ 0.5 * ((NC2 - 1.) / NC2)*(
					+(yd(s,v).real() * yu(p,w).real() - yd(s,v).imag() * yu(p,w).imag()) * WC5(Cquqd8R, v, r, w, t)
					-(yd(s,v).imag() * yu(p,w).real() + yd(s,v).real() * yu(p,w).imag()) * WC5(Cquqd8I, v, r, w, t)
					+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC5(Cquqd8R, v, p, w, s)
					+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC5(Cquqd8I, v, p, w, s)
					)
					- 2. * (
					+(yu(p,v).real() * yudag(w,r).real() - yu(p,v).imag() * yudag(w,r).imag()) * WC7R(Cqd1R, v, w, s, t)
					-(yu(p,v).imag() * yudag(w,r).real() + yu(p,v).real() * yudag(w,r).imag()) * WC7I(Cqd1I, v, w, s, t)
					)
					- 2. * (
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7R(Cqu1R, v, w, p, r)
					-(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7I(Cqu1I, v, w, p, r)
					)
					//RGE 2 
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//Cud1I
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			-2. * (gu(p,r).imag() * WC2R(CHdR, s, t) + gu(p,r).real() * WC2I(CHdI, s, t))
			+ 2. * (gd(s,t).imag() * WC2R(CHuR, p, r) + gd(s,t).real() * WC2I(CHuI, p, r))
			+(2. / NC)*
			(ydyud(s,r).imag() * WC1(CHudR, p, t) + ydyud(s,r).real() * WC1(CHudI, p, t))
			+(2. / NC)*
			(yuydd(p,t).imag() * WC1(CHudR, r, s) - yuydd(p,t).real() * WC1(CHudI, r, s))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYu * WC2I(CHdI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYd * WC2I(CHuI, p, r) * delta[s][t]
			+ 12. * YuYd * g12 * WC7I(Cud1I, p, r, s, t)
			+ 3. * ((NC2 - 1) / NC2) * g32 * WC7I(Cud8I, p, r, s, t)
			//RGE 3
			;
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gu(p,v).imag() * WC7R(Cud1R, v, r, s, t) + gu(p,v).real() * WC7I(Cud1I, v, r, s, t))
				+(gd(s,v).imag() * WC7R(Cud1R, p, r, v, t) + gd(s,v).real() * WC7I(Cud1I, p, r, v, t))
				+(WC7R(Cud1R, p, v, s, t) * gu(v,r).imag() + WC7I(Cud1I, p, v, s, t) * gu(v,r).real())
				+(WC7R(Cud1R, p, r, s, v) * gd(v,t).imag() + WC7I(Cud1I, p, r, s, v) * gd(v,t).real())
				//RGE 2
				+ FOUR_THIRDS * g12 * NC * YuYd * delta[s][t]*(
				WC6I(CuuI, p, r, w, w) + WC6I(CuuI, w, w, p, r))
				+ FOUR_THIRDS * g12 * YuYd * delta[s][t]*(
				WC6I(CuuI, p, w, w, r) + WC6I(CuuI, w, r, p, w))
				+ FOUR_THIRDS * g12 * NC * YuYd * delta[p][r]*(
				WC6I(CddI, s, t, w, w) + WC6I(CddI, w, w, s, t))
				+ FOUR_THIRDS * g12 * YuYd * delta[p][r]*(
				WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w))
				+ EIGHT_THIRDS * g12 * NC * YdYq * WC7I(Cqu1I, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YuYq * WC7I(Cqd1I, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * YdYl * WC7I(CluI, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * YuYl * WC7I(CldI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * Yd2 * WC7I(Cud1I, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * Yu2 * WC7I(Cud1I, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * YdYe * WC7I(CeuI, w, w, p, r) * delta[s][t]
				+ FOUR_THIRDS * g12 * YuYe * WC7I(CedI, w, w, s, t) * delta[p][r]
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					+(1. / NC)*(
					+(yd(s,v).real() * yu(p,w).real() - yd(s,v).imag() * yu(p,w).imag()) * WC5(Cquqd1I, v, r, w, t)
					+(yd(s,v).imag() * yu(p,w).real() + yd(s,v).real() * yu(p,w).imag()) * WC5(Cquqd1R, v, r, w, t)
					-(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC5(Cquqd1I, v, p, w, s)
					+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC5(Cquqd1R, v, p, w, s)
					)
					-(
					+(yd(s,w).real() * yu(p,v).real() - yd(s,w).imag() * yu(p,v).imag()) * WC5(Cquqd1I, v, r, w, t)
					+(yd(s,w).imag() * yu(p,v).real() + yd(s,w).real() * yu(p,v).imag()) * WC5(Cquqd1R, v, r, w, t)
					)
					-(
					-(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC5(Cquqd1I, v, p, w, s)
					+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC5(Cquqd1R, v, p, w, s)
					)
					+ 0.5 * ((NC2 - 1.) / NC2)*(
					+(yd(s,v).real() * yu(p,w).real() - yd(s,v).imag() * yu(p,w).imag()) * WC5(Cquqd8I, v, r, w, t)
					+(yd(s,v).imag() * yu(p,w).real() + yd(s,v).real() * yu(p,w).imag()) * WC5(Cquqd8R, v, r, w, t)
					-(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC5(Cquqd8I, v, p, w, s)
					+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC5(Cquqd8R, v, p, w, s)
					)
					- 2. * (
					+(yu(p,v).real() * yudag(w,r).real() - yu(p,v).imag() * yudag(w,r).imag()) * WC7I(Cqd1I, v, w, s, t)
					+(yu(p,v).imag() * yudag(w,r).real() + yu(p,v).real() * yudag(w,r).imag()) * WC7R(Cqd1R, v, w, s, t)
					)
					- 2. * (
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7I(Cqu1I, v, w, p, r)
					+(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7R(Cqu1R, v, w, p, r)
					)
					//RGE 2 
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//Cud8R
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			+4. * (ydyud(s,r).real() * WC1(CHudR, p, t) - ydyud(s,r).imag() * WC1(CHudI, p, t))
			+ 4. * (yuydd(p,t).real() * WC1(CHudR, r, s) + yuydd(p,t).imag() * WC1(CHudI, r, s))
			//RGE 2
			+ 12. * (YuYd * g12 - g32 / NC) * WC7R(Cud8R, p, r, s, t)
			+ 12. * g32 * WC7R(Cud1R, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gu(p,v).real() * WC7R(Cud8R, v, r, s, t) - gu(p,v).imag() * WC7I(Cud8I, v, r, s, t))
				+(gd(s,v).real() * WC7R(Cud8R, p, r, v, t) - gd(s,v).imag() * WC7I(Cud8I, p, r, v, t))
				+(WC7R(Cud8R, p, v, s, t) * gu(v,r).real() - WC7I(Cud8I, p, v, s, t) * gu(v,r).imag())
				+(WC7R(Cud8R, p, r, s, v) * gd(v,t).real() - WC7I(Cud8I, p, r, s, v) * gd(v,t).imag())
				//RGE 2
				+ FOUR_THIRDS * g32 * delta[s][t]*(
				+WC6R(CuuR, p, w, w, r) + WC6R(CuuR, w, r, p, w)
				)
				+ FOUR_THIRDS * g32 * delta[p][r]*(
				+WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w)
				)
				+ FOUR_THIRDS * g32 * WC7R(Cqu8R, w, w, p, r) * delta[s][t]
				+ FOUR_THIRDS * g32 * WC7R(Cqd8R, w, w, s, t) * delta[p][r]
				+ TWO_THIRDS * g32 * (
				+WC7R(Cud8R, p, r, w, w) * delta[s][t] + WC7R(Cud8R, w, w, s, t) * delta[p][r]
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += +2. * (
					+(yd(s,v).real() * yu(p,w).real() - yd(s,v).imag() * yu(p,w).imag()) * WC5(Cquqd1R, v, r, w, t)
					-(yd(s,v).imag() * yu(p,w).real() + yd(s,v).real() * yu(p,w).imag()) * WC5(Cquqd1I, v, r, w, t)
					+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC5(Cquqd1R, v, p, w, s)
					+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC5(Cquqd1I, v, p, w, s)
					)
					- 2. * (
					+(yu(p,v).real() * yudag(w,r).real() - yu(p,v).imag() * yudag(w,r).imag()) * WC7R(Cqd8R, v, w, s, t)
					-(yu(p,v).imag() * yudag(w,r).real() + yu(p,v).real() * yudag(w,r).imag()) * WC7I(Cqd8I, v, w, s, t)
					)
					- 2. * (
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7R(Cqu8R, v, w, p, r)
					-(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7I(Cqu8I, v, w, p, r)
					)
					-(1. / NC)*(
					+(yd(s,v).real() * yu(p,w).real() - yd(s,v).imag() * yu(p,w).imag()) * WC5(Cquqd8R, v, r, w, t)
					-(yd(s,v).imag() * yu(p,w).real() + yd(s,v).real() * yu(p,w).imag()) * WC5(Cquqd8I, v, r, w, t)
					+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC5(Cquqd8R, v, p, w, s)
					+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC5(Cquqd8I, v, p, w, s)
					)
					-(
					+(yd(s,w).real() * yu(p,v).real() - yd(s,w).imag() * yu(p,v).imag()) * WC5(Cquqd8R, v, r, w, t)
					- (yd(s,w).imag() * yu(p,v).real() + yd(s,w).real() * yu(p,v).imag()) * WC5(Cquqd8I, v, r, w, t)
					+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC5(Cquqd8R, v, p, w, s)
					+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC5(Cquqd8I, v, p, w, s)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//Cud8I
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			+4. * (ydyud(s,r).imag() * WC1(CHudR, p, t) + ydyud(s,r).real() * WC1(CHudI, p, t))
			+ 4. * (yuydd(p,t).imag() * WC1(CHudR, r, s) - yuydd(p,t).real() * WC1(CHudI, r, s))
			//RGE 2
			+ 12. * (YuYd * g12 - g32 / NC) * WC7I(Cud8I, p, r, s, t)
			+ 12. * g32 * WC7I(Cud1I, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gu(p,v).imag() * WC7R(Cud8R, v, r, s, t) + gu(p,v).real() * WC7I(Cud8I, v, r, s, t))
				+(gd(s,v).imag() * WC7R(Cud8R, p, r, v, t) + gd(s,v).real() * WC7I(Cud8I, p, r, v, t))
				+(WC7R(Cud8R, p, v, s, t) * gu(v,r).imag() + WC7I(Cud8I, p, v, s, t) * gu(v,r).real())
				+(WC7R(Cud8R, p, r, s, v) * gd(v,t).imag() + WC7I(Cud8I, p, r, s, v) * gd(v,t).real())
				//RGE 2
				+ FOUR_THIRDS * g32 * delta[s][t]*(
				+WC6I(CuuI, p, w, w, r) + WC6I(CuuI, w, r, p, w)
				)
				+ FOUR_THIRDS * g32 * delta[p][r]*(
				+WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w)
				)
				+ FOUR_THIRDS * g32 * WC7I(Cqu8I, w, w, p, r) * delta[s][t]
				+ FOUR_THIRDS * g32 * WC7I(Cqd8I, w, w, s, t) * delta[p][r]
				+ TWO_THIRDS * g32 * (
				+WC7I(Cud8I, p, r, w, w) * delta[s][t] + WC7I(Cud8I, w, w, s, t) * delta[p][r]
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += +2. * (
					+(yd(s,v).real() * yu(p,w).real() - yd(s,v).imag() * yu(p,w).imag()) * WC5(Cquqd1I, v, r, w, t)
					+(yd(s,v).imag() * yu(p,w).real() + yd(s,v).real() * yu(p,w).imag()) * WC5(Cquqd1R, v, r, w, t)
					-(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC5(Cquqd1I, v, p, w, s)
					+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC5(Cquqd1R, v, p, w, s)
					)
					- 2. * (
					+(yu(p,v).real() * yudag(w,r).real() - yu(p,v).imag() * yudag(w,r).imag()) * WC7I(Cqd8I, v, w, s, t)
					+(yu(p,v).imag() * yudag(w,r).real() + yu(p,v).real() * yudag(w,r).imag()) * WC7R(Cqd8R, v, w, s, t)
					)
					- 2. * (
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7I(Cqu8I, v, w, p, r)
					+(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7R(Cqu8R, v, w, p, r)
					)
					-(1. / NC)*(
					+(yd(s,v).real() * yu(p,w).real() - yd(s,v).imag() * yu(p,w).imag()) * WC5(Cquqd8I, v, r, w, t)
					+(yd(s,v).imag() * yu(p,w).real() + yd(s,v).real() * yu(p,w).imag()) * WC5(Cquqd8R, v, r, w, t)
					-(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC5(Cquqd8I, v, p, w, s)
					+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC5(Cquqd8R, v, p, w, s)
					)
					-(
					+(yd(s,w).real() * yu(p,v).real() - yd(s,w).imag() * yu(p,v).imag()) * WC5(Cquqd8I, v, r, w, t)
					+ (yd(s,w).imag() * yu(p,v).real() + yd(s,w).real() * yu(p,v).imag()) * WC5(Cquqd8R, v, r, w, t)
					-(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC5(Cquqd8I, v, p, w, s)
					+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC5(Cquqd8R, v, p, w, s)
					)
					; //RGE 2

			}
		}
		f[c] *= loop_factor;
		c++;
	}



	//----------------------RGE SMEFT class 8_LLRR----------------------

	//CleR
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			+(ye(s,r).real() * xie(p,t).real() - ye(s,r).imag() * xie(p,t).imag())
			+(yedag(p,t).real() * xie(r,s).real() + yedag(p,t).imag() * xie(r,s).imag())
			-(yedye(p,r).real() * WC2R(CHeR, s, t) - yedye(p,r).imag() * WC2I(CHeI, s, t))
			+ 2. * (ge(s,t).real() * WC2R(CHl1R, p, r) - ge(s,t).imag() * WC2I(CHl1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYl * WC2R(CHeR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYe * WC2R(CHl1R, p, r) * delta[s][t]
			- 12. * YeYl * g12 * WC7R(CleR, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).real() * WC7R(CleR, v, r, s, t) - gl(p,v).imag() * WC7I(CleI, v, r, s, t))
				+(ge(s,v).real() * WC7R(CleR, p, r, v, t) - ge(s,v).imag() * WC7I(CleI, p, r, v, t))
				+(WC7R(CleR, p, v, s, t) * gl(v,r).real() - WC7I(CleI, p, v, s, t) * gl(v,r).imag())
				+(WC7R(CleR, p, r, s, v) * ge(v,t).real() - WC7I(CleI, p, r, s, v) * ge(v,t).imag())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YeYl * delta[s][t]*(
				+WC6R(CllR, p, r, w, w) + WC6R(CllR, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YeYl * delta[s][t]*(
				+WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * NC * YeYq * WC7R(Clq1R, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YlYq * WC7R(CqeR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * Ye2 * WC7R(CleR, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * Yl2 * WC7R(CleR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYe * WC7R(CluR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YdYe * WC7R(CldR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYl * WC7R(CeuR, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYl * WC7R(CedR, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYl * delta[p][r]* (
				+WC8R(CeeR, s, t, w, w) + WC8R(CeeR, s, w, w, t)
				+ WC8R(CeeR, w, t, s, w) + WC8R(CeeR, w, w, s, t)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					-(
					+(yedag(p,v).real() * ye(w,r).real() - yedag(p,v).imag() * ye(w,r).imag()) * WC8R(CeeR, v, t, s, w)
					-(yedag(p,v).imag() * ye(w,r).real() + yedag(p,v).real() * ye(w,r).imag()) * WC8I(CeeI, v, t, s, w)
					)
					-(
					+(yedag(p,w).real() * ye(v,r).real() - yedag(p,w).imag() * ye(v,r).imag()) * WC8R(CeeR, w, t, s, v)
					-(yedag(p,w).imag() * ye(v,r).real() + yedag(p,w).real() * ye(v,r).imag()) * WC8I(CeeI, w, t, s, v)
					)
					- 2. * (
					+(yedag(p,v).real() * ye(w,r).real() - yedag(p,v).imag() * ye(w,r).imag()) * WC8R(CeeR, v, w, s, t)
					-(yedag(p,v).imag() * ye(w,r).real() + yedag(p,v).real() * ye(w,r).imag()) * WC8I(CeeI, v, w, s, t)
					)
					+(
					+(yedag(p,w).real() * ye(s,v).real() - yedag(p,w).imag() * ye(s,v).imag()) * WC7R(CleR, v, r, w, t)
					-(yedag(p,w).imag() * ye(s,v).real() + yedag(p,w).real() * ye(s,v).imag()) * WC7I(CleI, v, r, w, t)
					)
					-(
					+(yedag(w,t).real() * ye(s,v).real() - yedag(w,t).imag() * ye(s,v).imag()) * WC6R(CllR, p, w, v, r)
					-(yedag(w,t).imag() * ye(s,v).real() + yedag(w,t).real() * ye(s,v).imag()) * WC6I(CllI, p, w, v, r)
					)
					-(
					+(yedag(v,t).real() * ye(s,w).real() - yedag(v,t).imag() * ye(s,w).imag()) * WC6R(CllR, p, v, w, r)
					-(yedag(v,t).imag() * ye(s,w).real() + yedag(v,t).real() * ye(s,w).imag()) * WC6I(CllI, p, v, w, r)
					)
					- 4. * (
					+(yedag(w,t).real() * ye(s,v).real() - yedag(w,t).imag() * ye(s,v).imag()) * WC6R(CllR, p, r, v, w)
					-(yedag(w,t).imag() * ye(s,v).real() + yedag(w,t).real() * ye(s,v).imag()) * WC6I(CllI, p, r, v, w)
					)
					+(
					+(yedag(v,t).real() * ye(w,r).real() - yedag(v,t).imag() * ye(w,r).imag()) * WC7R(CleR, p, v, s, w)
					-(yedag(v,t).imag() * ye(w,r).real() + yedag(v,t).real() * ye(w,r).imag()) * WC7I(CleI, p, v, s, w)
					)
					//RGE 2
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CleI
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			+(ye(s,r).imag() * xie(p,t).real() + ye(s,r).real() * xie(p,t).imag())
			+(yedag(p,t).imag() * xie(r,s).real() - yedag(p,t).real() * xie(r,s).imag())
			-(yedye(p,r).imag() * WC2R(CHeR, s, t) + yedye(p,r).real() * WC2I(CHeI, s, t))
			+ 2. * (ge(s,t).imag() * WC2R(CHl1R, p, r) + ge(s,t).real() * WC2I(CHl1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYl * WC2I(CHeI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYe * WC2I(CHl1I, p, r) * delta[s][t]
			- 12. * YeYl * g12 * WC7I(CleI, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).imag() * WC7R(CleR, v, r, s, t) + gl(p,v).real() * WC7I(CleI, v, r, s, t))
				+(ge(s,v).imag() * WC7R(CleR, p, r, v, t) + ge(s,v).real() * WC7I(CleI, p, r, v, t))
				+(WC7R(CleR, p, v, s, t) * gl(v,r).imag() + WC7I(CleI, p, v, s, t) * gl(v,r).real())
				+(WC7R(CleR, p, r, s, v) * ge(v,t).imag() + WC7I(CleI, p, r, s, v) * ge(v,t).real())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YeYl * delta[s][t]*(
				+WC6I(CllI, p, r, w, w) + WC6I(CllI, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YeYl * delta[s][t]*(
				+WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * NC * YeYq * WC7I(Clq1I, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YlYq * WC7I(CqeI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * Ye2 * WC7I(CleI, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * Yl2 * WC7I(CleI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYe * WC7I(CluI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YdYe * WC7I(CldI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYl * WC7I(CeuI, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYl * WC7I(CedI, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYl * delta[p][r]* (
				+WC8I(CeeI, s, t, w, w) + WC8I(CeeI, s, w, w, t)
				+ WC8I(CeeI, w, t, s, w) + WC8I(CeeI, w, w, s, t)
				)
				; //RGE 3

			for (w = 0; w < NG; w++) {
				f[c] += -(
					+(yedag(p,v).real() * ye(w,r).real() - yedag(p,v).imag() * ye(w,r).imag()) * WC8I(CeeI, v, t, s, w)
					+(yedag(p,v).imag() * ye(w,r).real() + yedag(p,v).real() * ye(w,r).imag()) * WC8R(CeeR, v, t, s, w)
					)
					-(
					+(yedag(p,w).real() * ye(v,r).real() - yedag(p,w).imag() * ye(v,r).imag()) * WC8I(CeeI, w, t, s, v)
					+(yedag(p,w).imag() * ye(v,r).real() + yedag(p,w).real() * ye(v,r).imag()) * WC8R(CeeR, w, t, s, v)
					)
					- 2. * (
					+(yedag(p,v).real() * ye(w,r).real() - yedag(p,v).imag() * ye(w,r).imag()) * WC8I(CeeI, v, w, s, t)
					+(yedag(p,v).imag() * ye(w,r).real() + yedag(p,v).real() * ye(w,r).imag()) * WC8R(CeeR, v, w, s, t)
					)
					+(
					+(yedag(p,w).real() * ye(s,v).real() - yedag(p,w).imag() * ye(s,v).imag()) * WC7I(CleI, v, r, w, t)
					+(yedag(p,w).imag() * ye(s,v).real() + yedag(p,w).real() * ye(s,v).imag()) * WC7R(CleR, v, r, w, t)
					)
					-(
					+(yedag(w,t).real() * ye(s,v).real() - yedag(w,t).imag() * ye(s,v).imag()) * WC6I(CllI, p, w, v, r)
					+(yedag(w,t).imag() * ye(s,v).real() + yedag(w,t).real() * ye(s,v).imag()) * WC6R(CllR, p, w, v, r)
					)
					-(
					+(yedag(v,t).real() * ye(s,w).real() - yedag(v,t).imag() * ye(s,w).imag()) * WC6I(CllI, p, v, w, r)
					+(yedag(v,t).imag() * ye(s,w).real() + yedag(v,t).real() * ye(s,w).imag()) * WC6R(CllR, p, v, w, r)
					)
					- 4. * (
					+(yedag(w,t).real() * ye(s,v).real() - yedag(w,t).imag() * ye(s,v).imag()) * WC6I(CllI, p, r, v, w)
					+(yedag(w,t).imag() * ye(s,v).real() + yedag(w,t).real() * ye(s,v).imag()) * WC6R(CllR, p, r, v, w)
					)
					+(
					+(yedag(v,t).real() * ye(w,r).real() - yedag(v,t).imag() * ye(w,r).imag()) * WC7I(CleI, p, v, s, w)
					+(yedag(v,t).imag() * ye(w,r).real() + yedag(v,t).real() * ye(w,r).imag()) * WC7R(CleR, p, v, s, w)
					)
					//RGE 2
					;
			}
		}
		f[c] *= loop_factor;
		c++;
	}




	//CluR
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			-(yedye(p,r).real() * WC2R(CHuR, s, t) - yedye(p,r).imag() * WC2I(CHuI, s, t))
			- 2. * (gu(s,t).real() * WC2R(CHl1R, p, r) - gu(s,t).imag() * WC2I(CHl1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYl * WC2R(CHuR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYu * WC2R(CHl1R, p, r) * delta[s][t]
			- 12. * g12 * YuYl * WC7R(CluR, p, r, s, t)
			; //RGE 3

		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).real() * WC7R(CluR, v, r, s, t) - gl(p,v).imag() * WC7I(CluI, v, r, s, t))
				+(gu(s,v).real() * WC7R(CluR, p, r, v, t) - gu(s,v).imag() * WC7I(CluI, p, r, v, t))
				+(WC7R(CluR, p, v, s, t) * gl(v,r).real() - WC7I(CluI, p, v, s, t) * gl(v,r).imag())
				+(WC7R(CluR, p, r, s, v) * gu(v,t).real() - WC7I(CluI, p, r, s, v) * gu(v,t).imag())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YuYl * delta[s][t]*(
				+WC6R(CllR, p, r, w, w) + WC6R(CllR, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YuYl * delta[s][t]*(
				+WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * NC * YuYq * WC7R(Clq1R, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YlYq * WC7R(Cqu1R, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * Yu2 * WC7R(CluR, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * Yl2 * WC7R(CluR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7R(CldR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * YuYe * WC7R(CleR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YdYl * WC7R(Cud1R, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYl * WC7R(CeuR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYl * delta[p][r]* (
				+WC6R(CuuR, s, t, w, w) + WC6R(CuuR, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YuYl * delta[p][r]* (
				+WC6R(CuuR, s, w, w, t) + WC6R(CuuR, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					-0.5 * (
					+(ye(v,r).real() * yu(s,w).real() - ye(v,r).imag() * yu(s,w).imag()) * WC5(Clequ1R, p, v, w, t)
					-(ye(v,r).imag() * yu(s,w).real() + ye(v,r).real() * yu(s,w).imag()) * WC5(Clequ1I, p, v, w, t)
					+(yedag(p,v).real() * yudag(w,t).real() - yedag(p,v).imag() * yudag(w,t).imag()) * WC5(Clequ1R, r, v, w, s)
					+(yedag(p,v).imag() * yudag(w,t).real() + yedag(p,v).real() * yudag(w,t).imag()) * WC5(Clequ1I, r, v, w, s)
					)
					- 2. * (
					+(yu(s,v).real() * yudag(w,t).real() - yu(s,v).imag() * yudag(w,t).imag()) * WC7R(Clq1R, p, r, v, w)
					-(yu(s,v).imag() * yudag(w,t).real() + yu(s,v).real() * yudag(w,t).imag()) * WC7I(Clq1I, p, r, v, w)
					)
					- 6. * (
					+(ye(v,r).real() * yu(s,w).real() - ye(v,r).imag() * yu(s,w).imag()) * WC5(Clequ3R, p, v, w, t)
					-(ye(v,r).imag() * yu(s,w).real() + ye(v,r).real() * yu(s,w).imag()) * WC5(Clequ3I, p, v, w, t)
					+(yedag(p,v).real() * yudag(w,t).real() - yedag(p,v).imag() * yudag(w,t).imag()) * WC5(Clequ3R, r, v, w, s)
					+(yedag(p,v).imag() * yudag(w,t).real() + yedag(p,v).real() * yudag(w,t).imag()) * WC5(Clequ3I, r, v, w, s)
					)
					-(
					+(ye(w,r).real() * yedag(p,v).real() - ye(w,r).imag() * yedag(p,v).imag()) * WC7R(CeuR, v, w, s, t)
					-(ye(w,r).imag() * yedag(p,v).real() + ye(w,r).real() * yedag(p,v).imag()) * WC7I(CeuI, v, w, s, t)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CluI
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			-(yedye(p,r).imag() * WC2R(CHuR, s, t) + yedye(p,r).real() * WC2I(CHuI, s, t))
			- 2. * (gu(s,t).imag() * WC2R(CHl1R, p, r) + gu(s,t).real() * WC2I(CHl1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYl * WC2I(CHuI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYu * WC2I(CHl1I, p, r) * delta[s][t]
			- 12. * g12 * YuYl * WC7I(CluI, p, r, s, t)
			; //RGE 3

		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).imag() * WC7R(CluR, v, r, s, t) + gl(p,v).real() * WC7I(CluI, v, r, s, t))
				+(gu(s,v).imag() * WC7R(CluR, p, r, v, t) + gu(s,v).real() * WC7I(CluI, p, r, v, t))
				+(WC7R(CluR, p, v, s, t) * gl(v,r).imag() + WC7I(CluI, p, v, s, t) * gl(v,r).real())
				+(WC7R(CluR, p, r, s, v) * gu(v,t).imag() + WC7I(CluI, p, r, s, v) * gu(v,t).real())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YuYl * delta[s][t]*(
				+WC6I(CllI, p, r, w, w) + WC6I(CllI, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YuYl * delta[s][t]*(
				+WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * NC * YuYq * WC7I(Clq1I, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YlYq * WC7I(Cqu1I, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * Yu2 * WC7I(CluI, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * Yl2 * WC7I(CluI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7I(CldI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * YuYe * WC7I(CleI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YdYl * WC7I(Cud1I, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYl * WC7I(CeuI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYl * delta[p][r]* (
				+WC6I(CuuI, s, t, w, w) + WC6I(CuuI, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YuYl * delta[p][r]* (
				+WC6I(CuuI, s, w, w, t) + WC6I(CuuI, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					-0.5 * (
					+(ye(v,r).real() * yu(s,w).real() - ye(v,r).imag() * yu(s,w).imag()) * WC5(Clequ1I, p, v, w, t)
					+(ye(v,r).imag() * yu(s,w).real() + ye(v,r).real() * yu(s,w).imag()) * WC5(Clequ1R, p, v, w, t)
					-(yedag(p,v).real() * yudag(w,t).real() - yedag(p,v).imag() * yudag(w,t).imag()) * WC5(Clequ1I, r, v, w, s)
					+(yedag(p,v).imag() * yudag(w,t).real() + yedag(p,v).real() * yudag(w,t).imag()) * WC5(Clequ1R, r, v, w, s)
					)
					- 2. * (
					+(yu(s,v).real() * yudag(w,t).real() - yu(s,v).imag() * yudag(w,t).imag()) * WC7I(Clq1I, p, r, v, w)
					+(yu(s,v).imag() * yudag(w,t).real() + yu(s,v).real() * yudag(w,t).imag()) * WC7R(Clq1R, p, r, v, w)
					)
					- 6. * (
					+(ye(v,r).real() * yu(s,w).real() - ye(v,r).imag() * yu(s,w).imag()) * WC5(Clequ3I, p, v, w, t)
					+(ye(v,r).imag() * yu(s,w).real() + ye(v,r).real() * yu(s,w).imag()) * WC5(Clequ3R, p, v, w, t)
					-(yedag(p,v).real() * yudag(w,t).real() - yedag(p,v).imag() * yudag(w,t).imag()) * WC5(Clequ3I, r, v, w, s)
					+(yedag(p,v).imag() * yudag(w,t).real() + yedag(p,v).real() * yudag(w,t).imag()) * WC5(Clequ3R, r, v, w, s)
					)
					-(
					+(ye(w,r).real() * yedag(p,v).real() - ye(w,r).imag() * yedag(p,v).imag()) * WC7I(CeuI, v, w, s, t)
					+(ye(w,r).imag() * yedag(p,v).real() + ye(w,r).real() * yedag(p,v).imag()) * WC7R(CeuR, v, w, s, t)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}

	//CldR
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			-(yedye(p,r).real() * WC2R(CHdR, s, t) - yedye(p,r).imag() * WC2I(CHdI, s, t))
			+ 2. * (gd(s,t).real() * WC2R(CHl1R, p, r) - gd(s,t).imag() * WC2I(CHl1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYl * WC2R(CHdR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYd * WC2R(CHl1R, p, r) * delta[s][t]
			- 12. * g12 * YdYl * WC7R(CldR, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).real() * WC7R(CldR, v, r, s, t) - gl(p,v).imag() * WC7I(CldI, v, r, s, t))
				+(gd(s,v).real() * WC7R(CldR, p, r, v, t) - gd(s,v).imag() * WC7I(CldI, p, r, v, t))
				+(WC7R(CldR, p, v, s, t) * gl(v,r).real() - WC7I(CldI, p, v, s, t) * gl(v,r).imag())
				+(WC7R(CldR, p, r, s, v) * gd(v,t).real() - WC7I(CldI, p, r, s, v) * gd(v,t).imag())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YdYl * delta[s][t]*(
				+WC6R(CllR, p, r, w, w) + WC6R(CllR, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YdYl * delta[s][t]*(
				+WC6R(CllR, p, w, w, r) + WC6R(CllR, w, r, p, w)
				)
				+EIGHT_THIRDS * g12 * NC * YdYq * WC7R(Clq1R, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YlYq * WC7R(Cqd1R, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * Yd2 * WC7R(CldR, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * Yl2 * WC7R(CldR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7R(CluR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * YdYe * WC7R(CleR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYl * WC7R(Cud1R, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYl * WC7R(CedR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYl * delta[p][r]*(
				+WC6R(CddR, s, t, w, w) + WC6R(CddR, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YdYl * delta[p][r]*(
				+WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					-0.5 * (
					+(ye(v,r).real() * yddag(w,t).real() - ye(v,r).imag() * yddag(w,t).imag()) * WC5(CledqR, p, v, s, w)
					-(ye(v,r).imag() * yddag(w,t).real() + ye(v,r).real() * yddag(w,t).imag()) * WC5(CledqI, p, v, s, w)
					+(yedag(p,v).real() * yd(s,w).real() - yedag(p,v).imag() * yd(s,w).imag()) * WC5(CledqR, r, v, t, w)
					+(yedag(p,v).imag() * yd(s,w).real() + yedag(p,v).real() * yd(s,w).imag()) * WC5(CledqI, r, v, t, w)
					)
					- 2. * (
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7R(Clq1R, p, r, v, w)
					-(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7I(Clq1I, p, r, v, w)
					)
					- (
					+(ye(w,r).real() * yedag(p,v).real() - ye(w,r).imag() * yedag(p,v).imag()) * WC7R(CedR, v, w, s, t)
					-(ye(w,r).imag() * yedag(p,v).real() + ye(w,r).real() * yedag(p,v).imag()) * WC7I(CedI, v, w, s, t)
					)

					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//CldI
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			-(yedye(p,r).imag() * WC2R(CHdR, s, t) + yedye(p,r).real() * WC2I(CHdI, s, t))
			+ 2. * (gd(s,t).imag() * WC2R(CHl1R, p, r) + gd(s,t).real() * WC2I(CHl1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYl * WC2I(CHdI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYd * WC2I(CHl1I, p, r) * delta[s][t]
			- 12. * g12 * YdYl * WC7I(CldI, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gl(p,v).imag() * WC7R(CldR, v, r, s, t) + gl(p,v).real() * WC7I(CldI, v, r, s, t))
				+(gd(s,v).imag() * WC7R(CldR, p, r, v, t) + gd(s,v).real() * WC7I(CldI, p, r, v, t))
				+(WC7R(CldR, p, v, s, t) * gl(v,r).imag() + WC7I(CldI, p, v, s, t) * gl(v,r).real())
				+(WC7R(CldR, p, r, s, v) * gd(v,t).imag() + WC7I(CldI, p, r, s, v) * gd(v,t).real())
				//RGE 2
				+ EIGHT_THIRDS * g12 * YdYl * delta[s][t]*(
				+WC6I(CllI, p, r, w, w) + WC6I(CllI, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YdYl * delta[s][t]*(
				+WC6I(CllI, p, w, w, r) + WC6I(CllI, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * NC * YdYq * WC7I(Clq1I, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * YlYq * WC7I(Cqd1I, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * Yd2 * WC7I(CldI, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * Yl2 * WC7I(CldI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7I(CluI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * YdYe * WC7I(CleI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYl * WC7I(Cud1I, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYl * WC7I(CedI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYl * delta[p][r]*(
				+WC6I(CddI, s, t, w, w) + WC6I(CddI, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YdYl * delta[p][r]*(
				+WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					-0.5 * (
					+(ye(v,r).real() * yddag(w,t).real() - ye(v,r).imag() * yddag(w,t).imag()) * WC5(CledqI, p, v, s, w)
					+(ye(v,r).imag() * yddag(w,t).real() + ye(v,r).real() * yddag(w,t).imag()) * WC5(CledqR, p, v, s, w)
					-(yedag(p,v).real() * yd(s,w).real() - yedag(p,v).imag() * yd(s,w).imag()) * WC5(CledqI, r, v, t, w)
					+(yedag(p,v).imag() * yd(s,w).real() + yedag(p,v).real() * yd(s,w).imag()) * WC5(CledqR, r, v, t, w)
					)
					- 2. * (
					+(yd(s,v).real() * yddag(w,t).real() - yd(s,v).imag() * yddag(w,t).imag()) * WC7I(Clq1I, p, r, v, w)
					+(yd(s,v).imag() * yddag(w,t).real() + yd(s,v).real() * yddag(w,t).imag()) * WC7R(Clq1R, p, r, v, w)
					)
					- (
					+(ye(w,r).real() * yedag(p,v).real() - ye(w,r).imag() * yedag(p,v).imag()) * WC7I(CedI, v, w, s, t)
					+(ye(w,r).imag() * yedag(p,v).real() + ye(w,r).real() * yedag(p,v).imag()) * WC7R(CedR, v, w, s, t)
					)

					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}


	//CqeR
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] = (
			+(yudyu(p,r).real() - yddyd(p,r).real()) * WC2R(CHeR, s, t)
			-(yudyu(p,r).imag() - yddyd(p,r).imag()) * WC2I(CHeI, s, t)
			)
			+ 2. * (ge(s,t).real() * WC2R(CHq1R, p, r) - ge(s,t).imag() * WC2I(CHq1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYq * WC2R(CHeR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYe * WC2R(CHq1R, p, r) * delta[s][t]
			- 12. * YeYq * g12 * WC7R(CqeR, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).real() * WC7R(CqeR, v, r, s, t) - gq(p,v).imag() * WC7I(CqeI, v, r, s, t))
				+ (ge(s,v).real() * WC7R(CqeR, p, r, v, t) - ge(s,v).imag() * WC7I(CqeI, p, r, v, t))
				+(WC7R(CqeR, p, v, s, t) * gq(v,r).real() - WC7I(CqeI, p, v, s, t) * gq(v,r).imag())
				+(WC7R(CqeR, p, r, s, v) * ge(v,t).real() - WC7I(CqeI, p, r, s, v) * ge(v,t).imag())
				//RGE 2
				+ EIGHT_THIRDS * g12 * NC * YeYq * delta[s][t]*(
				+WC6R(Cqq1R, p, r, w, w) + WC6R(Cqq1R, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YeYq * delta[s][t]*(
				+WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
				)
				+ 4. * g12 * YeYq * delta[s][t]*(
				+WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * YeYl * WC7R(Clq1R, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * YlYq * WC7R(CleR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * Ye2 * WC7R(CqeR, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * Yq2 * WC7R(CqeR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYe * WC7R(Cqu1R, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YdYe * WC7R(Cqd1R, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYq * WC7R(CeuR, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYq * WC7R(CedR, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYq * delta[p][r]*(
				+WC8R(CeeR, s, t, w, w) + WC8R(CeeR, w, w, s, t)
				+ WC8R(CeeR, s, w, w, t) + WC8R(CeeR, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -0.5 * (
					+(yddag(p,w).real() * ye(s,v).real() - yddag(p,w).imag() * ye(s,v).imag()) * WC5(CledqR, v, t, w, r)
					-(yddag(p,w).imag() * ye(s,v).real() + yddag(p,w).real() * ye(s,v).imag()) * WC5(CledqI, v, t, w, r)
					+(yedag(v,t).real() * yd(w,r).real() - yedag(v,t).imag() * yd(w,r).imag()) * WC5(CledqR, v, s, w, p)
					+(yedag(v,t).imag() * yd(w,r).real() + yedag(v,t).real() * yd(w,r).imag()) * WC5(CledqI, v, s, w, p)
					)
					- 2. * (
					+(ye(s,v).real() * yedag(w,t).real() - ye(s,v).imag() * yedag(w,t).imag()) * WC7R(Clq1R, v, w, p, r)
					-(ye(s,v).imag() * yedag(w,t).real() + ye(s,v).real() * yedag(w,t).imag()) * WC7I(Clq1I, v, w, p, r)
					)
					- 0.5 * (
					+(yu(w,r).real() * ye(s,v).real() - yu(w,r).imag() * ye(s,v).imag()) * WC5(Clequ1R, v, t, p, w)
					-(yu(w,r).imag() * ye(s,v).real() + yu(w,r).real() * ye(s,v).imag()) * WC5(Clequ1I, v, t, p, w)
					+(yedag(v,t).real() * yudag(p,w).real() - yedag(v,t).imag() * yudag(p,w).imag()) * WC5(Clequ1R, v, s, r, w)
					+(yedag(v,t).imag() * yudag(p,w).real() + yedag(v,t).real() * yudag(p,w).imag()) * WC5(Clequ1I, v, s, r, w)
					)
					-(
					+(yd(w,r).real() * yddag(p,v).real() - yd(w,r).imag() * yddag(p,v).imag()) * WC7R(CedR, s, t, v, w)
					-(yd(w,r).imag() * yddag(p,v).real() + yd(w,r).real() * yddag(p,v).imag()) * WC7I(CedI, s, t, v, w)
					)
					- 6. * (
					+(yu(w,r).real() * ye(s,v).real() - yu(w,r).imag() * ye(s,v).imag()) * WC5(Clequ3R, v, t, p, w)
					-(yu(w,r).imag() * ye(s,v).real() + yu(w,r).real() * ye(s,v).imag()) * WC5(Clequ3I, v, t, p, w)
					+(yedag(v,t).real() * yudag(p,w).real() - yedag(v,t).imag() * yudag(p,w).imag()) * WC5(Clequ3R, v, s, r, w)
					+(yedag(v,t).imag() * yudag(p,w).real() + yedag(v,t).real() * yudag(p,w).imag()) * WC5(Clequ3I, v, s, r, w)
					)
					-(
					+(yu(w,r).real() * yudag(p,v).real() - yu(w,r).imag() * yudag(p,v).imag()) * WC7R(CeuR, s, t, v, w)
					-(yu(w,r).imag() * yudag(p,v).real() + yu(w,r).real() * yudag(p,v).imag()) * WC7I(CeuI, s, t, v, w)
					)
					; //RGE 2
			}
		}


		f[c] *= loop_factor;
		c++;
	}

	//CqeI
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] = (
			+(yudyu(p,r).imag() - yddyd(p,r).imag()) * WC2R(CHeR, s, t)
			+(yudyu(p,r).real() - yddyd(p,r).real()) * WC2I(CHeI, s, t)
			)
			+ 2. * (ge(s,t).imag() * WC2R(CHq1R, p, r) + ge(s,t).real() * WC2I(CHq1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYq * WC2I(CHeI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYe * WC2I(CHq1I, p, r) * delta[s][t]
			- 12. * YeYq * g12 * WC7I(CqeI, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).imag() * WC7R(CqeR, v, r, s, t) + gq(p,v).real() * WC7I(CqeI, v, r, s, t))
				+ (ge(s,v).imag() * WC7R(CqeR, p, r, v, t) + ge(s,v).real() * WC7I(CqeI, p, r, v, t))
				+(WC7R(CqeR, p, v, s, t) * gq(v,r).imag() + WC7I(CqeI, p, v, s, t) * gq(v,r).real())
				+(WC7R(CqeR, p, r, s, v) * ge(v,t).imag() + WC7I(CqeI, p, r, s, v) * ge(v,t).real())
				//RGE 2
				+ EIGHT_THIRDS * g12 * NC * YeYq * delta[s][t]*(
				+WC6I(Cqq1I, p, r, w, w) + WC6I(Cqq1I, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YeYq * delta[s][t]*(
				+WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
				)
				+ 4. * g12 * YeYq * delta[s][t]*(
				+WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * YeYl * WC7I(Clq1I, w, w, p, r) * delta[s][t]
				+ EIGHT_THIRDS * g12 * YlYq * WC7I(CleI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * Ye2 * WC7I(CqeI, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * Yq2 * WC7I(CqeI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYe * WC7I(Cqu1I, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YdYe * WC7I(Cqd1I, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYq * WC7I(CeuI, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYq * WC7I(CedI, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYq * delta[p][r]*(
				+WC8I(CeeI, s, t, w, w) + WC8I(CeeI, w, w, s, t)
				+ WC8I(CeeI, s, w, w, t) + WC8I(CeeI, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -0.5 * (
					+(yddag(p,w).real() * ye(s,v).real() - yddag(p,w).imag() * ye(s,v).imag()) * WC5(CledqI, v, t, w, r)
					+(yddag(p,w).imag() * ye(s,v).real() + yddag(p,w).real() * ye(s,v).imag()) * WC5(CledqR, v, t, w, r)
					-(yedag(v,t).real() * yd(w,r).real() - yedag(v,t).imag() * yd(w,r).imag()) * WC5(CledqI, v, s, w, p)
					+(yedag(v,t).imag() * yd(w,r).real() + yedag(v,t).real() * yd(w,r).imag()) * WC5(CledqR, v, s, w, p)
					)
					- 2. * (
					+(ye(s,v).real() * yedag(w,t).real() - ye(s,v).imag() * yedag(w,t).imag()) * WC7I(Clq1I, v, w, p, r)
					+(ye(s,v).imag() * yedag(w,t).real() + ye(s,v).real() * yedag(w,t).imag()) * WC7R(Clq1R, v, w, p, r)
					)
					- 0.5 * (
					+(yu(w,r).real() * ye(s,v).real() - yu(w,r).imag() * ye(s,v).imag()) * WC5(Clequ1I, v, t, p, w)
					+(yu(w,r).imag() * ye(s,v).real() + yu(w,r).real() * ye(s,v).imag()) * WC5(Clequ1R, v, t, p, w)
					-(yedag(v,t).real() * yudag(p,w).real() - yedag(v,t).imag() * yudag(p,w).imag()) * WC5(Clequ1I, v, s, r, w)
					+(yedag(v,t).imag() * yudag(p,w).real() + yedag(v,t).real() * yudag(p,w).imag()) * WC5(Clequ1R, v, s, r, w)
					)
					-(
					+(yd(w,r).real() * yddag(p,v).real() - yd(w,r).imag() * yddag(p,v).imag()) * WC7I(CedI, s, t, v, w)
					+(yd(w,r).imag() * yddag(p,v).real() + yd(w,r).real() * yddag(p,v).imag()) * WC7R(CedR, s, t, v, w)
					)
					- 6. * (
					+(yu(w,r).real() * ye(s,v).real() - yu(w,r).imag() * ye(s,v).imag()) * WC5(Clequ3I, v, t, p, w)
					+(yu(w,r).imag() * ye(s,v).real() + yu(w,r).real() * ye(s,v).imag()) * WC5(Clequ3R, v, t, p, w)
					-(yedag(v,t).real() * yudag(p,w).real() - yedag(v,t).imag() * yudag(p,w).imag()) * WC5(Clequ3I, v, s, r, w)
					+(yedag(v,t).imag() * yudag(p,w).real() + yedag(v,t).real() * yudag(p,w).imag()) * WC5(Clequ3R, v, s, r, w)
					)
					-(
					+(yu(w,r).real() * yudag(p,v).real() - yu(w,r).imag() * yudag(p,v).imag()) * WC7I(CeuI, s, t, v, w)
					+(yu(w,r).imag() * yudag(p,v).real() + yu(w,r).real() * yudag(p,v).imag()) * WC7R(CeuR, s, t, v, w)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}

	//Cqu1R
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			+(1. / NC)*(yu(s,r).real() * xiu(p,t).real() - yu(s,r).imag() * xiu(p,t).imag())
			+(1. / NC)*(yudag(p,t).real() * xiu(r,s).real() + yudag(p,t).imag() * xiu(r,s).imag())
			+(
			+(yudyu(p,r).real() - yddyd(p,r).real()) * WC2R(CHuR, s, t)
			-(yudyu(p,r).imag() - yddyd(p,r).imag()) * WC2I(CHuI, s, t)
			)
			- 2. * (gu(s,t).real() * WC2R(CHq1R, p, r) - gu(s,t).imag() * WC2I(CHq1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYq * WC2R(CHuR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYu * WC2R(CHq1R, p, r) * delta[s][t]
			- 12. * YuYq * g12 * WC7R(Cqu1R, p, r, s, t)
			- 3. * ((NC2 - 1.) / NC2) * g32 * WC7R(Cqu8R, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).real() * WC7R(Cqu1R, v, r, s, t) - gq(p,v).imag() * WC7I(Cqu1I, v, r, s, t))
				+(gu(s,v).real() * WC7R(Cqu1R, p, r, v, t) - gu(s,v).imag() * WC7I(Cqu1I, p, r, v, t))
				+(WC7R(Cqu1R, p, v, s, t) * gq(v,r).real() - WC7I(Cqu1I, p, v, s, t) * gq(v,r).imag())
				+(WC7R(Cqu1R, p, r, s, v) * gu(v,t).real() - WC7I(Cqu1I, p, r, s, v) * gu(v,t).imag())
				//RGE 2
				+ EIGHT_THIRDS * g12 * NC * YuYq * delta[s][t]*
				(+WC6R(Cqq1R, p, r, w, w) + WC6R(Cqq1R, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YuYq * delta[s][t]*
				(+WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
				)
				+ 4. * g12 * YuYq * delta[s][t]* (
				+WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * YuYl * WC7R(Clq1R, w, w, p, r) * delta[s][t]
				+ FOUR_THIRDS * g12 * YuYe * WC7R(CqeR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7R(Cqd1R, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * Yu2 * WC7R(Cqu1R, p, r, w, w) * delta[s][t]

				+ EIGHT_THIRDS * g12 * NC * Yq2 * WC7R(Cqu1R, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * YlYq * WC7R(CluR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYq * WC7R(CeuR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYq * WC7R(Cud1R, s, t, w, w) * delta[p][r]

				+ FOUR_THIRDS * g12 * NC * YuYq * delta[p][r]*(
				+WC6R(CuuR, s, t, w, w) + WC6R(CuuR, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YuYq * delta[p][r]*(
				+WC6R(CuuR, s, w, w, t) + WC6R(CuuR, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += +(1. / NC)*(
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7R(Cqu1R, v, r, w, t)
					-(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7I(Cqu1I, v, r, w, t)
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7R(Cqu1R, p, v, s, w)
					-(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7I(Cqu1I, p, v, s, w)
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd1R, p, t, v, w)
					-(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd1I, p, t, v, w)
					+(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd1R, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd1I, r, s, v, w)
					)
					-(0.5 / NC2)*(
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7R(Cqu8R, v, r, w, t)
					-(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7I(Cqu8I, v, r, w, t)
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7R(Cqu8R, p, v, s, w)
					-(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7I(Cqu8I, p, v, s, w)
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd8R, p, t, v, w)
					-(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd8I, p, t, v, w)
					+(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd8R, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd8I, r, s, v, w)
					)
					-(2. / NC)*(
					+(yudag(v,t).real() * yu(s,w).real() - yudag(v,t).imag() * yu(s,w).imag()) * WC6R(Cqq1R, p, v, w, r)
					-(yudag(v,t).imag() * yu(s,w).real() + yudag(v,t).real() * yu(s,w).imag()) * WC6I(Cqq1I, p, v, w, r)
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC6R(CuuR, v, t, s, w)
					-(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC6I(CuuI, v, t, s, w)
					)
					-(6. / NC)*(
					+(yudag(v,t).real() * yu(s,w).real() - yudag(v,t).imag() * yu(s,w).imag()) * WC6R(Cqq3R, p, v, w, r)
					-(yudag(v,t).imag() * yu(s,w).real() + yudag(v,t).real() * yu(s,w).imag()) * WC6I(Cqq3I, p, v, w, r)
					)
					+ 0.5 * (
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7R(Cqu8R, v, r, w, t)
					-(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7I(Cqu8I, v, r, w, t)
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7R(Cqu8R, p, v, s, w)
					-(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7I(Cqu8I, p, v, s, w)
					)
					+ 0.5 * (
					+(yu(s,v).real() * yd(w,r).real() - yu(s,v).imag() * yd(w,r).imag()) * WC5(Cquqd1R, v, t, p, w)
					-(yu(s,v).imag() * yd(w,r).real() + yu(s,v).real() * yd(w,r).imag()) * WC5(Cquqd1I, v, t, p, w)
					+(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd1R, v, s, r, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd1I, v, s, r, w)
					+(yu(s,v).real() * yd(w,r).real() - yu(s,v).imag() * yd(w,r).imag()) * WC5(Cquqd8R, p, t, v, w)
					-(yu(s,v).imag() * yd(w,r).real() + yu(s,v).real() * yd(w,r).imag()) * WC5(Cquqd8I, p, t, v, w)
					+(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd8R, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd8I, r, s, v, w)
					)

					- 4. * (
					+(yudag(w,t).real() * yu(s,v).real() - yudag(w,t).imag() * yu(s,v).imag()) * WC6R(Cqq1R, p, r, v, w)
					-(yudag(w,t).imag() * yu(s,v).real() + yudag(w,t).real() * yu(s,v).imag()) * WC6I(Cqq1I, p, r, v, w)
					)
					- 2. * (
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC6R(CuuR, v, w, s, t)
					-(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC6I(CuuI, v, w, s, t)
					)
					-(
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC7R(Cud1R, s, t, v, w)
					-(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC7I(Cud1I, s, t, v, w)
					)



					; //RGE 2
			}

		}
		f[c] *= loop_factor;
		c++;
	}

	//Cqu1I
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			+(1. / NC)*(yu(s,r).imag() * xiu(p,t).real() + yu(s,r).real() * xiu(p,t).imag())
			+(1. / NC)*(yudag(p,t).imag() * xiu(r,s).real() - yudag(p,t).real() * xiu(r,s).imag())
			+(
			+(yudyu(p,r).imag() - yddyd(p,r).imag()) * WC2R(CHuR, s, t)
			+(yudyu(p,r).real() - yddyd(p,r).real()) * WC2I(CHuI, s, t)
			)
			- 2. * (gu(s,t).imag() * WC2R(CHq1R, p, r) + gu(s,t).real() * WC2I(CHq1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYq * WC2I(CHuI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYu * WC2I(CHq1I, p, r) * delta[s][t]
			- 12. * YuYq * g12 * WC7I(Cqu1I, p, r, s, t)
			- 3. * ((NC2 - 1.) / NC2) * g32 * WC7I(Cqu8I, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).imag() * WC7R(Cqu1R, v, r, s, t) + gq(p,v).real() * WC7I(Cqu1I, v, r, s, t))
				+(gu(s,v).imag() * WC7R(Cqu1R, p, r, v, t) + gu(s,v).real() * WC7I(Cqu1I, p, r, v, t))
				+(WC7R(Cqu1R, p, v, s, t) * gq(v,r).imag() + WC7I(Cqu1I, p, v, s, t) * gq(v,r).real())
				+(WC7R(Cqu1R, p, r, s, v) * gu(v,t).imag() + WC7I(Cqu1I, p, r, s, v) * gu(v,t).real())
				//RGE 2
				+ EIGHT_THIRDS * g12 * NC * YuYq * delta[s][t]*
				(+WC6I(Cqq1I, p, r, w, w) + WC6I(Cqq1I, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YuYq * delta[s][t]*
				(+WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
				)
				+ 4. * g12 * YuYq * delta[s][t]* (
				+WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * YuYl * WC7I(Clq1I, w, w, p, r) * delta[s][t]
				+ FOUR_THIRDS * g12 * YuYe * WC7I(CqeI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7I(Cqd1I, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * Yu2 * WC7I(Cqu1I, p, r, w, w) * delta[s][t]

				+ EIGHT_THIRDS * g12 * NC * Yq2 * WC7I(Cqu1I, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * YlYq * WC7I(CluI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYq * WC7I(CeuI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YdYq * WC7I(Cud1I, s, t, w, w) * delta[p][r]

				+ FOUR_THIRDS * g12 * NC * YuYq * delta[p][r]*(
				+WC6I(CuuI, s, t, w, w) + WC6I(CuuI, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YuYq * delta[p][r]*(
				+WC6I(CuuI, s, w, w, t) + WC6I(CuuI, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += +(1. / NC)*(
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7I(Cqu1I, v, r, w, t)
					+(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7R(Cqu1R, v, r, w, t)
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7I(Cqu1I, p, v, s, w)
					+(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7R(Cqu1R, p, v, s, w)
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd1I, p, t, v, w)
					+(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd1R, p, t, v, w)
					-(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd1I, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd1R, r, s, v, w)
					)
					-(0.5 / NC2)*(
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7I(Cqu8I, v, r, w, t)
					+(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7R(Cqu8R, v, r, w, t)
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7I(Cqu8I, p, v, s, w)
					+(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7R(Cqu8R, p, v, s, w)
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd8I, p, t, v, w)
					+(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd8R, p, t, v, w)
					-(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd8I, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd8R, r, s, v, w)
					)
					-(2. / NC)*(
					+(yudag(v,t).real() * yu(s,w).real() - yudag(v,t).imag() * yu(s,w).imag()) * WC6I(Cqq1I, p, v, w, r)
					+(yudag(v,t).imag() * yu(s,w).real() + yudag(v,t).real() * yu(s,w).imag()) * WC6R(Cqq1R, p, v, w, r)
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC6I(CuuI, v, t, s, w)
					+(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC6R(CuuR, v, t, s, w)
					)
					-(6. / NC)*(
					+(yudag(v,t).real() * yu(s,w).real() - yudag(v,t).imag() * yu(s,w).imag()) * WC6I(Cqq3I, p, v, w, r)
					+(yudag(v,t).imag() * yu(s,w).real() + yudag(v,t).real() * yu(s,w).imag()) * WC6R(Cqq3R, p, v, w, r)
					)
					+ 0.5 * (
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7I(Cqu8I, v, r, w, t)
					+(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7R(Cqu8R, v, r, w, t)
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7I(Cqu8I, p, v, s, w)
					+(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7R(Cqu8R, p, v, s, w)
					)
					+ 0.5 * (
					+(yu(s,v).real() * yd(w,r).real() - yu(s,v).imag() * yd(w,r).imag()) * WC5(Cquqd1I, v, t, p, w)
					+(yu(s,v).imag() * yd(w,r).real() + yu(s,v).real() * yd(w,r).imag()) * WC5(Cquqd1R, v, t, p, w)
					-(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd1I, v, s, r, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd1R, v, s, r, w)
					+(yu(s,v).real() * yd(w,r).real() - yu(s,v).imag() * yd(w,r).imag()) * WC5(Cquqd8I, p, t, v, w)
					+(yu(s,v).imag() * yd(w,r).real() + yu(s,v).real() * yd(w,r).imag()) * WC5(Cquqd8R, p, t, v, w)
					-(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd8I, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd8R, r, s, v, w)
					)

					- 4. * (
					+(yudag(w,t).real() * yu(s,v).real() - yudag(w,t).imag() * yu(s,v).imag()) * WC6I(Cqq1I, p, r, v, w)
					+(yudag(w,t).imag() * yu(s,v).real() + yudag(w,t).real() * yu(s,v).imag()) * WC6R(Cqq1R, p, r, v, w)
					)
					- 2. * (
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC6I(CuuI, v, w, s, t)
					+(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC6R(CuuR, v, w, s, t)
					)
					-(
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC7I(Cud1I, s, t, v, w)
					+(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC7R(Cud1R, s, t, v, w)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}

	//Cqu8R
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			+2. * (yu(s,r).real() * xiu(p,t).real() - yu(s,r).imag() * xiu(p,t).imag())
			+ 2. * (yudag(p,t).real() * xiu(r,s).real() + yudag(p,t).imag() * xiu(r,s).imag())
			//RGE 2
			-(12. * YuYq * g12 + 6. * (NC - (2. / NC)) * g32) * WC7R(Cqu8R, p, r, s, t)
			- 12. * g32 * WC7R(Cqu1R, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).real() * WC7R(Cqu8R, v, r, s, t) - gq(p,v).imag() * WC7I(Cqu8I, v, r, s, t))
				+(gu(s,v).real() * WC7R(Cqu8R, p, r, v, t) - gu(s,v).imag() * WC7I(Cqu8I, p, r, v, t))
				+(WC7R(Cqu8R, p, v, s, t) * gq(v,r).real() - WC7I(Cqu8I, p, v, s, t) * gq(v,r).imag())
				+(WC7R(Cqu8R, p, r, s, v) * gu(v,t).real() - WC7I(Cqu8I, p, r, s, v) * gu(v,t).imag())
				//RGE 2
				+ FOUR_THIRDS * g32 * delta[s][t]*(
				+WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
				)
				+ 4. * g32 * delta[s][t]*(
				+WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
				)
				+ TWO_THIRDS * g32 * WC7R(Cqu8R, p, r, w, w) * delta[s][t]
				+ TWO_THIRDS * g32 * WC7R(Cqd8R, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g32 * WC7R(Cqu8R, w, w, s, t) * delta[p][r]
				+ TWO_THIRDS * g32 * WC7R(Cud8R, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g32 * delta[p][r]*(
				+WC6R(CuuR, s, w, w, t) + WC6R(CuuR, w, t, s, w)
				)
				; // RGE 3
			for (w = 0; w < NG; w++) {
				f[c] +=
					-(1. / NC)*(
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7R(Cqu8R, v, r, w, t)
					-(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7I(Cqu8I, v, r, w, t)
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7R(Cqu8R, p, v, s, w)
					-(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7I(Cqu8I, p, v, s, w)
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd8R, p, t, v, w)
					-(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd8I, p, t, v, w)
					+(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd8R, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd8I, r, s, v, w)
					)

					+ 2. * (
					(
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd1R, p, t, v, w)
					-(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd1I, p, t, v, w)
					+(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd1R, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd1I, r, s, v, w)
					)
					+ 0.25 * (
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd8R, v, t, p, w)
					-(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd8I, v, t, p, w)
					+(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd8R, v, s, r, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd8I, v, s, r, w)
					)
					)

					- 2. * (
					+2. * (
					+(yudag(v,t).real() * yu(s,w).real() - yudag(v,t).imag() * yu(s,w).imag()) * WC6R(Cqq1R, p, v, w, r)
					-(yudag(v,t).imag() * yu(s,w).real() + yudag(v,t).real() * yu(s,w).imag()) * WC6I(Cqq1I, p, v, w, r)
					)
					-(
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7R(Cqu1R, v, r, w, t)
					-(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7I(Cqu1I, v, r, w, t)
					)
					-(
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7R(Cqu1R, p, v, s, w)
					-(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7I(Cqu1I, p, v, s, w)
					)
					)

					- 4. * (
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC6R(CuuR, v, t, s, w)
					-(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC6I(CuuI, v, t, s, w)
					)
					- 12. * (
					+(yudag(v,t).real() * yu(s,w).real() - yudag(v,t).imag() * yu(s,w).imag()) * WC6R(Cqq3R, p, v, w, r)
					-(yudag(v,t).imag() * yu(s,w).real() + yudag(v,t).real() * yu(s,w).imag()) * WC6I(Cqq3I, p, v, w, r)
					)
					-(
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC7R(Cud8R, s, t, v, w)
					-(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC7I(Cud8I, s, t, v, w)
					)

					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}

	//Cqu8I
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			+2. * (yu(s,r).imag() * xiu(p,t).real() + yu(s,r).real() * xiu(p,t).imag())
			+ 2. * (yudag(p,t).imag() * xiu(r,s).real() - yudag(p,t).real() * xiu(r,s).imag())
			//RGE 2
			-(12. * YuYq * g12 + 6. * (NC - (2. / NC)) * g32) * WC7I(Cqu8I, p, r, s, t)
			- 12. * g32 * WC7I(Cqu1I, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).imag() * WC7R(Cqu8R, v, r, s, t) + gq(p,v).real() * WC7I(Cqu8I, v, r, s, t))
				+(gu(s,v).imag() * WC7R(Cqu8R, p, r, v, t) + gu(s,v).real() * WC7I(Cqu8I, p, r, v, t))
				+(WC7R(Cqu8R, p, v, s, t) * gq(v,r).imag() + WC7I(Cqu8I, p, v, s, t) * gq(v,r).real())
				+(WC7R(Cqu8R, p, r, s, v) * gu(v,t).imag() + WC7I(Cqu8I, p, r, s, v) * gu(v,t).real())
				//RGE 2
				+ FOUR_THIRDS * g32 * delta[s][t]*(
				+WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
				)
				+ 4. * g32 * delta[s][t]*(
				+WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
				)
				+ TWO_THIRDS * g32 * WC7I(Cqu8I, p, r, w, w) * delta[s][t]
				+ TWO_THIRDS * g32 * WC7I(Cqd8I, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g32 * WC7I(Cqu8I, w, w, s, t) * delta[p][r]
				+ TWO_THIRDS * g32 * WC7I(Cud8I, s, t, w, w) * delta[p][r]
				+ FOUR_THIRDS * g32 * delta[p][r]*(
				+WC6I(CuuI, s, w, w, t) + WC6I(CuuI, w, t, s, w)
				)
				; // RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -(1. / NC)*(
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7I(Cqu8I, v, r, w, t)
					+(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7R(Cqu8R, v, r, w, t)
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7I(Cqu8I, p, v, s, w)
					+(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7R(Cqu8R, p, v, s, w)
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd8I, p, t, v, w)
					+(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd8R, p, t, v, w)
					-(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd8I, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd8R, r, s, v, w)
					)

					+ 2. * (
					(
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd1I, p, t, v, w)
					+(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd1R, p, t, v, w)
					-(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd1I, r, s, v, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd1R, r, s, v, w)
					)
					+ 0.25 * (
					+(yd(w,r).real() * yu(s,v).real() - yd(w,r).imag() * yu(s,v).imag()) * WC5(Cquqd8I, v, t, p, w)
					+(yd(w,r).imag() * yu(s,v).real() + yd(w,r).real() * yu(s,v).imag()) * WC5(Cquqd8R, v, t, p, w)
					-(yddag(p,w).real() * yudag(v,t).real() - yddag(p,w).imag() * yudag(v,t).imag()) * WC5(Cquqd8I, v, s, r, w)
					+(yddag(p,w).imag() * yudag(v,t).real() + yddag(p,w).real() * yudag(v,t).imag()) * WC5(Cquqd8R, v, s, r, w)
					)
					)

					- 2. * (
					+2. * (
					+(yudag(v,t).real() * yu(s,w).real() - yudag(v,t).imag() * yu(s,w).imag()) * WC6I(Cqq1I, p, v, w, r)
					+(yudag(v,t).imag() * yu(s,w).real() + yudag(v,t).real() * yu(s,w).imag()) * WC6R(Cqq1R, p, v, w, r)
					)
					-(
					+(yudag(p,w).real() * yu(s,v).real() - yudag(p,w).imag() * yu(s,v).imag()) * WC7I(Cqu1I, v, r, w, t)
					+(yudag(p,w).imag() * yu(s,v).real() + yudag(p,w).real() * yu(s,v).imag()) * WC7R(Cqu1R, v, r, w, t)
					)
					-(
					+(yudag(v,t).real() * yu(w,r).real() - yudag(v,t).imag() * yu(w,r).imag()) * WC7I(Cqu1I, p, v, s, w)
					+(yudag(v,t).imag() * yu(w,r).real() + yudag(v,t).real() * yu(w,r).imag()) * WC7R(Cqu1R, p, v, s, w)
					)
					)

					- 4. * (
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC6I(CuuI, v, t, s, w)
					+(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC6R(CuuR, v, t, s, w)
					)
					- 12. * (
					+(yudag(v,t).real() * yu(s,w).real() - yudag(v,t).imag() * yu(s,w).imag()) * WC6I(Cqq3I, p, v, w, r)
					+(yudag(v,t).imag() * yu(s,w).real() + yudag(v,t).real() * yu(s,w).imag()) * WC6R(Cqq3R, p, v, w, r)
					)
					-(
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC7I(Cud8I, s, t, v, w)
					+(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC7R(Cud8R, s, t, v, w)
					)

					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}


	//Cqd1R
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			+(1. / NC)*(yd(s,r).real() * xid(p,t).real() - yd(s,r).imag() * xid(p,t).imag())
			+(1. / NC)*(yddag(p,t).real() * xid(r,s).real() + yddag(p,t).imag() * xid(r,s).imag())
			+(
			+(yudyu(p,r).real() - yddyd(p,r).real()) * WC2R(CHdR, s, t)
			-(yudyu(p,r).imag() - yddyd(p,r).imag()) * WC2I(CHdI, s, t)
			)
			+ 2. * (gd(s,t).real() * WC2R(CHq1R, p, r) - gd(s,t).imag() * WC2I(CHq1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYq * WC2R(CHdR, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYd * WC2R(CHq1R, p, r) * delta[s][t]
			- 12. * YdYq * g12 * WC7R(Cqd1R, p, r, s, t)
			- 3. * ((NC2 - 1.) / NC2) * g32 * WC7R(Cqd8R, p, r, s, t)
			; //RGE 3


		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).real() * WC7R(Cqd1R, v, r, s, t) - gq(p,v).imag() * WC7I(Cqd1I, v, r, s, t))
				+(gd(s,v).real() * WC7R(Cqd1R, p, r, v, t) - gd(s,v).imag() * WC7I(Cqd1I, p, r, v, t))
				+(WC7R(Cqd1R, p, v, s, t) * gq(v,r).real() - WC7I(Cqd1I, p, v, s, t) * gq(v,r).imag())
				+(WC7R(Cqd1R, p, r, s, v) * gd(v,t).real() - WC7I(Cqd1I, p, r, s, v) * gd(v,t).imag())
				//RGE 2
				+ EIGHT_THIRDS * g12 * NC * YdYq * delta[s][t]*(
				+WC6R(Cqq1R, p, r, w, w) + WC6R(Cqq1R, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YdYq * delta[s][t]*(
				+WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
				)
				+ 4. * g12 * YdYq * delta[s][t]*(
				+WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * YdYl * WC7R(Clq1R, w, w, p, r) * delta[s][t]
				+ FOUR_THIRDS * g12 * YdYe * WC7R(CqeR, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7R(Cqu1R, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * Yd2 * WC7R(Cqd1R, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * Yq2 * WC7R(Cqd1R, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * YlYq * WC7R(CldR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYq * WC7R(CedR, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYq * WC7R(Cud1R, w, w, s, t) * delta[p][r]

				+ FOUR_THIRDS * g12 * NC * YdYq * delta[p][r]*(
				+WC6R(CddR, s, t, w, w) + WC6R(CddR, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YdYq * delta[p][r]*(
				+WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w)
				)
				; //RGE 3


			for (w = 0; w < NG; w++) {
				f[c] += +(1. / NC)*(
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7R(Cqd1R, v, r, w, t)
					-(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7I(Cqd1I, v, r, w, t)
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7R(Cqd1R, p, v, s, w)
					-(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7I(Cqd1I, p, v, s, w)
					+(yu(w,r).real() * yd(s,v).real() - yu(w,r).imag() * yd(s,v).imag()) * WC5(Cquqd1R, v, w, p, t)
					-(yu(w,r).imag() * yd(s,v).real() + yu(w,r).real() * yd(s,v).imag()) * WC5(Cquqd1I, v, w, p, t)
					+(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd1R, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd1I, v, w, r, s)
					)

					-(0.5 / NC2)*(
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7R(Cqd8R, v, r, w, t)
					-(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7I(Cqd8I, v, r, w, t)
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7R(Cqd8R, p, v, s, w)
					-(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7I(Cqd8I, p, v, s, w)
					+(yu(w,r).real() * yd(s,v).real() - yu(w,r).imag() * yd(s,v).imag()) * WC5(Cquqd8R, v, w, p, t)
					-(yu(w,r).imag() * yd(s,v).real() + yu(w,r).real() * yd(s,v).imag()) * WC5(Cquqd8I, v, w, p, t)
					+(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd8R, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd8I, v, w, r, s)
					)

					-(2. / NC)*(
					+(yddag(v,t).real() * yd(s,w).real() - yddag(v,t).imag() * yd(s,w).imag()) * WC6R(Cqq1R, p, v, w, r)
					-(yddag(v,t).imag() * yd(s,w).real() + yddag(v,t).real() * yd(s,w).imag()) * WC6I(Cqq1I, p, v, w, r)
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC6R(CddR, v, t, s, w)
					-(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC6I(CddI, v, t, s, w)
					)

					-(6. / NC)*(
					+(yddag(v,t).real() * yd(s,w).real() - yddag(v,t).imag() * yd(s,w).imag()) * WC6R(Cqq3R, p, v, w, r)
					-(yddag(v,t).imag() * yd(s,w).real() + yddag(v,t).real() * yd(s,w).imag()) * WC6I(Cqq3I, p, v, w, r)
					)

					+ 0.5 * (
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7R(Cqd8R, v, r, w, t)
					-(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7I(Cqd8I, v, r, w, t)
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7R(Cqd8R, p, v, s, w)
					-(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7I(Cqd8I, p, v, s, w)
					)

					+ 0.5 * (
					+(yd(s,w).real() * yu(v,r).real() - yd(s,w).imag() * yu(v,r).imag()) * WC5(Cquqd1R, p, v, w, t)
					-(yd(s,w).imag() * yu(v,r).real() + yd(s,w).real() * yu(v,r).imag()) * WC5(Cquqd1I, p, v, w, t)
					+(yudag(p,v).real() * yddag(w,t).real() - yudag(p,v).imag() * yddag(w,t).imag()) * WC5(Cquqd1R, r, v, w, s)
					+(yudag(p,v).imag() * yddag(w,t).real() + yudag(p,v).real() * yddag(w,t).imag()) * WC5(Cquqd1I, r, v, w, s)
					+(yd(s,v).real() * yu(w,r).real() - yd(s,v).imag() * yu(w,r).imag()) * WC5(Cquqd8R, v, w, p, t)
					-(yd(s,v).imag() * yu(w,r).real() + yd(s,v).real() * yu(w,r).imag()) * WC5(Cquqd8I, v, w, p, t)
					+(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd8R, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd8I, v, w, r, s)
					)

					- 4. * (
					+(yddag(w,t).real() * yd(s,v).real() - yddag(w,t).imag() * yd(s,v).imag()) * WC6R(Cqq1R, p, r, v, w)
					-(yddag(w,t).imag() * yd(s,v).real() + yddag(w,t).real() * yd(s,v).imag()) * WC6I(Cqq1I, p, r, v, w)
					)
					- 2. * (
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC6R(CddR, v, w, s, t)
					-(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC6I(CddI, v, w, s, t)
					)

					-(
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC7R(Cud1R, v, w, s, t)
					-(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC7I(Cud1I, v, w, s, t)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}


	//Cqd1I
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			+(1. / NC)*(yd(s,r).imag() * xid(p,t).real() + yd(s,r).real() * xid(p,t).imag())
			+(1. / NC)*(yddag(p,t).imag() * xid(r,s).real() - yddag(p,t).real() * xid(r,s).imag())
			+(
			+(yudyu(p,r).imag() - yddyd(p,r).imag()) * WC2R(CHdR, s, t)
			+(yudyu(p,r).real() - yddyd(p,r).real()) * WC2I(CHdI, s, t)
			)
			+ 2. * (gd(s,t).imag() * WC2R(CHq1R, p, r) + gd(s,t).real() * WC2I(CHq1I, p, r))
			//RGE 2
			+ FOUR_THIRDS * g12 * YhYq * WC2I(CHdI, s, t) * delta[p][r]
			+ FOUR_THIRDS * g12 * YhYd * WC2I(CHq1I, p, r) * delta[s][t]
			- 12. * YdYq * g12 * WC7I(Cqd1I, p, r, s, t)
			- 3. * ((NC2 - 1.) / NC2) * g32 * WC7I(Cqd8I, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).imag() * WC7R(Cqd1R, v, r, s, t) + gq(p,v).real() * WC7I(Cqd1I, v, r, s, t))
				+(gd(s,v).imag() * WC7R(Cqd1R, p, r, v, t) + gd(s,v).real() * WC7I(Cqd1I, p, r, v, t))
				+(WC7R(Cqd1R, p, v, s, t) * gq(v,r).imag() + WC7I(Cqd1I, p, v, s, t) * gq(v,r).real())
				+(WC7R(Cqd1R, p, r, s, v) * gd(v,t).imag() + WC7I(Cqd1I, p, r, s, v) * gd(v,t).real())
				//RGE 2
				+ EIGHT_THIRDS * g12 * NC * YdYq * delta[s][t]*(
				+WC6I(Cqq1I, p, r, w, w) + WC6I(Cqq1I, w, w, p, r)
				)
				+ FOUR_THIRDS * g12 * YdYq * delta[s][t]*(
				+WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
				)
				+ 4. * g12 * YdYq * delta[s][t]*(
				+WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
				)
				+ EIGHT_THIRDS * g12 * YdYl * WC7I(Clq1I, w, w, p, r) * delta[s][t]
				+ FOUR_THIRDS * g12 * YdYe * WC7I(CqeI, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * YuYd * WC7I(Cqu1I, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g12 * NC * Yd2 * WC7I(Cqd1I, p, r, w, w) * delta[s][t]
				+ EIGHT_THIRDS * g12 * NC * Yq2 * WC7I(Cqd1I, w, w, s, t) * delta[p][r]
				+ EIGHT_THIRDS * g12 * YlYq * WC7I(CldI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * YeYq * WC7I(CedI, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g12 * NC * YuYq * WC7I(Cud1I, w, w, s, t) * delta[p][r]

				+ FOUR_THIRDS * g12 * NC * YdYq * delta[p][r]*(
				+WC6I(CddI, s, t, w, w) + WC6I(CddI, w, w, s, t)
				)
				+ FOUR_THIRDS * g12 * YdYq * delta[p][r]*(
				+WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += +(1. / NC)*(
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7I(Cqd1I, v, r, w, t)
					+(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7R(Cqd1R, v, r, w, t)
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7I(Cqd1I, p, v, s, w)
					+(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7R(Cqd1R, p, v, s, w)
					+(yu(w,r).real() * yd(s,v).real() - yu(w,r).imag() * yd(s,v).imag()) * WC5(Cquqd1I, v, w, p, t)
					+(yu(w,r).imag() * yd(s,v).real() + yu(w,r).real() * yd(s,v).imag()) * WC5(Cquqd1R, v, w, p, t)
					-(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd1I, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd1R, v, w, r, s)
					)

					-(0.5 / NC2)*(
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7I(Cqd8I, v, r, w, t)
					+(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7R(Cqd8R, v, r, w, t)
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7I(Cqd8I, p, v, s, w)
					+(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7R(Cqd8R, p, v, s, w)
					+(yu(w,r).real() * yd(s,v).real() - yu(w,r).imag() * yd(s,v).imag()) * WC5(Cquqd8I, v, w, p, t)
					+(yu(w,r).imag() * yd(s,v).real() + yu(w,r).real() * yd(s,v).imag()) * WC5(Cquqd8R, v, w, p, t)
					-(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd8I, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd8R, v, w, r, s)
					)

					-(2. / NC)*(
					+(yddag(v,t).real() * yd(s,w).real() - yddag(v,t).imag() * yd(s,w).imag()) * WC6I(Cqq1I, p, v, w, r)
					+(yddag(v,t).imag() * yd(s,w).real() + yddag(v,t).real() * yd(s,w).imag()) * WC6R(Cqq1R, p, v, w, r)
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC6I(CddI, v, t, s, w)
					+(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC6R(CddR, v, t, s, w)
					)

					-(6. / NC)*(
					+(yddag(v,t).real() * yd(s,w).real() - yddag(v,t).imag() * yd(s,w).imag()) * WC6I(Cqq3I, p, v, w, r)
					+(yddag(v,t).imag() * yd(s,w).real() + yddag(v,t).real() * yd(s,w).imag()) * WC6R(Cqq3R, p, v, w, r)
					)

					+ 0.5 * (
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7I(Cqd8I, v, r, w, t)
					+(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7R(Cqd8R, v, r, w, t)
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7I(Cqd8I, p, v, s, w)
					+(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7R(Cqd8R, p, v, s, w)
					)

					+ 0.5 * (
					+(yd(s,w).real() * yu(v,r).real() - yd(s,w).imag() * yu(v,r).imag()) * WC5(Cquqd1I, p, v, w, t)
					+(yd(s,w).imag() * yu(v,r).real() + yd(s,w).real() * yu(v,r).imag()) * WC5(Cquqd1R, p, v, w, t)
					-(yudag(p,v).real() * yddag(w,t).real() - yudag(p,v).imag() * yddag(w,t).imag()) * WC5(Cquqd1I, r, v, w, s)
					+(yudag(p,v).imag() * yddag(w,t).real() + yudag(p,v).real() * yddag(w,t).imag()) * WC5(Cquqd1R, r, v, w, s)
					+(yd(s,v).real() * yu(w,r).real() - yd(s,v).imag() * yu(w,r).imag()) * WC5(Cquqd8I, v, w, p, t)
					+(yd(s,v).imag() * yu(w,r).real() + yd(s,v).real() * yu(w,r).imag()) * WC5(Cquqd8R, v, w, p, t)
					-(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd8I, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd8R, v, w, r, s)
					)

					- 4. * (
					+(yddag(w,t).real() * yd(s,v).real() - yddag(w,t).imag() * yd(s,v).imag()) * WC6I(Cqq1I, p, r, v, w)
					+(yddag(w,t).imag() * yd(s,v).real() + yddag(w,t).real() * yd(s,v).imag()) * WC6R(Cqq1R, p, r, v, w)
					)
					- 2. * (
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC6I(CddI, v, w, s, t)
					+(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC6R(CddR, v, w, s, t)
					)

					-(
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC7I(Cud1I, v, w, s, t)
					+(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC7R(Cud1R, v, w, s, t)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}




	//Cqd8R
	for (d = 0; d < DWC7R; d++) {
		p = WC7R_indices[d][0];
		r = WC7R_indices[d][1];
		s = WC7R_indices[d][2];
		t = WC7R_indices[d][3];
		f[c] =
			+2. * (yd(s,r).real() * xid(p,t).real() - yd(s,r).imag() * xid(p,t).imag())
			+ 2. * (yddag(p,t).real() * xid(r,s).real() + yddag(p,t).imag() * xid(r,s).imag())
			//RGE 2
			-(12. * YdYq * g12 + 6. * (NC - (2. / NC)) * g32) * WC7R(Cqd8R, p, r, s, t)
			- 12. * g32 * WC7R(Cqd1R, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).real() * WC7R(Cqd8R, v, r, s, t) - gq(p,v).imag() * WC7I(Cqd8I, v, r, s, t))
				+(gd(s,v).real() * WC7R(Cqd8R, p, r, v, t) - gd(s,v).imag() * WC7I(Cqd8I, p, r, v, t))
				+(WC7R(Cqd8R, p, v, s, t) * gq(v,r).real() - WC7I(Cqd8I, p, v, s, t) * gq(v,r).imag())
				+(WC7R(Cqd8R, p, r, s, v) * gd(v,t).real() - WC7I(Cqd8I, p, r, s, v) * gd(v,t).imag())
				//RGE 2
				+ FOUR_THIRDS * g32 * delta[s][t]*(
				+WC6R(Cqq1R, p, w, w, r) + WC6R(Cqq1R, w, r, p, w)
				)
				+ 4. * g32 * delta[s][t]*(
				+WC6R(Cqq3R, p, w, w, r) + WC6R(Cqq3R, w, r, p, w)
				)
				+ TWO_THIRDS * g32 * WC7R(Cqu8R, p, r, w, w) * delta[s][t]
				+ TWO_THIRDS * g32 * WC7R(Cqd8R, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g32 * WC7R(Cqd8R, w, w, s, t) * delta[p][r]
				+ TWO_THIRDS * g32 * WC7R(Cud8R, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g32 * delta[p][r]*(
				+WC6R(CddR, s, w, w, t) + WC6R(CddR, w, t, s, w)
				)
				; //RGE 3
			for (w = 0; w < NG; w++) {
				f[c] += -(1. / NC)*(
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7R(Cqd8R, v, r, w, t)
					-(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7I(Cqd8I, v, r, w, t)
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7R(Cqd8R, p, v, s, w)
					-(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7I(Cqd8I, p, v, s, w)
					+(yu(w,r).real() * yd(s,v).real() - yu(w,r).imag() * yd(s,v).imag()) * WC5(Cquqd8R, v, w, p, t)
					-(yu(w,r).imag() * yd(s,v).real() + yu(w,r).real() * yd(s,v).imag()) * WC5(Cquqd8I, v, w, p, t)
					+(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd8R, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd8I, v, w, r, s)
					)
					+ 2. * (
					+(
					+(yu(w,r).real() * yd(s,v).real() - yu(w,r).imag() * yd(s,v).imag()) * WC5(Cquqd1R, v, w, p, t)
					-(yu(w,r).imag() * yd(s,v).real() + yu(w,r).real() * yd(s,v).imag()) * WC5(Cquqd1I, v, w, p, t)
					+(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd1R, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd1I, v, w, r, s)
					)
					+ 0.25 * (
					+(yu(v,r).real() * yd(s,w).real() - yu(v,r).imag() * yd(s,w).imag()) * WC5(Cquqd8R, p, v, w, t)
					-(yu(v,r).imag() * yd(s,w).real() + yu(v,r).real() * yd(s,w).imag()) * WC5(Cquqd8I, p, v, w, t)
					+(yudag(p,v).real() * yddag(w,t).real() - yudag(p,v).imag() * yddag(w,t).imag()) * WC5(Cquqd8R, r, v, w, s)
					+(yudag(p,v).imag() * yddag(w,t).real() + yudag(p,v).real() * yddag(w,t).imag()) * WC5(Cquqd8I, r, v, w, s)
					)
					)

					- 2. * (
					+2. * (
					+(yddag(v,t).real() * yd(s,w).real() - yddag(v,t).imag() * yd(s,w).imag()) * WC6R(Cqq1R, p, v, w, r)
					-(yddag(v,t).imag() * yd(s,w).real() + yddag(v,t).real() * yd(s,w).imag()) * WC6I(Cqq1I, p, v, w, r)
					)
					-(
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7R(Cqd1R, v, r, w, t)
					-(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7I(Cqd1I, v, r, w, t)
					)
					-(
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7R(Cqd1R, p, v, s, w)
					-(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7I(Cqd1I, p, v, s, w)
					)
					)
					- 4. * (
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC6R(CddR, v, t, s, w)
					-(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC6I(CddI, v, t, s, w)
					)
					- 12. * (
					+(yddag(v,t).real() * yd(s,w).real() - yddag(v,t).imag() * yd(s,w).imag()) * WC6R(Cqq3R, p, v, w, r)
					-(yddag(v,t).imag() * yd(s,w).real() + yddag(v,t).real() * yd(s,w).imag()) * WC6I(Cqq3I, p, v, w, r)
					)
					-(
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC7R(Cud8R, v, w, s, t)
					-(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC7I(Cud8I, v, w, s, t)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}
	//Cqd8I
	for (d = 0; d < DWC7I; d++) {
		p = WC7I_indices[d][0];
		r = WC7I_indices[d][1];
		s = WC7I_indices[d][2];
		t = WC7I_indices[d][3];
		f[c] =
			+2. * (yd(s,r).imag() * xid(p,t).real() + yd(s,r).real() * xid(p,t).imag())
			+ 2. * (yddag(p,t).imag() * xid(r,s).real() - yddag(p,t).real() * xid(r,s).imag())
			//RGE 2
			-(12. * YdYq * g12 + 6. * (NC - (2. / NC)) * g32) * WC7I(Cqd8I, p, r, s, t)
			- 12. * g32 * WC7I(Cqd1I, p, r, s, t)
			; //RGE 3
		for (v = 0; v < NG; v++) {
			w = v; //Necessary since RGE 3 uses w as summed index
			f[c] +=
				+(gq(p,v).imag() * WC7R(Cqd8R, v, r, s, t) + gq(p,v).real() * WC7I(Cqd8I, v, r, s, t))
				+(gd(s,v).imag() * WC7R(Cqd8R, p, r, v, t) + gd(s,v).real() * WC7I(Cqd8I, p, r, v, t))
				+(WC7R(Cqd8R, p, v, s, t) * gq(v,r).imag() + WC7I(Cqd8I, p, v, s, t) * gq(v,r).real())
				+(WC7R(Cqd8R, p, r, s, v) * gd(v,t).imag() + WC7I(Cqd8I, p, r, s, v) * gd(v,t).real())
				//RGE 2
				+ FOUR_THIRDS * g32 * delta[s][t]*(
				+WC6I(Cqq1I, p, w, w, r) + WC6I(Cqq1I, w, r, p, w)
				)
				+ 4. * g32 * delta[s][t]*(
				+WC6I(Cqq3I, p, w, w, r) + WC6I(Cqq3I, w, r, p, w)
				)
				+ TWO_THIRDS * g32 * WC7I(Cqu8I, p, r, w, w) * delta[s][t]
				+ TWO_THIRDS * g32 * WC7I(Cqd8I, p, r, w, w) * delta[s][t]
				+ FOUR_THIRDS * g32 * WC7I(Cqd8I, w, w, s, t) * delta[p][r]
				+ TWO_THIRDS * g32 * WC7I(Cud8I, w, w, s, t) * delta[p][r]
				+ FOUR_THIRDS * g32 * delta[p][r]*(
				+WC6I(CddI, s, w, w, t) + WC6I(CddI, w, t, s, w)
				)
				;
			for (w = 0; w < NG; w++) {
				f[c] += -(1. / NC)*(
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7I(Cqd8I, v, r, w, t)
					+(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7R(Cqd8R, v, r, w, t)
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7I(Cqd8I, p, v, s, w)
					+(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7R(Cqd8R, p, v, s, w)
					+(yu(w,r).real() * yd(s,v).real() - yu(w,r).imag() * yd(s,v).imag()) * WC5(Cquqd8I, v, w, p, t)
					+(yu(w,r).imag() * yd(s,v).real() + yu(w,r).real() * yd(s,v).imag()) * WC5(Cquqd8R, v, w, p, t)
					-(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd8I, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd8R, v, w, r, s)
					)
					+ 2. * (
					+(
					+(yu(w,r).real() * yd(s,v).real() - yu(w,r).imag() * yd(s,v).imag()) * WC5(Cquqd1I, v, w, p, t)
					+(yu(w,r).imag() * yd(s,v).real() + yu(w,r).real() * yd(s,v).imag()) * WC5(Cquqd1R, v, w, p, t)
					-(yudag(p,w).real() * yddag(v,t).real() - yudag(p,w).imag() * yddag(v,t).imag()) * WC5(Cquqd1I, v, w, r, s)
					+(yudag(p,w).imag() * yddag(v,t).real() + yudag(p,w).real() * yddag(v,t).imag()) * WC5(Cquqd1R, v, w, r, s)
					)
					+ 0.25 * (
					+(yu(v,r).real() * yd(s,w).real() - yu(v,r).imag() * yd(s,w).imag()) * WC5(Cquqd8I, p, v, w, t)
					+(yu(v,r).imag() * yd(s,w).real() + yu(v,r).real() * yd(s,w).imag()) * WC5(Cquqd8R, p, v, w, t)
					-(yudag(p,v).real() * yddag(w,t).real() - yudag(p,v).imag() * yddag(w,t).imag()) * WC5(Cquqd8I, r, v, w, s)
					+(yudag(p,v).imag() * yddag(w,t).real() + yudag(p,v).real() * yddag(w,t).imag()) * WC5(Cquqd8R, r, v, w, s)
					)
					)

					- 2. * (
					+2. * (
					+(yddag(v,t).real() * yd(s,w).real() - yddag(v,t).imag() * yd(s,w).imag()) * WC6I(Cqq1I, p, v, w, r)
					+(yddag(v,t).imag() * yd(s,w).real() + yddag(v,t).real() * yd(s,w).imag()) * WC6R(Cqq1R, p, v, w, r)
					)
					-(
					+(yddag(p,w).real() * yd(s,v).real() - yddag(p,w).imag() * yd(s,v).imag()) * WC7I(Cqd1I, v, r, w, t)
					+(yddag(p,w).imag() * yd(s,v).real() + yddag(p,w).real() * yd(s,v).imag()) * WC7R(Cqd1R, v, r, w, t)
					)
					-(
					+(yddag(v,t).real() * yd(w,r).real() - yddag(v,t).imag() * yd(w,r).imag()) * WC7I(Cqd1I, p, v, s, w)
					+(yddag(v,t).imag() * yd(w,r).real() + yddag(v,t).real() * yd(w,r).imag()) * WC7R(Cqd1R, p, v, s, w)
					)
					)

					- 4. * (
					+(yddag(p,v).real() * yd(w,r).real() - yddag(p,v).imag() * yd(w,r).imag()) * WC6I(CddI, v, t, s, w)
					+(yddag(p,v).imag() * yd(w,r).real() + yddag(p,v).real() * yd(w,r).imag()) * WC6R(CddR, v, t, s, w)
					)
					- 12. * (
					+(yddag(v,t).real() * yd(s,w).real() - yddag(v,t).imag() * yd(s,w).imag()) * WC6I(Cqq3I, p, v, w, r)
					+(yddag(v,t).imag() * yd(s,w).real() + yddag(v,t).real() * yd(s,w).imag()) * WC6R(Cqq3R, p, v, w, r)
					)
					-(
					+(yudag(p,v).real() * yu(w,r).real() - yudag(p,v).imag() * yu(w,r).imag()) * WC7I(Cud8I, v, w, s, t)
					+(yudag(p,v).imag() * yu(w,r).real() + yudag(p,v).real() * yu(w,r).imag()) * WC7R(Cud8R, v, w, s, t)
					)
					; //RGE 2
			}
		}
		f[c] *= loop_factor;
		c++;
	}


	//----------------------RGE SMEFT class 8_LRRL----------------------
	for (p = 0; p < NG; p++) {
		for (r = 0; r < NG; r++) {
			for (s = 0; s < NG; s++) {
				for (t = 0; t < NG; t++) {
					//CledqR
					f[c] = -2. * (
						+(yd(s,t).real() * xie(p,r).real() - yd(s,t).imag() * xie(p,r).imag())
						+(yedag(p,r).real() * xid(t,s).real() + yedag(p,r).imag() * xid(t,s).imag())
						) //RGE 2
						-(
						6. * (Yd * (Yq - Ye) + Ye * (Ye + Yq)) * g12
						+ 3. * (NC - (1. / NC)) * g32
						) * WC5(CledqR, p, r, s, t); //RGE 3
					//CledqI
					f[c + NG * NG * NG * NG] = -2. * (
						+(yd(s,t).imag() * xie(p,r).real() + yd(s,t).real() * xie(p,r).imag())
						+(yedag(p,r).imag() * xid(t,s).real() - yedag(p,r).real() * xid(t,s).imag())
						) //RGE 2
						-(
						6. * (Yd * (Yq - Ye) + Ye * (Ye + Yq)) * g12
						+ 3. * (NC - (1. / NC)) * g32
						) * WC5(CledqI, p, r, s, t); //RGE 3

					//Entries with 1 summed index
					for (v = 0; v < NG; v++) {
						//CledqR
						f[c] +=
							+(gl(p,v).real() * WC5(CledqR, v, r, s, t) - gl(p,v).imag() * WC5(CledqI, v, r, s, t))
							+(gd(s,v).real() * WC5(CledqR, p, r, v, t) - gd(s,v).imag() * WC5(CledqI, p, r, v, t))
							+(WC5(CledqR, p, v, s, t) * ge(v,r).real() - WC5(CledqI, p, v, s, t) * ge(v,r).imag())
							+(WC5(CledqR, p, r, s, v) * gq(v,t).real() - WC5(CledqI, p, r, s, v) * gq(v,t).imag())
							; //RGE 2
						//CledqI
						f[c + NG * NG * NG * NG] +=
							+(gl(p,v).imag() * WC5(CledqR, v, r, s, t) + gl(p,v).real() * WC5(CledqI, v, r, s, t))
							+(gd(s,v).imag() * WC5(CledqR, p, r, v, t) + gd(s,v).real() * WC5(CledqI, p, r, v, t))
							+(WC5(CledqR, p, v, s, t) * ge(v,r).imag() + WC5(CledqI, p, v, s, t) * ge(v,r).real())
							+(WC5(CledqR, p, r, s, v) * gq(v,t).imag() + WC5(CledqI, p, r, s, v) * gq(v,t).real())
							; //RGE 2
						for (w = 0; w < NG; w++) {
							//CledqR
							f[c] += 2. * (
								+(yedag(p,v).real() * yd(w,t).real() - yedag(p,v).imag() * yd(w,t).imag()) * WC7R(CedR, v, r, s, w)
								-(yedag(p,v).imag() * yd(w,t).real() + yedag(p,v).real() * yd(w,t).imag()) * WC7I(CedI, v, r, s, w))
								- 2. * (
								+(yedag(v,r).real() * yd(w,t).real() - yedag(v,r).imag() * yd(w,t).imag()) * WC7R(CldR, p, v, s, w)
								-(yedag(v,r).imag() * yd(w,t).real() + yedag(v,r).real() * yd(w,t).imag()) * WC7I(CldI, p, v, s, w))
								+ 2. * (
								+(yedag(v,r).real() * yd(s,w).real() - yedag(v,r).imag() * yd(s,w).imag()) * WC7R(Clq1R, p, v, w, t)
								-(yedag(v,r).imag() * yd(s,w).real() + yedag(v,r).real() * yd(s,w).imag()) * WC7I(Clq1I, p, v, w, t))
								+ 6. * (
								+(yedag(v,r).real() * yd(s,w).real() - yedag(v,r).imag() * yd(s,w).imag()) * WC7R(Clq3R, p, v, w, t)
								-(yedag(v,r).imag() * yd(s,w).real() + yedag(v,r).real() * yd(s,w).imag()) * WC7I(Clq3I, p, v, w, t))
								- 2. * (
								+(yedag(p,w).real() * yd(s,v).real() - yedag(p,w).imag() * yd(s,v).imag()) * WC7R(CqeR, v, t, w, r)
								-(yedag(p,w).imag() * yd(s,v).real() + yedag(p,w).real() * yd(s,v).imag()) * WC7I(CqeI, v, t, w, r))
								+ 2. * (
								+(yd(s,v).real() * yu(w,t).real() - yd(s,v).imag() * yu(w,t).imag()) * WC5(Clequ1R, p, r, v, w)
								-(yd(s,v).imag() * yu(w,t).real() + yd(s,v).real() * yu(w,t).imag()) * WC5(Clequ1I, p, r, v, w))
								; //RGE 2
							//CledqI
							f[c + NG * NG * NG * NG] += 2. * (
								+(yedag(p,v).real() * yd(w,t).real() - yedag(p,v).imag() * yd(w,t).imag()) * WC7I(CedI, v, r, s, w)
								+(yedag(p,v).imag() * yd(w,t).real() + yedag(p,v).real() * yd(w,t).imag()) * WC7R(CedR, v, r, s, w))
								- 2. * (
								+(yedag(v,r).real() * yd(w,t).real() - yedag(v,r).imag() * yd(w,t).imag()) * WC7I(CldI, p, v, s, w)
								+(yedag(v,r).imag() * yd(w,t).real() + yedag(v,r).real() * yd(w,t).imag()) * WC7R(CldR, p, v, s, w))
								+ 2. * (
								+(yedag(v,r).real() * yd(s,w).real() - yedag(v,r).imag() * yd(s,w).imag()) * WC7I(Clq1I, p, v, w, t)
								+(yedag(v,r).imag() * yd(s,w).real() + yedag(v,r).real() * yd(s,w).imag()) * WC7R(Clq1R, p, v, w, t))
								+ 6. * (
								+(yedag(v,r).real() * yd(s,w).real() - yedag(v,r).imag() * yd(s,w).imag()) * WC7I(Clq3I, p, v, w, t)
								+(yedag(v,r).imag() * yd(s,w).real() + yedag(v,r).real() * yd(s,w).imag()) * WC7R(Clq3R, p, v, w, t))
								- 2. * (
								+(yedag(p,w).real() * yd(s,v).real() - yedag(p,w).imag() * yd(s,v).imag()) * WC7I(CqeI, v, t, w, r)
								+(yedag(p,w).imag() * yd(s,v).real() + yedag(p,w).real() * yd(s,v).imag()) * WC7R(CqeR, v, t, w, r))
								+ 2. * (
								+(yd(s,v).real() * yu(w,t).real() - yd(s,v).imag() * yu(w,t).imag()) * WC5(Clequ1I, p, r, v, w)
								+(yd(s,v).imag() * yu(w,t).real() + yd(s,v).real() * yu(w,t).imag()) * WC5(Clequ1R, p, r, v, w))
								; //RGE 2
						}
					}

					f[c] *= loop_factor;
					f[c + NG * NG * NG * NG] *= loop_factor;
					c++;
				}
			}
		}
	}
	c += NG * NG * NG * NG;


	//----------------------RGE SMEFT class 8_LRLR----------------------
	for (p = 0; p < NG; p++) {
		for (r = 0; r < NG; r++) {
			for (s = 0; s < NG; s++) {
				for (t = 0; t < NG; t++) {
					count = 0;
					//Cquqd1R
					f[c + count * NG * NG * NG * NG] =
						-2. * (yudag(p,r).real() * xid(s,t).real() - yudag(p,r).imag() * xid(s,t).imag())
						- 2. * (yddag(s,t).real() * xiu(p,r).real() - yddag(s,t).imag() * xiu(p,r).imag())
						//RGE 2
						+ 4. * g1 * (Yq + Yu)*
						(WC1(CdBR, s, t) * yudag(p,r).real() - WC1(CdBI, s, t) * yudag(p,r).imag())
						- 6. * g2 *
						(WC1(CdWR, s, t) * yudag(p,r).real() - WC1(CdWI, s, t) * yudag(p,r).imag())
						-(8. / NC) * g1 * (Yq + Yu)*
						(WC1(CdBR, p, t) * yudag(s,r).real() - WC1(CdBI, p, t) * yudag(s,r).imag())
						+(12. / NC) * g2 *
						(WC1(CdWR, p, t) * yudag(s,r).real() - WC1(CdWI, p, t) * yudag(s,r).imag())
						- 8. * ((NC2 - 1.) / NC2) * g3 *
						(WC1(CdGR, p, t) * yudag(s,r).real() - WC1(CdGI, p, t) * yudag(s,r).imag())

						+ 4. * g1 * (Yq + Yd)*
						(WC1(CuBR, p, r) * yddag(s,t).real() - WC1(CuBI, p, r) * yddag(s,t).imag())
						- 6. * g2 *
						(WC1(CuWR, p, r) * yddag(s,t).real() - WC1(CuWI, p, r) * yddag(s,t).imag())
						-(8. / NC) * g1 * (Yq + Yd)*
						(WC1(CuBR, s, r) * yddag(p,t).real() - WC1(CuBI, s, r) * yddag(p,t).imag())
						+(12. / NC) * g2 *
						(WC1(CuWR, s, r) * yddag(p,t).real() - WC1(CuWI, s, r) * yddag(p,t).imag())
						- 8. * ((NC2 - 1.) / NC2) * g3 *
						(WC1(CuGR, s, r) * yddag(p,t).real() - WC1(CuGI, s, r) * yddag(p,t).imag())

						- 0.5 * (
						(3. * Yd2 + 2. * YuYd + 3. * Yu2) * g12 + 3. * g22
						+ 12. * (NC - (1. / NC)) * g32) * WC5(Cquqd1R, p, r, s, t)

						-(1. / NC)*(
						+(3. * Yd2 + 10 * YuYd + 3. * Yu2) * g12 - 3. * g22
						+ 8. * (NC - (1. / NC)) * g32
						) * WC5(Cquqd1R, s, r, p, t)

						- 0.5 * (1. - (1. / NC2))*(
						(3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12
						- 3. * g22 + 4. * (NC - (2. / NC)) * g32) * WC5(Cquqd8R, s, r, p, t)
						+ 2. * (1. - (1. / NC2)) * g32 * WC5(Cquqd8R, p, r, s, t)
						//RGE 3
						;
					count++;
					//Cquqd1I
					f[c + count * NG * NG * NG * NG] =
						-2. * (yudag(p,r).imag() * xid(s,t).real() + yudag(p,r).real() * xid(s,t).imag())
						- 2. * (yddag(s,t).imag() * xiu(p,r).real() + yddag(s,t).real() * xiu(p,r).imag())
						//RGE 2
						+ 4. * g1 * (Yq + Yu)*
						(WC1(CdBR, s, t) * yudag(p,r).imag() + WC1(CdBI, s, t) * yudag(p,r).real())
						- 6. * g2 *
						(WC1(CdWR, s, t) * yudag(p,r).imag() + WC1(CdWI, s, t) * yudag(p,r).real())
						-(8. / NC) * g1 * (Yq + Yu)*
						(WC1(CdBR, p, t) * yudag(s,r).imag() + WC1(CdBI, p, t) * yudag(s,r).real())
						+(12. / NC) * g2 *
						(WC1(CdWR, p, t) * yudag(s,r).imag() + WC1(CdWI, p, t) * yudag(s,r).real())
						- 8. * ((NC2 - 1.) / NC2) * g3 *
						(WC1(CdGR, p, t) * yudag(s,r).imag() + WC1(CdGI, p, t) * yudag(s,r).real())

						+ 4. * g1 * (Yq + Yd)*
						(WC1(CuBR, p, r) * yddag(s,t).imag() + WC1(CuBI, p, r) * yddag(s,t).real())
						- 6. * g2 *
						(WC1(CuWR, p, r) * yddag(s,t).imag() + WC1(CuWI, p, r) * yddag(s,t).real())
						-(8. / NC) * g1 * (Yq + Yd)*
						(WC1(CuBR, s, r) * yddag(p,t).imag() + WC1(CuBI, s, r) * yddag(p,t).real())
						+(12. / NC) * g2 *
						(WC1(CuWR, s, r) * yddag(p,t).imag() + WC1(CuWI, s, r) * yddag(p,t).real())
						- 8. * ((NC2 - 1.) / NC2) * g3 *
						(WC1(CuGR, s, r) * yddag(p,t).imag() + WC1(CuGI, s, r) * yddag(p,t).real())

						- 0.5 * ((3. * Yd2 + 2. * YuYd + 3. * Yu2) * g12 + 3. * g22
						+ 12. * (NC - (1. / NC)) * g32) * WC5(Cquqd1I, p, r, s, t)
						-(1. / NC)*(
						+(3. * Yd2 + 10 * YuYd + 3. * Yu2) * g12 - 3. * g22
						+ 8. * (NC - (1. / NC)) * g32
						) * WC5(Cquqd1I, s, r, p, t)
						- 0.5 * (1 - (1. / NC2))*((3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12
						- 3. * g22 + 4. * (NC - (2. / NC)) * g32) * WC5(Cquqd8I, s, r, p, t)
						+ 2. * (1. - (1. / NC2)) * g32 * WC5(Cquqd8I, p, r, s, t)
						//RGE 3
						;
					count++;


					//Cquqd8R
					f[c + count * NG * NG * NG * NG] = +8. * g3 *
						(WC1(CdGR, s, t) * yudag(p,r).real() - WC1(CdGI, s, t) * yudag(p,r).imag())
						- 16. * g1 * (Yq + Yu)*
						(WC1(CdBR, p, t) * yudag(s,r).real() - WC1(CdBI, p, t) * yudag(s,r).imag())
						+ 24. * g2 *
						(WC1(CdWR, p, t) * yudag(s,r).real() - WC1(CdWI, p, t) * yudag(s,r).imag())
						+(16. / NC) * g3 *
						(WC1(CdGR, p, t) * yudag(s,r).real() - WC1(CdGI, p, t) * yudag(s,r).imag())

						+ 8. * g3 *
						(WC1(CuGR, p, r) * yddag(s,t).real() - WC1(CuGI, p, r) * yddag(s,t).imag())
						- 16. * g1 * (Yq + Yd)*
						(WC1(CuBR, s, r) * yddag(p,t).real() - WC1(CuBI, s, r) * yddag(p,t).imag())
						+ 24. * g2 *
						(WC1(CuWR, s, r) * yddag(p,t).real() - WC1(CuWI, s, r) * yddag(p,t).imag())
						+(16. / NC) * g3 *
						(WC1(CuGR, s, r) * yddag(p,t).real() - WC1(CuGI, s, r) * yddag(p,t).imag())

						+ 8. * g32 * WC5(Cquqd1R, p, r, s, t)
						+(-2. * (3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12 + 6. * g22
						+ (16. / NC) * g32) * WC5(Cquqd1R, s, r, p, t)

						+((-1.5 * Yd2 - YuYd - 1.5 * Yu2) * g12 - 1.5 * g22
						+ 2. * (NC - (1. / NC)) * g32) * WC5(Cquqd8R, p, r, s, t)

						+(1 / NC)*(
						(3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12 - 3. * g22
						+ 4. * (-NC - 2. / NC) * g32) * WC5(Cquqd8R, s, r, p, t)
						; //RGE 3
					count++;
					//Cquqd8I
					f[c + count * NG * NG * NG * NG] = +8. * g3 *
						(WC1(CdGR, s, t) * yudag(p,r).imag() + WC1(CdGI, s, t) * yudag(p,r).real())
						- 16. * g1 * (Yq + Yu)*
						(WC1(CdBR, p, t) * yudag(s,r).imag() + WC1(CdBI, p, t) * yudag(s,r).real())
						+ 24. * g2 *
						(WC1(CdWR, p, t) * yudag(s,r).imag() + WC1(CdWI, p, t) * yudag(s,r).real())
						+(16. / NC) * g3 *
						(WC1(CdGR, p, t) * yudag(s,r).imag() + WC1(CdGI, p, t) * yudag(s,r).real())

						+ 8. * g3 *
						(WC1(CuGR, p, r) * yddag(s,t).imag() + WC1(CuGI, p, r) * yddag(s,t).real())
						- 16. * g1 * (Yq + Yd)*
						(WC1(CuBR, s, r) * yddag(p,t).imag() + WC1(CuBI, s, r) * yddag(p,t).real())
						+ 24. * g2 *
						(WC1(CuWR, s, r) * yddag(p,t).imag() + WC1(CuWI, s, r) * yddag(p,t).real())
						+(16. / NC) * g3 *
						(WC1(CuGR, s, r) * yddag(p,t).imag() + WC1(CuGI, s, r) * yddag(p,t).real())

						+ 8. * g32 * WC5(Cquqd1I, p, r, s, t)
						+(-2. * (3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12 + 6. * g22
						+ (16. / NC) * g32) * WC5(Cquqd1I, s, r, p, t)
						+((-1.5 * Yd2 - YuYd - 1.5 * Yu2) * g12 - 1.5 * g22
						+ 2. * (NC - (1. / NC)) * g32) * WC5(Cquqd8I, p, r, s, t)

						+(1. / NC)*((3. * Yd2 + 10. * YuYd + 3. * Yu2) * g12 - 3. * g22
						+ 4. * (-NC - 2. / NC) * g32) * WC5(Cquqd8I, s, r, p, t)
						; //RGE 3
					count++;
					//Clequ1R
					f[c + count * NG * NG * NG * NG] =
						+2. * (yudag(s,t).real() * xie(p,r).real() - yudag(s,t).imag() * xie(p,r).imag())
						+ 2. * (yedag(p,r).real() * xiu(s,t).real() - yedag(p,r).imag() * xiu(s,t).imag())
						//RGE 2
						-(6. * (Ye2 + Ye * (Yu - Yq) + YuYq) * g12
						+ 3. * (NC - (1. / NC)) * g32) * WC5(Clequ1R, p, r, s, t)
						-(24. * (Yu + Yq)*(2. * Ye - Yq + Yu) * g12 - 18. * g22) * WC5(Clequ3R, p, r, s, t)
						; //RGE 3
					count++;
					//Clequ1I
					f[c + count * NG * NG * NG * NG] =
						+2. * (yudag(s,t).imag() * xie(p,r).real() + yudag(s,t).real() * xie(p,r).imag())
						+ 2. * (yedag(p,r).imag() * xiu(s,t).real() + yedag(p,r).real() * xiu(s,t).imag())
						//RGE 2
						-(6. * (Ye2 + Ye * (Yu - Yq) + YuYq) * g12
						+ 3. * (NC - (1 / NC)) * g32) * WC5(Clequ1I, p, r, s, t)
						-(24. * (Yu + Yq)*(2. * Ye - Yq + Yu) * g12 - 18. * g22) * WC5(Clequ3I, p, r, s, t)
						; //RGE 3
					count++;
					//Clequ3R
					f[c + count * NG * NG * NG * NG] = g1 * (Yq + Yu)*
						(WC1(CeBR, p, r) * yudag(s,t).real() - WC1(CeBI, p, r) * yudag(s,t).imag())
						- 1.5 * g2 *
						(WC1(CuWR, s, t) * yedag(p,r).real() - WC1(CuWI, s, t) * yedag(p,r).imag())
						+ g1 * (Yl + Ye)*
						(WC1(CuBR, s, t) * yedag(p,r).real() - WC1(CuBI, s, t) * yedag(p,r).imag())
						- 1.5 * g2 *
						(WC1(CeWR, p, r) * yudag(s,t).real() - WC1(CeWI, p, r) * yudag(s,t).imag())

						+((2. * (Ye2 - YeYq + YuYe - 2. * Yq2 + 5. * YuYq - 2. * Yu2) * g12 - 3. * g22)
						+(NC - (1 / NC)) * g32) * WC5(Clequ3R, p, r, s, t)

						+ 0.125 * (-4. * (Yq + Yu)*(2. * Ye - Yq + Yu) * g12 + 3. * g22) * WC5(Clequ1R, p, r, s, t)

						; //RGE 3
					count++;
					//Clequ3I
					f[c + count * NG * NG * NG * NG] = g1 * (Yq + Yu)*
						(WC1(CeBR, p, r) * yudag(s,t).imag() + WC1(CeBI, p, r) * yudag(s,t).real())
						- 1.5 * g2 *
						(WC1(CuWR, s, t) * yedag(p,r).imag() + WC1(CuWI, s, t) * yedag(p,r).real())
						+ g1 * (Yl + Ye)*
						(WC1(CuBR, s, t) * yedag(p,r).imag() + WC1(CuBI, s, t) * yedag(p,r).real())
						- 1.5 * g2 *
						(WC1(CeWR, p, r) * yudag(s,t).imag() + WC1(CeWI, p, r) * yudag(s,t).real())

						+((2. * (Ye2 - YeYq + YuYe - 2. * Yq2 + 5. * YuYq - 2. * Yu2) * g12 - 3. * g22)
						+(NC - (1 / NC)) * g32) * WC5(Clequ3I, p, r, s, t)

						+ 0.125 * (-4. * (Yq + Yu)*(2. * Ye - Yq + Yu) * g12 + 3. * g22) * WC5(Clequ1I, p, r, s, t)

						; //RGE 3
					count++;
					//Entries with 1 summed index
					for (v = 0; v < NG; v++) {
						count = 0;
						//Cquqd1R
						f[c + count * NG * NG * NG * NG] +=
							+(gq(p,v).real() * WC5(Cquqd1R, v, r, s, t) - gq(p,v).imag() * WC5(Cquqd1I, v, r, s, t))
							+(gq(s,v).real() * WC5(Cquqd1R, p, r, v, t) - gq(s,v).imag() * WC5(Cquqd1I, p, r, v, t))
							+(WC5(Cquqd1R, p, v, s, t) * gu(v,r).real() - WC5(Cquqd1I, p, v, s, t) * gu(v,r).imag())
							+(WC5(Cquqd1R, p, r, s, v) * gd(v,t).real() - WC5(Cquqd1I, p, r, s, v) * gd(v,t).imag())
							; //RGE 2
						count++;
						//Cquqd1I
						f[c + count * NG * NG * NG * NG] +=
							+(gq(p,v).imag() * WC5(Cquqd1R, v, r, s, t) + gq(p,v).real() * WC5(Cquqd1I, v, r, s, t))
							+(gq(s,v).imag() * WC5(Cquqd1R, p, r, v, t) + gq(s,v).real() * WC5(Cquqd1I, p, r, v, t))
							+(WC5(Cquqd1R, p, v, s, t) * gu(v,r).imag() + WC5(Cquqd1I, p, v, s, t) * gu(v,r).real())
							+(WC5(Cquqd1R, p, r, s, v) * gd(v,t).imag() + WC5(Cquqd1I, p, r, s, v) * gd(v,t).real())
							; //RGE 2
						count++;
						//Cquqd8R
						f[c + count * NG * NG * NG * NG] +=
							+(gq(p,v).real() * WC5(Cquqd8R, v, r, s, t) - gq(p,v).imag() * WC5(Cquqd8I, v, r, s, t))
							+(gq(s,v).real() * WC5(Cquqd8R, p, r, v, t) - gq(s,v).imag() * WC5(Cquqd8I, p, r, v, t))
							+(WC5(Cquqd8R, p, v, s, t) * gu(v,r).real() - WC5(Cquqd8I, p, v, s, t) * gu(v,r).imag())
							+(WC5(Cquqd8R, p, r, s, v) * gd(v,t).real() - WC5(Cquqd8I, p, r, s, v) * gd(v,t).imag())
							; //RGE 2;
						count++;
						//Cquqd8I
						f[c + count * NG * NG * NG * NG] +=
							+(gq(p,v).imag() * WC5(Cquqd8R, v, r, s, t) + gq(p,v).real() * WC5(Cquqd8I, v, r, s, t))
							+(gq(s,v).imag() * WC5(Cquqd8R, p, r, v, t) + gq(s,v).real() * WC5(Cquqd8I, p, r, v, t))
							+(WC5(Cquqd8R, p, v, s, t) * gu(v,r).imag() + WC5(Cquqd8I, p, v, s, t) * gu(v,r).real())
							+(WC5(Cquqd8R, p, r, s, v) * gd(v,t).imag() + WC5(Cquqd8I, p, r, s, v) * gd(v,t).real())
							; //RGE 2
						count++;
						//Clequ1R
						f[c + count * NG * NG * NG * NG] +=
							+(gl(p,v).real() * WC5(Clequ1R, v, r, s, t) - gl(p,v).imag() * WC5(Clequ1I, v, r, s, t))
							+(gq(s,v).real() * WC5(Clequ1R, p, r, v, t) - gq(s,v).imag() * WC5(Clequ1I, p, r, v, t))
							+(WC5(Clequ1R, p, v, s, t) * ge(v,r).real() - WC5(Clequ1I, p, v, s, t) * ge(v,r).imag())
							+(WC5(Clequ1R, p, r, s, v) * gu(v,t).real() - WC5(Clequ1I, p, r, s, v) * gu(v,t).imag())
							; //RGE 2
						count++;
						//Clequ1I
						f[c + count * NG * NG * NG * NG] +=
							+(gl(p,v).imag() * WC5(Clequ1R, v, r, s, t) + gl(p,v).real() * WC5(Clequ1I, v, r, s, t))
							+(gq(s,v).imag() * WC5(Clequ1R, p, r, v, t) + gq(s,v).real() * WC5(Clequ1I, p, r, v, t))
							+(WC5(Clequ1R, p, v, s, t) * ge(v,r).imag() + WC5(Clequ1I, p, v, s, t) * ge(v,r).real())
							+(WC5(Clequ1R, p, r, s, v) * gu(v,t).imag() + WC5(Clequ1I, p, r, s, v) * gu(v,t).real())
							; //RGE 2
						count++;
						//Clequ3R
						f[c + count * NG * NG * NG * NG] +=
							+(gl(p,v).real() * WC5(Clequ3R, v, r, s, t) - gl(p,v).imag() * WC5(Clequ3I, v, r, s, t))
							+(gq(s,v).real() * WC5(Clequ3R, p, r, v, t) - gq(s,v).imag() * WC5(Clequ3I, p, r, v, t))
							+(WC5(Clequ3R, p, v, s, t) * ge(v,r).real() - WC5(Clequ3I, p, v, s, t) * ge(v,r).imag())
							+(WC5(Clequ3R, p, r, s, v) * gu(v,t).real() - WC5(Clequ3I, p, r, s, v) * gu(v,t).imag())
							; //RGE 2
						count++;
						//Clequ1I
						f[c + count * NG * NG * NG * NG] +=
							+(gl(p,v).imag() * WC5(Clequ3R, v, r, s, t) + gl(p,v).real() * WC5(Clequ3I, v, r, s, t))
							+(gq(s,v).imag() * WC5(Clequ3R, p, r, v, t) + gq(s,v).real() * WC5(Clequ3I, p, r, v, t))
							+(WC5(Clequ3R, p, v, s, t) * ge(v,r).imag() + WC5(Clequ3I, p, v, s, t) * ge(v,r).real())
							+(WC5(Clequ3R, p, r, s, v) * gu(v,t).imag() + WC5(Clequ3I, p, r, s, v) * gu(v,t).real())
							; //RGE 2
						count++;
						//Entries with 2 summed indices
						for (w = 0; w < NG; w++) {
							count = 0;
							//Cquqd1R
							f[c + count * NG * NG * NG * NG] += -(2. / NC2)*(
								+(yudag(v,r).real() * yddag(p,w).real() - yudag(v,r).imag() * yddag(p,w).imag()) * WC7R(Cqd8R, s, v, w, t)
								-(yudag(v,r).imag() * yddag(p,w).real() + yudag(v,r).real() * yddag(p,w).imag()) * WC7I(Cqd8I, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7R(Cqu8R, p, v, w, r)
								- (yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7I(Cqu8I, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7R(Cud8R, v, r, w, t)
								- (yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7I(Cud8I, v, r, w, t)
								)

								+(4. / NC)*(
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6R(Cqq1R, s, v, p, w)
								-(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6I(Cqq1I, s, v, p, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6R(Cqq1R, p, v, s, w)
								-(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6I(Cqq1I, p, v, s, w)
								- 3. * ((yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6R(Cqq3R, s, v, p, w)
								-(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6I(Cqq3I, s, v, p, w))
								- 3. * ((yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6R(Cqq3R, p, v, s, w)
								-(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6I(Cqq3I, p, v, s, w))
								)

								+(4. / NC)*(
								+(yddag(p,w).real() * yudag(v,r).real() - yddag(p,w).imag() * yudag(v,r).imag()) * WC7R(Cqd1R, s, v, w, t)
								-(yddag(p,w).imag() * yudag(v,r).real() + yddag(p,w).real() * yudag(v,r).imag()) * WC7I(Cqd1I, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7R(Cqu1R, p, v, w, r)
								-(yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7I(Cqu1I, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7R(Cud1R, v, r, w, t)
								-(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7I(Cud1I, v, r, w, t)
								)

								- 4. * (
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6R(Cqq1R, p, v, s, w)
								-(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6I(Cqq1I, p, v, s, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6R(Cqq1R, s, v, p, w)
								-(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6I(Cqq1I, s, v, p, w)
								)

								+ 12. * (
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6R(Cqq3R, p, v, s, w)
								-(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6I(Cqq3I, p, v, s, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6R(Cqq3R, s, v, p, w)
								-(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6I(Cqq3I, s, v, p, w)
								)

								+ 2. * (
								+(yddag(p,w).real() * yudag(v,r).real() - yddag(p,w).imag() * yudag(v,r).imag()) * WC7R(Cqd8R, s, v, w, t)
								- (yddag(p,w).imag() * yudag(v,r).real() + yddag(p,w).real() * yudag(v,r).imag()) * WC7I(Cqd8I, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7R(Cqu8R, p, v, w, r)
								- (yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7I(Cqu8I, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7R(Cud8R, v, r, w, t)
								- (yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7I(Cud8I, v, r, w, t)
								)
								- 4. * (
								+(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC7R(Cud1R, v, r, w, t)
								-(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC7I(Cud1I, v, r, w, t)
								)
								; //RGE 2
							count++;
							//Cquqd1I
							f[c + count * NG * NG * NG * NG] += -(2. / NC2)*(
								+(yudag(v,r).real() * yddag(p,w).real() - yudag(v,r).imag() * yddag(p,w).imag()) * WC7I(Cqd8I, s, v, w, t)
								+(yudag(v,r).imag() * yddag(p,w).real() + yudag(v,r).real() * yddag(p,w).imag()) * WC7R(Cqd8R, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7I(Cqu8I, p, v, w, r)
								+ (yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7R(Cqu8R, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7I(Cud8I, v, r, w, t)
								+ (yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7R(Cud8R, v, r, w, t)
								)

								+(4. / NC)*(
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6I(Cqq1I, s, v, p, w)
								+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6R(Cqq1R, s, v, p, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6I(Cqq1I, p, v, s, w)
								+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6R(Cqq1R, p, v, s, w)
								- 3. * ((yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6I(Cqq3I, s, v, p, w)
								+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6R(Cqq3R, s, v, p, w))
								- 3. * ((yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6I(Cqq3I, p, v, s, w)
								+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6R(Cqq3R, p, v, s, w))
								)

								+(4. / NC)*(
								+(yddag(p,w).real() * yudag(v,r).real() - yddag(p,w).imag() * yudag(v,r).imag()) * WC7I(Cqd1I, s, v, w, t)
								+(yddag(p,w).imag() * yudag(v,r).real() + yddag(p,w).real() * yudag(v,r).imag()) * WC7R(Cqd1R, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7I(Cqu1I, p, v, w, r)
								+(yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7R(Cqu1R, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7I(Cud1I, v, r, w, t)
								+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7R(Cud1R, v, r, w, t)
								)

								- 4. * (
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6I(Cqq1I, p, v, s, w)
								+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6R(Cqq1R, p, v, s, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6I(Cqq1I, s, v, p, w)
								+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6R(Cqq1R, s, v, p, w)
								)

								+ 12. * (
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6I(Cqq3I, p, v, s, w)
								+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6R(Cqq3R, p, v, s, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6I(Cqq3I, s, v, p, w)
								+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6R(Cqq3R, s, v, p, w)
								)

								+ 2. * (
								+(yddag(p,w).real() * yudag(v,r).real() - yddag(p,w).imag() * yudag(v,r).imag()) * WC7I(Cqd8I, s, v, w, t)
								+ (yddag(p,w).imag() * yudag(v,r).real() + yddag(p,w).real() * yudag(v,r).imag()) * WC7R(Cqd8R, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7I(Cqu8I, p, v, w, r)
								+(yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7R(Cqu8R, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7I(Cud8I, v, r, w, t)
								+ (yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7R(Cud8R, v, r, w, t)
								)

								- 4. * (
								+(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC7I(Cud1I, v, r, w, t)
								+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC7R(Cud1R, v, r, w, t)
								)
								; //RGE 2
							count++;
							//Cquqd8R
							f[c + count * NG * NG * NG * NG] += -(4. / NC) * (
								+(yddag(p,w).real() * yudag(v,r).real() - yddag(p,w).imag() * yudag(v,r).imag()) * WC7R(Cqd8R, s, v, w, t)
								-(yddag(p,w).imag() * yudag(v,r).real() + yddag(p,w).real() * yudag(v,r).imag()) * WC7I(Cqd8I, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7R(Cqu8R, p, v, w, r)
								-(yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7I(Cqu8I, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7R(Cud8R, v, r, w, t)
								-(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7I(Cud8I, v, r, w, t)
								)

								+ 8. * (
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6R(Cqq1R, s, v, p, w)
								-(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6I(Cqq1I, s, v, p, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6R(Cqq1R, p, v, s, w)
								-(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6I(Cqq1I, p, v, s, w)
								)

								- 24. * (
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6R(Cqq3R, s, v, p, w)
								-(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6I(Cqq3I, s, v, p, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6R(Cqq3R, p, v, s, w)
								-(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6I(Cqq3I, p, v, s, w)
								)

								+ 8. * (
								+(yddag(p,w).real() * yudag(v,r).real() - yddag(p,w).imag() * yudag(v,r).imag()) * WC7R(Cqd1R, s, v, w, t)
								-(yddag(p,w).imag() * yudag(v,r).real() + yddag(p,w).real() * yudag(v,r).imag()) * WC7I(Cqd1I, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7R(Cqu1R, p, v, w, r)
								-(yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7I(Cqu1I, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7R(Cud1R, v, r, w, t)
								-(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7I(Cud1I, v, r, w, t)
								)

								- 4. * (
								+(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC7R(Cud8R, v, r, w, t)
								-(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC7I(Cud8I, v, r, w, t)
								)
								; //RGE 2
							count++;
							//Cquqd8I
							f[c + count * NG * NG * NG * NG] += -(4. / NC) * (
								+(yddag(p,w).real() * yudag(v,r).real() - yddag(p,w).imag() * yudag(v,r).imag()) * WC7I(Cqd8I, s, v, w, t)
								+(yddag(p,w).imag() * yudag(v,r).real() + yddag(p,w).real() * yudag(v,r).imag()) * WC7R(Cqd8R, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7I(Cqu8I, p, v, w, r)
								+(yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7R(Cqu8R, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7I(Cud8I, v, r, w, t)
								+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7R(Cud8R, v, r, w, t)
								)

								+ 8. * (
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6I(Cqq1I, s, v, p, w)
								+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6R(Cqq1R, s, v, p, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6I(Cqq1I, p, v, s, w)
								+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6R(Cqq1R, p, v, s, w)
								)

								- 24. * (
								+(yddag(w,t).real() * yudag(v,r).real() - yddag(w,t).imag() * yudag(v,r).imag()) * WC6I(Cqq3I, s, v, p, w)
								+(yddag(w,t).imag() * yudag(v,r).real() + yddag(w,t).real() * yudag(v,r).imag()) * WC6R(Cqq3R, s, v, p, w)
								+(yddag(v,t).real() * yudag(w,r).real() - yddag(v,t).imag() * yudag(w,r).imag()) * WC6I(Cqq3I, p, v, s, w)
								+(yddag(v,t).imag() * yudag(w,r).real() + yddag(v,t).real() * yudag(w,r).imag()) * WC6R(Cqq3R, p, v, s, w)
								)

								+ 8. * (
								+(yddag(p,w).real() * yudag(v,r).real() - yddag(p,w).imag() * yudag(v,r).imag()) * WC7I(Cqd1I, s, v, w, t)
								+(yddag(p,w).imag() * yudag(v,r).real() + yddag(p,w).real() * yudag(v,r).imag()) * WC7R(Cqd1R, s, v, w, t)
								+(yddag(v,t).real() * yudag(s,w).real() - yddag(v,t).imag() * yudag(s,w).imag()) * WC7I(Cqu1I, p, v, w, r)
								+(yddag(v,t).imag() * yudag(s,w).real() + yddag(v,t).real() * yudag(s,w).imag()) * WC7R(Cqu1R, p, v, w, r)
								+(yddag(p,w).real() * yudag(s,v).real() - yddag(p,w).imag() * yudag(s,v).imag()) * WC7I(Cud1I, v, r, w, t)
								+(yddag(p,w).imag() * yudag(s,v).real() + yddag(p,w).real() * yudag(s,v).imag()) * WC7R(Cud1R, v, r, w, t)
								)

								- 4. * (
								+(yddag(s,w).real() * yudag(p,v).real() - yddag(s,w).imag() * yudag(p,v).imag()) * WC7I(Cud8I, v, r, w, t)
								+(yddag(s,w).imag() * yudag(p,v).real() + yddag(s,w).real() * yudag(p,v).imag()) * WC7R(Cud8R, v, r, w, t)
								)
								; //RGE 2
							count++;
							//Clequ1R
							f[c + count * NG * NG * NG * NG] += 2. * (
								+(yddag(s,v).real() * yudag(w,t).real() - yddag(s,v).imag() * yudag(w,t).imag()) * WC5(CledqR, p, r, v, w)
								-(yddag(s,v).imag() * yudag(w,t).real() + yddag(s,v).real() * yudag(w,t).imag()) * WC5(CledqI, p, r, v, w)
								)

								+ 2. * (
								+(yedag(p,v).real() * yudag(s,w).real() - yedag(p,v).imag() * yudag(s,w).imag()) * WC7R(CeuR, v, r, w, t)
								-(yedag(p,v).imag() * yudag(s,w).real() + yedag(p,v).real() * yudag(s,w).imag()) * WC7I(CeuI, v, r, w, t)
								)

								+ 2. * (
								+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC7R(Clq1R, p, v, s, w)
								-(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC7I(Clq1I, p, v, s, w)
								)

								- 6. * (
								+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC7R(Clq3R, p, v, s, w)
								-(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC7I(Clq3I, p, v, s, w)
								)

								- 2. * (
								+(yedag(v,r).real() * yudag(s,w).real() - yedag(v,r).imag() * yudag(s,w).imag()) * WC7R(CluR, p, v, w, t)
								-(yedag(v,r).imag() * yudag(s,w).real() + yedag(v,r).real() * yudag(s,w).imag()) * WC7I(CluI, p, v, w, t)
								)

								- 2. * (
								+(yedag(p,w).real() * yudag(v,t).real() - yedag(p,w).imag() * yudag(v,t).imag()) * WC7R(CqeR, s, v, w, r)
								-(yedag(p,w).imag() * yudag(v,t).real() + yedag(p,w).real() * yudag(v,t).imag()) * WC7I(CqeI, s, v, w, r)
								)
								; //RGE 2
							count++;
							//Clequ1I
							f[c + count * NG * NG * NG * NG] += 2. * (
								+(yddag(s,v).real() * yudag(w,t).real() - yddag(s,v).imag() * yudag(w,t).imag()) * WC5(CledqI, p, r, v, w)
								+(yddag(s,v).imag() * yudag(w,t).real() + yddag(s,v).real() * yudag(w,t).imag()) * WC5(CledqR, p, r, v, w)
								)

								+ 2. * (
								+(yedag(p,v).real() * yudag(s,w).real() - yedag(p,v).imag() * yudag(s,w).imag()) * WC7I(CeuI, v, r, w, t)
								+(yedag(p,v).imag() * yudag(s,w).real() + yedag(p,v).real() * yudag(s,w).imag()) * WC7R(CeuR, v, r, w, t)
								)

								+ 2. * (
								+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC7I(Clq1I, p, v, s, w)
								+(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC7R(Clq1R, p, v, s, w)
								)

								- 6. * (
								+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC7I(Clq3I, p, v, s, w)
								+(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC7R(Clq3R, p, v, s, w)
								)

								- 2. * (
								+(yedag(v,r).real() * yudag(s,w).real() - yedag(v,r).imag() * yudag(s,w).imag()) * WC7I(CluI, p, v, w, t)
								+(yedag(v,r).imag() * yudag(s,w).real() + yedag(v,r).real() * yudag(s,w).imag()) * WC7R(CluR, p, v, w, t)
								)

								- 2. * (
								+(yedag(p,w).real() * yudag(v,t).real() - yedag(p,w).imag() * yudag(v,t).imag()) * WC7I(CqeI, s, v, w, r)
								+(yedag(p,w).imag() * yudag(v,t).real() + yedag(p,w).real() * yudag(v,t).imag()) * WC7R(CqeR, s, v, w, r)
								)
								; //RGE 2	
							count++;
							//Clequ3R
							f[c + count * NG * NG * NG * NG] += -0.5 * (
								+(yudag(s,w).real() * yedag(p,v).real() - yudag(s,w).imag() * yedag(p,v).imag()) * WC7R(CeuR, v, r, w, t)
								-(yudag(s,w).imag() * yedag(p,v).real() + yudag(s,w).real() * yedag(p,v).imag()) * WC7I(CeuI, v, r, w, t)
								)
								- 0.5 * (
								+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC7R(Clq1R, p, v, s, w)
								-(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC7I(Clq1I, p, v, s, w)
								)
								+ 1.5 * (
								+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC7R(Clq3R, p, v, s, w)
								-(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC7I(Clq3I, p, v, s, w)
								)
								- 0.5 * (
								+(yedag(v,r).real() * yudag(s,w).real() - yedag(v,r).imag() * yudag(s,w).imag()) * WC7R(CluR, p, v, w, t)
								-(yedag(v,r).imag() * yudag(s,w).real() + yedag(v,r).real() * yudag(s,w).imag()) * WC7I(CluI, p, v, w, t)
								)
								- 0.5 * (
								+(yedag(p,w).real() * yudag(v,t).real() - yedag(p,w).imag() * yudag(v,t).imag()) * WC7R(CqeR, s, v, w, r)
								-(yedag(p,w).imag() * yudag(v,t).real() + yedag(p,w).real() * yudag(v,t).imag()) * WC7I(CqeI, s, v, w, r)
								)
								; //RGE 2
							count++;
							//Clequ3I
							f[c + count * NG * NG * NG * NG] += -0.5 * (
								+(yudag(s,w).real() * yedag(p,v).real() - yudag(s,w).imag() * yedag(p,v).imag()) * WC7I(CeuI, v, r, w, t)
								+(yudag(s,w).imag() * yedag(p,v).real() + yudag(s,w).real() * yedag(p,v).imag()) * WC7R(CeuR, v, r, w, t)
								)
								- 0.5 * (
								+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC7I(Clq1I, p, v, s, w)
								+(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC7R(Clq1R, p, v, s, w)
								)
								+ 1.5 * (
								+(yedag(v,r).real() * yudag(w,t).real() - yedag(v,r).imag() * yudag(w,t).imag()) * WC7I(Clq3I, p, v, s, w)
								+(yedag(v,r).imag() * yudag(w,t).real() + yedag(v,r).real() * yudag(w,t).imag()) * WC7R(Clq3R, p, v, s, w)
								)
								- 0.5 * (
								+(yedag(v,r).real() * yudag(s,w).real() - yedag(v,r).imag() * yudag(s,w).imag()) * WC7I(CluI, p, v, w, t)
								+(yedag(v,r).imag() * yudag(s,w).real() + yedag(v,r).real() * yudag(s,w).imag()) * WC7R(CluR, p, v, w, t)
								)
								- 0.5 * (
								+(yedag(p,w).real() * yudag(v,t).real() - yedag(p,w).imag() * yudag(v,t).imag()) * WC7I(CqeI, s, v, w, r)
								+(yedag(p,w).imag() * yudag(v,t).real() + yedag(p,w).real() * yudag(v,t).imag()) * WC7R(CqeR, s, v, w, r)
								)
								; //RGE 2
							count++;
						}
					}

					for (count = 0; count < N8_LRLR * 2; count++) {
						f[c + count * NG * NG * NG * NG] *= loop_factor;
					}
					c++;
				}
			}
		}
	}
	c += NG * NG * NG * NG * (2 * N8_LRLR - 1);


//INDICES CHECK
	{
		 
		//In fact decommenting this piece sets all the entries of the ADM to 0
		//so one should verify that all the initial condition
		//stay the same under the evolution.
		//If this happen, the switch from different structures 
		//in main (gslpp::complex numbers, matrices etc) to an array of doubles
		//has been correctly managed  
		/*
				for (i = 0; i < dim; i++) {
					f[i]=0.;
				}
		 */
	}



	return GSL_SUCCESS;
}
