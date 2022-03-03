#ifndef RGESolver_h
#define RGESolver_h
//RGESolver.h

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
//#include "IndependentIndices.h"

#include <unordered_map>
//#include <functional>
#include <boost/function.hpp>
//#include <boost/bind/bind.hpp>

#include "gslpp/src/gslpp.h"

/** 
 * @brief A class that performs renormalization group evolution in the context of the SMEFT
 * @details The class solves the Renormalization Group Equations (RGEs) both numerically and 
 * in the leading-log approximations. 
 * Only operators up to dimension six that preserve  lepton and baryon numbers 
 * are considered. The operator basis is the Warsaw basis, defined in https://arxiv.org/abs/1008.4884. @n 
 * The user must set separately real and imaginary part of each complex parameter.  
 * In tables  @latexonly \ref{SM}, \ref{0F}, \ref{2F} and \ref{4F}@endlatexonly are listed all the 
 * parameters, together with their name (that must be used to correctly invoke getter and setter functions). @n
 * The numerical integration is performed with an adaptive step-size routine 
 * (the Explicit embedded Runge-Kutta-Fehlberg method), using the
 * tools in the GNU Scientific Library. @n 
 * See https://www.gnu.org/software/gsl/doc/html/ode-initval.html for all the details. @n
 * The accuracy level of the numerical integration can be tuned selecting the parameters 
 * @f$\epsilon_{rel}, \epsilon_{abs}@f$ and the integration step using the dedicated getter functions. 
 * 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */


//Latex table with all the names of the coefficients.

/**
\latexonly
 ***************
 \begin{table}[H]
        \centering
        \renewcommand{\arraystretch}{1.5} % Default value: 1
\begin{tabular}{cc}
\begin{tabular}[t]{c|c} 
Parameter     & \multicolumn{1}{c}{Name}      \\ \hline 
$g_1$         & \texttt{g1}      \\
$g_2$         & \texttt{g2}      \\
$g_3$         & \texttt{g3}      \\
$\lambda$         & \texttt{lambda}      \\
$m_h^2$ $[\mathrm{GeV}^2]$  & \texttt{mh2}      \\
$\Re(\mathcal{Y}_u)$         & \texttt{YuR}      \\
$\Im(\mathcal{Y}_u)$         & \texttt{YuI}      \\
$\Re(\mathcal{Y}_d)$         & \texttt{YdR}      \\
$\Im(\mathcal{Y}_d)$         & \texttt{YdI}      \\
 $\Re(\mathcal{Y}_e)$         & \texttt{YeR}      \\
$\Im(\mathcal{Y}_e)$         & \texttt{YeI}      
\end{tabular} &
\begin{tabular}[t]{c|c} 
Parameter     & \multicolumn{1}{c}{Name}      \\ \hline 
$\theta_{12}$         & \texttt{CKM\_theta12}      \\
$\theta_{13}$         & \texttt{CKM\_theta13}      \\
$\theta_{23}$         & \texttt{CKM\_theta23}      \\
%$\delta$         & \texttt{CKM\_delta}      \\
$m_u$ $[\mathrm{GeV}]$  & \texttt{mu}      \\
$m_c$ $[\mathrm{GeV}]$  & \texttt{mc}      \\
$m_t$ $[\mathrm{GeV}]$  & \texttt{mt}      \\
 $m_d$ $[\mathrm{GeV}]$  & \texttt{md}      \\
$m_s$ $[\mathrm{GeV}]$  & \texttt{ms}      \\
$m_b$ $[\mathrm{GeV}]$  & \texttt{mb}      \\
$m_{e}$ $[\mathrm{GeV}]$  & \texttt{mel}      \\
$m_{\mu}$ $[\mathrm{GeV}]$  & \texttt{mmu}      \\
$m_{\tau}$ $[\mathrm{GeV}]$  & \texttt{mtau}      \\
\end{tabular}
\end{tabular} 
\vspace{3 pt}
\caption{Standard Model parameters. The parameters in the left column
 must be set (and accessed) with SetCoefficient (and GetCoefficient) methods.
 The ones in the right column must be set and accessed using 
 other dedicated methods (see the specific documentation for input/output) 
 }
\label{SM}
\end{table}
  
 * 
 * 
 * 
 * 
 * 

 
\begin{table}[H]
        \centering
        \renewcommand{\arraystretch}{1.5} % Default value: 1
\begin{tabular}{cc}
\begin{tabular}[t]{c|c}
 \multicolumn{2}{c}{Classes 1-3} \\ \hline 
Coefficient     & \multicolumn{1}{c}{Name}      \\ \hline 
$C_{G}$         & \texttt{CG}      \\
$C_{\tilde{G}}$ & \texttt{CGtilde} \\
$C_{W}$         & \texttt{CW}      \\
$C_{\tilde{W}}$ & \texttt{CWtilde} \\
$C_H$           & \texttt{CH}      \\
$C_{H \Box} $   & \texttt{CHbox}   \\
$C_{HD}$        & \texttt{CHD}    
\end{tabular} & 
 \begin{tabular}[t]{c|c}
\multicolumn{2}{c}{Class 4} \\ \hline 
Coefficient     & \multicolumn{1}{c}{Name}         \\ \hline 
$C_{HG}$         & \texttt{CHG}      \\
$C_{H\tilde{G}}$ & \texttt{CHGtilde} \\
$C_{HW}$         & \texttt{CHW}      \\
$C_{H\tilde{W}}$ & \texttt{CHWtilde} \\
$C_{HB}$         & \texttt{CHB}      \\
$C_{H\tilde{B}}$ & \texttt{CHBtilde} \\
$C_{HWB}$         & \texttt{CHWB}      \\
$C_{H\tilde{W}B}$ & \texttt{CHWtildeB} 
\end{tabular}
\end{tabular}
\vspace{3 pt}
\caption{Scalar (and real) SMEFT operators. They must be set and accessed 
using SetCoefficient and GetCoefficient.}
\label{0F}
\end{table}
  
 * 
 * 
 * 
 * 
 * 

\begin{table}[H]
        \centering
        \renewcommand{\arraystretch}{1.5} % Default value: 1
\begin{tabular}{ccc}
\begin{tabular}[t]{c|c|c}
 \multicolumn{3}{c}{Class 5} \\ \hline 
Coefficient     & Name   & Symmetry      \\ \hline 
$\Re(C_{eH})$ & \texttt{CeHR} & \multicolumn{1}{c}{WC1}   \\ 
$\Im(C_{eH})$ & \texttt{CeHI}  & \multicolumn{1}{c}{WC1}\\ 
$\Re(C_{uH})$ & \texttt{CuHR}  & \multicolumn{1}{c}{WC1}\\ 
$\Im(C_{uH})$ & \texttt{CuHI}  & \multicolumn{1}{c}{WC1}\\ 
$\Re(C_{dH})$ & \texttt{CdHR}  & \multicolumn{1}{c}{WC1}\\ 
$\Im(C_{dH})$ & \texttt{CdHI}  & \multicolumn{1}{c}{WC1}
\end{tabular} & 
 \begin{tabular}[t]{c|c|c}
\multicolumn{3}{c}{Class 6} \\ \hline 
Coefficient     & Name  & Symmetry \\ \hline 
$\Re(C_{eW})$ & \texttt{CeWR} & WC1 \\ 
$\Im(C_{eW})$ & \texttt{CeWI} & WC1\\ 
$\Re(C_{eB})$ & \texttt{CeBR} & WC1\\ 
$\Im(C_{eB})$ & \texttt{CeBI} & WC1\\ 
$\Re(C_{uG})$ & \texttt{CuGR} & WC1\\ 
$\Im(C_{uG})$ & \texttt{CuGI} & WC1\\ 
$\Re(C_{uW})$ & \texttt{CuWR} & WC1\\ 
$\Im(C_{uW})$ & \texttt{CuWI} & WC1\\ 
$\Re(C_{uB})$ & \texttt{CuBR} & WC1\\ 
$\Im(C_{uB})$ & \texttt{CuBI} & WC1\\ 
$\Re(C_{dG})$ & \texttt{CdGR} & WC1\\ 
$\Im(C_{dG})$ & \texttt{CdGI} & WC1\\ 
$\Re(C_{dW})$ & \texttt{CdWR} & WC1\\ 
$\Im(C_{dW})$ & \texttt{CdWI} & WC1\\ 
$\Re(C_{dB})$ & \texttt{CdBR} & WC1\\ 
$\Im(C_{dB})$ & \texttt{CdBI} & WC1
\end{tabular} &  
 \begin{tabular}[t]{c|c|c}
\multicolumn{3}{c}{Class 7} \\ \hline 
Coefficient     & Name & Symmetry \\ \hline
$\Re(C_{Hl1})$ & \texttt{CHl1R}& WC2R \\ 
$\Im(C_{Hl1})$ & \texttt{CHl1I} & WC2I\\ 
$\Re(C_{Hl3})$ & \texttt{CHl3R} & WC2R\\ 
$\Im(C_{Hl3})$ & \texttt{CHl3I} & WC2I\\ 
$\Re(C_{He})$ & \texttt{CHeR} & WC2R\\ 
$\Im(C_{He})$ & \texttt{CHeI} & WC2I\\ 
$\Re(C_{Hq1})$ & \texttt{CHq1R} & WC2R\\ 
$\Im(C_{Hq1})$ & \texttt{CHq1I} & WC2I\\ 
$\Re(C_{Hq3})$ & \texttt{CHq3R}& WC2R \\ 
$\Im(C_{Hq3})$ & \texttt{CHq3I} & WC2I\\ 
$\Re(C_{Hu})$ & \texttt{CHuR} & WC2R\\ 
$\Im(C_{Hu})$ & \texttt{CHuI} & WC2I\\ 
$\Re(C_{Hd})$ & \texttt{CHdR} & WC2R\\ 
$\Im(C_{Hd})$ & \texttt{CHdI} & WC2I\\
$\Re(C_{Hud})$ & \texttt{CHudR} & WC1\\ 
$\Im(C_{Hud})$ & \texttt{CHudI} & WC1  
\end{tabular}
\end{tabular}
\vspace{3 pt}
\caption{2F SMEFT operators. They must be set and accessed 
using SetCoefficient and GetCoefficient.}
\label{2F}
\end{table}

 * 
 * 
 * 
 *  

 
 \begin{table}[H]
\centering
        \renewcommand{\arraystretch}{1.5} % Default value: 1
 \begin{tabular}[t]{ccc}
\begin{tabular}[t]{c|c|c}
\multicolumn{3}{c}{Class 8 $(\bar{L}L)(\bar{L}L)$} \\ \hline
Coefficient     & Name & Symmetry \\ \hline
$\Re(C_{ll})$ & \texttt{CllR} & WC6R \\
$\Im(C_{ll})$ & \texttt{CllI} & WC6I \\	
$\Re(C_{qq1})$ & \texttt{Cqq1R} & WC6R \\
$\Im(C_{qq1})$ & \texttt{Cqq1I} & WC6I \\	
$\Re(C_{qq3})$ & \texttt{Cqq3R} & WC6R \\
$\Im(C_{qq3})$ & \texttt{Cqq3I} & WC6I \\
$\Re(C_{lq1})$ & \texttt{Clq1R} & WC7R \\
$\Im(C_{lq1})$ & \texttt{Clq1I} & WC7I \\
$\Re(C_{lq3})$ & \texttt{Clq3R} & WC7R \\
$\Im(C_{lq3})$ & \texttt{Clq3I} & WC7I  \\ [7 pt]             
\multicolumn{3}{c}{Class 8 $(\bar{L}R)(\bar{L}R)$} \\ \hline
Coefficient     & Name & Symmetry \\ \hline
$\Re(C_{quqd1})$ & \texttt{Cquqd1R} & WC5 \\
$\Im(C_{quqd1})$ & \texttt{Cquqd1I} & WC5 \\
$\Re(C_{quqd8})$ & \texttt{Cquqd8R} & WC5 \\
$\Im(C_{quqd8})$ & \texttt{Cquqs8I} & WC5 \\
$\Re(C_{lequ1})$ & \texttt{Clequ1R} & WC5 \\
$\Im(C_{lequ1})$ & \texttt{Clequ1I} & WC5 \\	
$\Re(C_{lequ3})$ & \texttt{Clequ3R} & WC5 \\
$\Im(C_{lequ3})$ & \texttt{Clequ3I} & WC5		
\end{tabular} &
\begin{tabular}[t]{c|c|c} 
\multicolumn{3}{c}{Class 8 $(\bar{R}R)(\bar{R}R)$} \\ \hline
Coefficient     & Name & Symmetry \\ \hline
$\Re(C_{ee})$ & \texttt{CeeR} & WC8R \\
$\Im(C_{ee})$ & \texttt{CeeI} & WC8I \\
$\Re(C_{uu})$ & \texttt{CuuR} & WC6R \\
$\Im(C_{uu})$ & \texttt{CuuI} & WC6I \\
$\Re(C_{dd})$ & \texttt{CddR} & WC6R \\
$\Im(C_{dd})$ & \texttt{CddI} & WC6I \\	
$\Re(C_{eu})$ & \texttt{CeuR} & WC7R \\
$\Im(C_{eu})$ & \texttt{CeuI} & WC7I \\
$\Re(C_{ed})$ & \texttt{CedR} & WC7R \\
$\Im(C_{ed})$ & \texttt{CedI} & WC7I \\
$\Re(C_{ud1})$ & \texttt{Cud1R} & WC7R \\
$\Im(C_{ud1})$ & \texttt{Cud1I} & WC7I \\
$\Re(C_{ud8})$ & \texttt{Cud8R} & WC7R \\
$\Im(C_{ud8})$ & \texttt{Cud8I} & WC7I \\ [7 pt]             
\multicolumn{3}{c}{Class 8 $(\bar{L}R)(\bar{R}L)$} \\ \hline
Coefficient     & Name & Symmetry \\ \hline
$\Re(C_{ledq})$ & \texttt{CledqR} & WC5 \\
$\Im(C_{ledq})$ & \texttt{CledqI} & WC5 
\end{tabular} & 
\begin{tabular}[t]{c|c|c} 
\multicolumn{3}{c}{Class 8 $(\bar{L}L)(\bar{R}R)$} \\ \hline
Coefficient     & Name & Symmetry \\ \hline
$\Re(C_{le})$ & \texttt{CleR} & WC7R \\
$\Im(C_{le})$ & \texttt{CleI} & WC7I \\
$\Re(C_{lu})$ & \texttt{CluR} & WC7R \\
$\Im(C_{lu})$ & \texttt{CluI} & WC7I \\
$\Re(C_{ld})$ & \texttt{CldR} & WC7R \\
$\Im(C_{ld})$ & \texttt{CldI} & WC7I \\
$\Re(C_{qe})$ & \texttt{CqeR} & WC7R \\
$\Im(C_{qe})$ & \texttt{CqeI} & WC7I \\
$\Re(C_{qu1})$ & \texttt{Cqu1R} & WC7R \\
$\Im(C_{qu1})$ & \texttt{Cqu1I} & WC7I \\
$\Re(C_{qu8})$ & \texttt{Cqu8R} & WC7R \\
$\Im(C_{qu8})$ & \texttt{Cqu8I} & WC7I \\
$\Re(C_{qd1})$ & \texttt{Cqd1R} & WC7R \\
$\Im(C_{qd1})$ & \texttt{Cqd1I} & WC7I \\
$\Re(C_{qd8})$ & \texttt{Cqd8R} & WC7R \\
$\Im(C_{qd8})$ & \texttt{Cqd8I} & WC7I \\
\end{tabular} 
\end{tabular}

 \caption{4F SMEFT Operators. They must be set and accessed 
using SetCoefficient and GetCoefficient.}
 \label{4F}
 \end{table} 
  
 
 ***************
\endlatexonly
 */

class RGESolver {
public:
    /**
     * @brief The default constructor
     */
    RGESolver();

    /**
     * @brief The default destructor.
     */
    ~RGESolver() {
    }


    /**@name Parameters related to the numeric integration. */

    /**
     * @brief Getter for the relative error used in the numerical integration
     */
    double epsrel() {
        return epsrel_;
    }

    /**
     * @brief Getter for the absolute error used in the numerical integration
     */
    double epsabs() {
        return epsabs_;
    }

    /**
     * @brief Getter for the step used in the numerical integration
     */
    double step() {
        return step_;
    }

    /**
     * @brief Setter for the relative error used in the numerical integration
     */
    void Set_epsrel(double epsrel) {
        epsrel_ = epsrel;
    }

    /**
     * @brief Setter for the absolute error used in the numerical integration
     */
    void Set_epsabs(double epsabs) {
        epsabs_ = epsabs;
    }

    /**
     * @brief Setter  for the step used in the numerical integration
     */
    void Set_step(double step) {
        step_ = step;
    }
    /** @name Evolution */
    /**
     * @brief Performs the RGE evolution
     * @details RGEs are solved with the chosen method from @p muI to @p muF.
     * Currently, the available methods are "Numeric" and "Leading-Log". @n 
     * The evolutor takes as initial values the current values of the parameters, 
     * set with the <tt>SetCoefficient(...)</tt> function. After completing the evolution
     * the values of the parameters are updated and are accessible with the
     * <tt>GetCoefficient(...)</tt> function.
     * @param method resolution method
     * @param muI initial energy scale 
     * @param muF final energy scale 
     */
    void Evolve(std::string method, double muI, double muF);
    //Evolution from muI to muF with method = 1 (Num), 2 (LL)



    /**
     * @brief Generates the initial conditions 
     * for Standard Model's parameters (gauge couplings,
     * Yukawa coupling, quartic coupling and Higgs' boson mass). 
     *
     * @details After evolving the SM parameters up to the scale 
     * mu, CKM parameters and fermion masses are updated with the 
     * new values. 
     * If the flag CKMinput is set to <tt>true</tt> (default), the input 
     * for the Yukawa matrices will be generated from the current value 
     * of the CKM matrix and the masses of the fermions. If 
     * set to false, the current values of the Yukawa matrices 
     * will be used to generate the SM initial conditions at the chosen scale.  
     *  
     * @param mu Scale (in GeV) at which the initial conditions 
     * are generated. If <tt>mu</tt> is different from
     * the scale at which the input is given (SMInputScale), <tt>RGESolver</tt>
     * will use the pure SM RGEs (at one-loop level) to run the parameters 
     * to the scale <tt>mu</tt>.
     * @param basis Flavour basis (
     * <tt>"UP"<tt> or <tt>"DOWN"</tt>)
     * @param method Method used by <tt>RGESolver</tt>
     * to run the SM parameters to the scale <tt>mu</tt> 
     * (<tt>"Numeric"</tt> or <tt>"Leading-Log"</tt>) 
     * @param inputCKM If set to <tt>true</tt> (default), the input 
     * for the Yukawa matrices will be generated from the current value 
     * of the CKM matrix and the masses of the fermions.
     */
    void GenerateSMInitialConditions(double mu, std::string basis, std::string method, bool inputCKM = true);

    /**
     * @brief Same as \ref Evolve, but only for the SM parameters. 
     * The user should use this method instead of \ref Evolve when 
     * interested in pure SM running.  
     * @param method
     * @param muI
     * @param muF
     */
    void EvolveSMOnly(std::string method, double muI, double muF);



    /**@name Input/output  
     * 
     * @brief Documentation for the input/output handling.
     * 
     * @details All the SMEFT coefficients are set using the  
     * SetCoefficient methods and 
     * accessed with the GetCoefficient methods. 
     * There exist three different signatures for each method, 
     * depending on the number of flavour indices of the 
     * parameter (0,2,4). @n
     * These two routines must be used also for the SM parameters 
     * @latexonly $g1,g2,g3,\lambda,m_h^2,$ @endlatexonly
     *  @latexonly $\Re(\mathcal{Y}_u),\Im(\mathcal{Y}_u),$ @endlatexonly
     *  @latexonly $\Re(\mathcal{Y}_d),\Im(\mathcal{Y}_d),$ @endlatexonly
     *  @latexonly $\Re(\mathcal{Y}_e),\Im(\mathcal{Y}_e)$ @endlatexonly
     *  (we follow https://arxiv.org/abs/1308.2627 for what concerns
     * the conventions in the Higgs' sector). @n
     * If the user is interested in using the \ref GenerateSMInitialConditions
     * method, the input for the CKM matrix parameters and 
     * the fermion masses must be given with the 
     * methods 
     * SetCKMAngle(std::string name, double val),
     * SetCKMPhase(double val)
     * SetFermionMass(std::string name, double val). @n
     * A complete list of the keys that must be used to 
     * correctly invoke setter/getter methods are given in 
     * tables  @latexonly \ref{SM}, \ref{0F}, \ref{2F} and \ref{4F}@endlatexonly
     */


    /**
     * @brief Compute CKM matrix and the mass of the fermions. 
     * 
     * @details The methods \ref Evolve and \ref EvolveSMOnly 
     * do not updates the value of CKM parameters and fermion masses after the evolution.
     * This process require the diagonalization of the Yukawa matrices and 
     * may slow the evolution. @n
     * If the user is interested in these parameters (accessible with 
     * \ref GetCKMAngle, \ref GetCKMPhase, \ref GetFermionMass)
     * must invoke this method after the evolution.
     */
    void ComputeCKMAndFermionMasses();


    /**
     * @brief Setter function for the mass of the 
     * fermions. 
     * Assignation is allowed only if the inserted 
     * value is not negative.
     * @param val
     */
    void SetFermionMass(std::string name, double val);
    /**
     * @brief Getter function for the mass of the 
     * fermions.
     */
    double GetFermionMass(std::string name);
    /**
     * @brief Setter function for the CKM matrix angles
     * @latexonly $\theta_{12},\theta_{13},\theta_{23}$@endlatexonly. 
     * The assignation is completed only if the inserted
     * angle is @latexonly $\in [0,\frac{\pi}{2}]$@endlatexonly
     * @param val
     */
    void SetCKMAngle(std::string name, double val);


    /**
     * @brief Getter function for the CKM matrix angles 
     * @latexonly $\theta_{12},\theta_{13},\theta_{23}$@endlatexonly. 
     * @return The selected CKM angle. 
     */
    double GetCKMAngle(std::string name);

    /**
     * @brief Setter function for the CKM matrix phase 
     * @latexonly $\delta$@endlatexonly. 
     * The assignation is completed only if 
     * @latexonly $\delta\in (-\pi,\pi]$@endlatexonly
     * @param val
     */
    void SetCKMPhase(double val);


    /**
     * @brief Getter function for the CKM matrix phase 
     * @latexonly $\delta$@endlatexonly. 
     * @return @latexonly $\delta$@endlatexonly. 
     */
    double GetCKMPhase();

    /**
     * @brief Setter method for the scale at which the method 
     * \ref GenerateSMInitialConditions
     *  takes the input values 
     * for SM parameters.
     * @param mu
     */
    void SetSMInputScale(double mu) {
        InputScale_SM = mu;
    };

    /**
     * @brief Getter method for the scale at which the method 
     * \ref GenerateSMInitialConditions 
     * takes the input values 
     * for SM parameters.
     */
    double GetSMInputScale() {
        return InputScale_SM;
    }

    //Setters for 0F,2F,4F

    /**
     * @brief Setter function for scalar/0F parameters (no flavour indices).
     * @details If the parameter name does not match with any of the parameters, 
     * an error message is printed and no assignation is performed.
     * @param name name of the parameter
     * @param val its value
     */
    void SetCoefficient(std::string name, double val);
    /**
     * @brief Setter function for 2F parameters (2 flavour indices).
     * @details If the parameter name does not match with any of the parameters or if at least 
     * one of the inserted indices is outside the [0:2] range, 
     * @param name name of the parameter
     * @param val its value
     * @param i first flavour index
     * @param j second flavour index
     */
    void SetCoefficient(std::string name, double val, int i, int j);
    /**
     * @brief Setter function for 4F parameters (4 flavour indices).
     * @details If the parameter name does not match with any of the parameters or if at least 
     * one of the inserted indices is outside the [0:2] range, 
     * an error message is printed and no assignation is performed.
     * 
     * @param name name of the parameter
     * @param val its value
     * @param i first flavour index
     * @param j second flavour index
     * @param k third flavour index
     * @param l fourth flavour index
     */
    void SetCoefficient(std::string name, double val, int i, int j,
            int k, int l);

    //Setters passing directly the array (for 2F and 4F)
    //void SetCoefficient(std::string name, double C[]);


    //Getters for 0F,2F,4F
    /**
     * @brief Getter function for scalar/0F parameters (no flavour indices).
     * @details If the parameter name does not match with any of the parameters, 
     * an error message is printed and the value 0 is returned.
     * @param name name of the parameter
     * @return the requested parameter (if it exists), otherwise returns 0. 
     */
    double GetCoefficient(std::string name);
    /**
     * @brief Getter function for 2F parameters (2 flavour indices).
     * @details If the parameter name does not match with any of the parameters 
     * or if at least one of the inserted indices is outside the [0:2] range,  
     * an error message is printed and the value 0 is returned. 
     * @param name name of the parameter
     * @param i first flavour index
     * @param j second flavour index
     * @return the requested parameter (if it exists), otherwise returns 0. 
     */
    double GetCoefficient(std::string name, int i, int j);

    /**
     * @brief Getter function for 4F parameters (4 flavour indices).
     * @details If the parameter name does not match with any of the parameters 
     * or if at least one of the inserted indices is outside the [0:2] range,  
     * an error message is printed and the value 0 is returned. 
     * @param name name of the parameter
     * @param i first flavour index
     * @param j second flavour index
     * @param k third flavour index
     * @param l fourth flavour index
     * @return the requested parameter (if it exists), otherwise returns 0. 
     */
    double GetCoefficient(std::string name, int i, int j,
            int k, int l);


    /**
     * @brief Resets all the SMEFT coefficients to 0 and the 
     * SM parameters to their default value.
     * @details 
     */
    void Reset();

    /**
     * @brief Saves the current values of parameters in a file
     * @details Currently, only "SLHA" format is implemented
     * @param filename Name of the output file 
     * @param format Format of the output file
     */
    void SaveOutputFile(std::string filename,
            std::string format);











private:




    /** @name Flavour */
    /**
     * @brief Goes into the chosen basis. No back-rotation
     * for SMEFT coefficients is performed. 
     * @param basis : allowed options are "UP","DOWN"
     */
    void GoToBasisSMOnly(std::string basis);
    //void GoToBasis(std::string basis);

    /**
     * @brief Extracts from the CKM matrix the 4 
     * physical parameters. 
     */
    void ExtractParametersFromCKM();
    /**
     * @brief Starting from the current values of the 
     * masses, Yukawa matrices are generated in the chosen basis. 
     * @param basis
     */
    void FromMassesToYukawas(std::string basis);
    /**
     * @brief Computes the CKM matrix with the 
     * current values of the angles and phase
     */
    void UpdateCKM();
    /**
     * @brief Inserts the initial values of the SM parameters in the array @p x 
     * @details Only used in @p EvolveSMOnly
     */
    /** @name Private functions to perform the evolution*/

    void InitSMOnly();
    /**
     * @brief Saves the evolved values of the SM parameters from @p x  
     * @details Only used in @p EvolveSMOnly
     */
    void UpdateSMOnly();


    /**
     * @brief Sets all the SM parameters (and the SMInputScale)
     * at the default value and all the SMEFT coefficients to 0.
     */
    void SetSMDefaultInput();

    /**
     * @brief Inserts the initial values of the parameters in the array @p x 
     * @details Only used in @p Evolve
     */
    void Init();
    /**
     * @brief Saves the evolved values of the coefficients from @p x 
     * @details Only used in @p Evolve
     */

    void Update();

    /**
     * @brief Computes the beta functions for the SMEFT.
     * @param logmu value of the logarithm of the energy scale at which the beta functions are computed 
     * @param y 1D array in which are stored the current values of the parameters
     * @param f 1D array in which the beta functions for each parameters are saved
     * @param params eventual additional parameters (not used)
     * @return @p GSL_SUCCESS
     */
    static int func(double logmu, const double y[],
            double f[], void* params);

    /**
     * @brief Computes the beta functions for the SM only.
     * @param logmu value of the logarithm of the energy scale at which the beta functions are computed 
     * @param y 1D array in which are stored the current values of the parameters
     * @param f 1D array in which the beta functions for each parameters are saved
     * @param params eventual additional parameters (not used)
     * @return @p GSL_SUCCESS
     */
    static int funcSMOnly(double logmu, const double y[],
            double f[], void* params);

    /**@name GSL Objects */
    ///@{ 

    /**
     * @brief Relative error used in the integrator with its default 
     * value
     * 
     * ().   */
    double epsrel_ = 0.0000000000000001;

    /**
     * @brief Absolute error used in the integrator with its default 
     * value
     */
    double epsabs_ = 0.0001;

    /**
     * @brief Last step used in the integrator 
     */


    double step_;




    gsl_odeiv2_system sys_ = {func, NULL, 3};
    //gsl_odeiv2_driver * d_ =
    //		  gsl_odeiv2_driver_alloc_y_new(&sys_, gsl_odeiv2_step_rk8pd,
    //		  epsrel_, epsabs_, 0.0);


    gsl_odeiv2_step * s_ = gsl_odeiv2_step_alloc(
            gsl_odeiv2_step_rkf45, 2558);
    gsl_odeiv2_control * con_ = gsl_odeiv2_control_standard_new(
            epsabs_, epsrel_, 1, 1);
    gsl_odeiv2_evolve* evo_ = gsl_odeiv2_evolve_alloc(2558);

    /** @brief 1D array for the integration */
    double x[2558];

    /// @}



    //-----------------------------------------------------------------------------------------

    /** @name Independent entries
     * We follow https://arxiv.org/abs/2010.16341 tab. 15, 16 for the symmetry class
     *  of the operators. @n 
     * Notice that operators in classes WC1 and WC5 have no restrictions
     *  for neither real nor imaginary part, each having @f$N_G ^2@f$ (for WC1) 
     * or @f$N_G ^4@f$ (for WC5).  */

    ///@{

    /**
     *@brief Number of independent entries of the real part of operators in symmetry class WC2 */
    static const int DWC2R = 6;
    /**
     *@brief Number of independent entries of the imaginary part of operators in symmetry class WC2 */
    static const int DWC2I = 3;
    /**
     *@brief Number of independent entries of the real part of operators in symmetry class WC6 */
    static const int DWC6R = 27;
    /**
     *@brief Number of independent entries of the imaginary part of operators in symmetry class WC6 */
    static const int DWC6I = 18;
    /**
     *@brief Number of independent entries of the real part of operators in symmetry class WC7 */
    static const int DWC7R = 45;
    /**
     *@brief Number of independent entries of the imaginary part of operators in symmetry class WC7 */
    static const int DWC7I = 36;
    /**
     *@brief Number of independent entries of the real part of operators in symmetry class WC8 */
    static const int DWC8R = 21;
    /**
     *@brief Number of independent entries of the imaginary part of operators in symmetry class WC8 */
    static const int DWC8I = 15;

    ///@}


    /** @name Indices chosen as independent entries
     * the element <tt>WCn_indices[a][b]</tt> must be interpretated as :
     * the b-th index (there are 4 in 4F operators and 2 in 2F) of the a-th 
     * independent entry for the WCn category, where n = 1,2R,2I...
     */
    ///@{
    /** @brief Independent indices for WC2R */
    static const int WC2R_indices[DWC2R][2];
    /** @brief Independent indices for WC2I */
    static const int WC2I_indices[DWC2I][2];
    /** @brief Independent indices for WC6R */
    static const int WC6R_indices[DWC6R][4];
    /** @brief Independent indices for WC6I */
    static const int WC6I_indices[DWC6I][4];
    /** @brief Independent indices for WC7R */
    static const int WC7R_indices[DWC7R][4];
    /** @brief Independent indices for WC7I */
    static const int WC7I_indices[DWC7I][4];
    /** @brief Independent indices for WC8R */
    static const int WC8R_indices[DWC8R][4];
    /** @brief Independent indices for WC8I */
    static const int WC8I_indices[DWC8I][4];


    ///@}







    /**@name Number of operators for each class
     * See https://arxiv.org/abs/1308.2627 tab. 1 for the full list of operators. 
     * There are 8 different classes, depending on the field content. This classification 
     * is different from the symmetry classification WC1, WC2R, ... @n
     * For each class  @p N is the number of operators and 
     * @p E the number of independent entries. When E is not explicitely
     * defined is understood  <tt> E=N </tt> (no flavour structure)
     */

    ///@{


    /**
     * @brief Number of fermion flavours 
     */
    static const int NG = 3;
    /**
     * @brief Dimension of matrices in flavour space
     */
    static const int DF = 9;
    /**
     * @brief Independent entries of a @f$N_G \times N_G@f$ real symmetric matrix
     */
    static const int DFs = (NG*NG + NG) / 2;
    /**
     * @brief Independent entries of a @f$N_G \times N_G@f$ real anti-symmetric matrix 
     */
    static const int DFa = (NG*NG - NG) / 2;

    /** @brief Number of gauge couplings 
     * @ingroup Miscellaneous*/
    static const int Ngauge = 3;
    /** @brief Number of Higgs' sector parameters
     * @ingroup Miscellaneous */
    static const int Nh = 2;
    /** @brief Number of Yukawa matrices 
     * @ingroup Miscellaneous*/
    static const int Nyukawa = 3;
    /** @brief Number of real parameters for each Yukawa matrix
     *  */
    static const int Eyuk = (Nyukawa * 2 * DF);
    static const int N1 = 4;
    static const int N23 = 3;
    static const int N4 = 8;

    static const int N5 = 3;
    static const int E5 = (N5 * 2 * DF);

    static const int N6 = 8;
    static const int E6 = (N6 * 2 * DF);

    static const int N7 = 8;
    /** @brief Number of Hermitian operators in class 7*/
    static const int N7H = 7;
    /** @brief Number of non-Hermitian operators in class 7*/
    static const int N7NH = 1;
    static const int E7 = (N7H*(DWC2R + DWC2I) + N7NH * 2 * DF);

    static const int N8_LLLL = 5;
    static const int E8_LLLL = 2 * (DWC7R + DWC7I) + 3 * (DWC6R + DWC6I);

    static const int N8_RRRR = 7;
    static const int E8_RRRR = 4 * (DWC7R + DWC7I) + 2 * (DWC6R + DWC6I) + 1 * (DWC8R + DWC8I);

    static const int N8_LLRR = 8;
    static const int E8_LLRR = 8 * (DWC7R + DWC7I);

    static const int N8_LRRL = 1;
    static const int E8_LRRL = 2 * NG*NG*NG*NG*N8_LRRL;
    static const int N8_LRLR = 4;
    static const int E8_LRLR = 2 * NG*NG*NG*NG*N8_LRLR;

    ///@}



    /** @name  Fractions 
     * Recurring fractions defined in order to increase efficiency 
     */

    ///@{
    static const double TWO_THIRDS;
    static const double FOUR_THIRDS;
    static const double EIGHT_THIRDS;
    static const double ONE_THIRD;
    static const double ONE_SIXTH;
    static const double TEN_THIRDS;
    ///@}


    /** @name Miscellaneous parameters
     */

    ///@{
    /**
     * @brief Kroenecker delta in flavour space 
     */
    static const double delta[3][3];


    /** @brief Dimension of the system */
    static const int dim = (Ngauge + Nh + Eyuk + N1 + N23 + N4 + E5 + E6 + E7 + E8_LLLL + E8_RRRR + E8_LLRR + E8_LRRL + E8_LRLR);

    /**@brief Number of colors */
    static const double NC;
    /**@brief Number of colors squared */
    static const double NC2;

    /** @brief Leading-order @f$g_1@f$  beta function
     * ( with @f$g_1@f$ normalized as usual and not as in GUT theories)*/
    static const double b01;
    /** @brief Leading-order @f$g_2@f$  beta function */
    static const double b02;
    /** @brief Leading-order @f$g_3@f$  beta function */
    static const double b03;

    //Casimirs
    /** @brief @f$\mathbf{SU(2)}@f$ adjoint Casimir */
    static const double cA2;
    /** @brief @f$\mathbf{SU(3)}@f$ adjoint Casimir */
    static const double cA3;
    /** @brief @f$\mathbf{SU(2)}@f$ fundamental Casimir */
    static const double cF2;
    /** @brief @f$\mathbf{SU(3)}@f$ fundamental Casimir */
    static const double cF3;
    ///@}

    //Hypercharges (and product of hypercharges)
    //

    /** @name Hypercharges (and products of hypercharges)
     * We use @f$Q = T_3 + Y @f$, being @f$T_3@f$ the z-component of the weak isospin 
     * and @f$Q@f$ the electric charge. The other possibility 
     * is @f$Q = T_3 + \frac{Y}{2} @f$ as in https://arxiv.org/abs/1308.2627 
     */
    ///@{
    static const double Yh;
    static const double Yh2;

    static const double Yq;
    static const double Yq2;
    static const double Yl;
    static const double Yl2;

    static const double Yu;
    static const double Yu2;
    static const double Yd;
    static const double Yd2;
    static const double Ye;
    static const double Ye2;

    static const double YhYu;
    static const double YhYd;
    static const double YhYe;
    static const double YhYq;
    static const double YhYl;

    static const double YuYd;
    static const double YuYe;
    static const double YuYq;
    static const double YuYl;

    static const double YdYe;
    static const double YdYq;
    static const double YdYl;

    static const double YeYq;
    static const double YeYl;

    static const double YlYq;
    ///@}




    /** @name Standard Model parameters 
     */

    ///@{
    /** @brief @f$g_1@f$  */
    double g1, /** @brief @f$g_2@f$  */ g2, /**@brief @f$g_3@f$  */ g3,
    /** @brief @f$m_h ^2 @f$ (Higgs boson mass squared) 
      @details See https://arxiv.org/pdf/1308.2627.pdf for the normalization */ mh2,
    /** @brief @f$ \lambda @f$ (Higgs quartic coupling)
     *  @details See https://arxiv.org/pdf/1308.2627.pdf for the normalization  */ lambda;
    double yuR[3][3], yuI[3][3], ydR[3][3],
    ydI[3][3], yeR[3][3], yeI[3][3];
    ///@}
    /**
     * @brief The CKM matrix
     */
    gslpp::matrix<gslpp::complex> CKM = gslpp::matrix<gslpp::complex>(3, 3, 0.);

    /**
     * @brief the scale at which the method 
     * \ref GenerateSMInitialConditions
     * for SM parameters.
     */
    double InputScale_SM; //GeV

    //CKM parameters
    /**
     * @brief @latexonly $\theta_{12}$ @endlatexonly of the CKM matrix (in radians). 
     * The default value is ??? 
     */
    double CKM_theta12;
    double CKM_theta13;
    double CKM_theta23;
    double CKM_delta;

    /**
     * @brief @latexonly $\cos \theta_{12}$ @endlatexonly 
     */
    double c12;
    /**
     * @brief @latexonly $\sin \theta_{12}$ @endlatexonly 
     */
    double s12;
    /**
     * @brief @latexonly $\cos \theta_{13}$ @endlatexonly 
     */
    double c13;
    /**
     * @brief @latexonly $\sin \theta_{13}$ @endlatexonly 
     */
    double s13;
    /**
     * @brief @latexonly $\cos \theta_{23}$ @endlatexonly 
     */
    double c23;
    /**
     * @brief @latexonly $\sin \theta_{23}$ @endlatexonly 
     */
    double s23;


    /**
     * @brief @latexonly $m_u$ @endlatexonly, the mass of up quark in GeV 
     * (default value ?????).
     */
    double mu;
    double mc;
    double mt;

    double md;
    double ms;
    double mb;

    double mel;
    double mmu;
    double mtau;

    /**@name SMEFT dimension-six operators 
     * By default, all SMEFT dimension-six operators' coefficients are set to 0. 
     * See https://arxiv.org/abs/1308.2627 tab. 1 for the full list of operators. 
     * Each member has a class from 1 to 8 depending on its field contents, as well 
     * as a flavour symmetry classification (WC1, WC2R, WC2I...).  
     */

    ///@{

    /**  @brief @f$ C_G@f$ (class 1, scalar) */ double cG = 0.;
    /**  @brief @f$ C_{\tilde{G}}@f$ (class 1, scalar)*/double cGT = 0.;
    /**  @brief @f$ C_W@f$ (class 1, scalar)*/ double cW = 0.;
    /**  @brief @f$ C_{\tilde{W}} @f$ (class 1, scalar)*/double cWT = 0.;

    /** @brief @f$C_H@f$ (class 2, scalar)*/double cH = 0.;

    /** @brief @f$C_{H \Box} @f$ (class 3, scalar)*/ double cHBOX = 0.;
    /** @brief @f$C_{H D} @f$ (class 3, scalar)*/ double cHD = 0.;

    /** @brief @f$C_{H G} @f$ (class 4, scalar)*/ double cHG = 0.;
    /** @brief @f$C_{H \tilde{G}} @f$ (class 4, scalar)*/ double cHGT = 0.;
    /** @brief @f$C_{H W} @f$ (class 4, scalar)*/ double cHW = 0.;
    /** @brief @f$C_{H \tilde{W}} @f$ (class 4, scalar)*/ double cHWT = 0.;
    /** @brief @f$C_{H B} @f$ (class 4, scalar)*/ double cHB = 0.;
    /** @brief @f$C_{H \tilde{B}} @f$ (class 4, scalar)*/ double cHBT = 0.;
    /** @brief @f$C_{H WB} @f$ (class 4, scalar)*/ double cHWB = 0.;
    /** @brief @f$C_{H \tilde{W} B} @f$ (class 4, scalar)*/double cHWBT = 0.;



    /** @brief @f$ \Re \left[ C_{eH } \right]@f$ (class 5, WC1) */
    double ceHR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{eH } \right]@f$ (class 5, WC1) */
    double ceHI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{uH } \right]@f$ (class 5, WC1) */
    double cuHR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{uH } \right]@f$ (class 5, WC1) */
    double cuHI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{dH } \right]@f$ (class 5, WC1) */
    double cdHR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{dH } \right]@f$ (class 5, WC1) */
    double cdHI[3 * 3] = {0.};


    /** @brief @f$ \Re \left[ C_{eW } \right]@f$ (class 6, WC1) */ double ceWR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{eW } \right]@f$ (class 6, WC1) */ double ceWI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{eB } \right]@f$ (class 6, WC1) */ double ceBR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{eB } \right]@f$ (class 6, WC1) */ double ceBI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{uG } \right]@f$ (class 6, WC1) */double cuGR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{uG } \right]@f$ (class 6, WC1) */double cuGI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{uW } \right]@f$ (class 6, WC1) */double cuWR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{uW } \right]@f$ (class 6, WC1) */double cuWI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{uB } \right]@f$ (class 6, WC1) */double cuBR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{uB } \right]@f$ (class 6, WC1) */double cuBI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{dG } \right]@f$ (class 6, WC1) */double cdGR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{dG } \right]@f$ (class 6, WC1) */double cdGI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{dW } \right]@f$ (class 6, WC1) */double cdWR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{dW } \right]@f$ (class 6, WC1) */double cdWI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{dB } \right]@f$ (class 6, WC1) */double cdBR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{dW } \right]@f$ (class 6, WC1) */double cdBI[3 * 3];




    /** @brief @f$ \Re \left[ C_{Hl1} \right]@f$ (class 7, WC2R) */ double cHl1R[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{Hl1} \right]@f$ (class 7, WC2I) */ double cHl1I[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{Hl3} \right]@f$ (class 7, WC2R) */double cHl3R[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{Hl3} \right]@f$ (class 7, WC2I) */double cHl3I[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{He} \right]@f$ (class 7, WC2R) */double cHeR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{He} \right]@f$ (class 7, WC2I) */double cHeI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{Hq1} \right]@f$ (class 7, WC2R) */double cHq1R[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{Hq1} \right]@f$ (class 7, WC2I) */double cHq1I[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{Hq3} \right]@f$ (class 7, WC2R) */double cHq3R[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{Hq3} \right]@f$ (class 7, WC2I) */double cHq3I[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{Hu} \right]@f$ (class 7, WC2R) */double cHuR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{Hu} \right]@f$ (class 7, WC2I) */double cHuI[3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{Hd} \right]@f$ (class 7, WC2R) */double cHdR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{Hd} \right]@f$ (class 7, WC2I) */double cHdI[3 * 3] = {0.};

    /** @brief @f$ \Re \left[ C_{Hud} \right]@f$ (class 7, WC1) */double cHudR[3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{Hud} \right]@f$ (class 7, WC1) */double cHudI[3 * 3] = {0.};




    /** @brief @f$ \Re \left[ C_{ll} \right]@f$ (class 8-[LL][LL], WC6R) */double cllR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{ll} \right]@f$ (class 8-[LL][LL], WC6I) */double cllI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{qq1} \right]@f$ (class 8-[LL][LL], WC6R) */double cqq1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{qq1} \right]@f$ (class 8-[LL][LL], WC6I) */double cqq1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{qq3} \right]@f$ (class 8-[LL][LL], WC6R) */double cqq3R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{qq3} \right]@f$ (class 8-[LL][LL], WC6I) */double cqq3I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{lq1} \right]@f$ (class 8-[LL][LL], WC7R) */double clq1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{lq1} \right]@f$ (class 8-[LL][LL], WC7I) */double clq1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{lq3} \right]@f$ (class 8-[LL][LL], WC7R) */double clq3R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{lq3} \right]@f$ (class 8-[LL][LL], WC7I) */double clq3I[3 * 3 * 3 * 3] = {0.};


    /** @brief @f$ \Re \left[ C_{uu} \right]@f$ (class 8-[RR][RR], WC6R) */double cuuR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{uu} \right]@f$ (class 8-[RR][RR], WC6I) */double cuuI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{dd} \right]@f$ (class 8-[RR][RR], WC6R) */double cddR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{dd} \right]@f$ (class 8-[RR][RR], WC6I) */double cddI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{ee} \right]@f$ (class 8-[RR][RR], WC8R) */double ceeR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{ee} \right]@f$ (class 8-[RR][RR], WC8I) */double ceeI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{eu} \right]@f$ (class 8-[RR][RR], WC7R) */double ceuR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{eu} \right]@f$ (class 8-[RR][RR], WC7I) */double ceuI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{ed} \right]@f$ (class 8-[RR][RR], WC7R) */double cedR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{ed} \right]@f$ (class 8-[RR][RR], WC7I) */ double cedI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{ud1} \right]@f$ (class 8-[RR][RR], WC7R) */double cud1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{ud1} \right]@f$ (class 8-[RR][RR], WC7I) */double cud1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{ud8} \right]@f$ (class 8-[RR][RR], WC7R) */double cud8R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{ud1} \right]@f$ (class 8-[RR][RR], WC7I) */double cud8I[3 * 3 * 3 * 3] = {0.};



    /** @brief @f$ \Re \left[ C_{le} \right]@f$ (class 8-[LL][RR], WC7R) */double cleR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{le} \right]@f$ (class 8-[LL][RR], WC7I) */double cleI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{lu} \right]@f$ (class 8-[LL][RR], WC7R) */double cluR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{lu} \right]@f$ (class 8-[LL][RR], WC7I) */double cluI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{ld} \right]@f$ (class 8-[LL][RR], WC7R) */double cldR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{ld} \right]@f$ (class 8-[LL][RR], WC7I) */double cldI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{qe} \right]@f$ (class 8-[LL][RR], WC7R) */double cqeR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{qe} \right]@f$ (class 8-[LL][RR], WC7I) */double cqeI[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{qu1} \right]@f$ (class 8-[LL][RR], WC7R) */double cqu1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{qu1} \right]@f$ (class 8-[LL][RR], WC7I) */double cqu1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{qu8} \right]@f$ (class 8-[LL][RR], WC7R) */double cqu8R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{qu8} \right]@f$ (class 8-[LL][RR], WC7I) */double cqu8I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{qd1} \right]@f$ (class 8-[LL][RR], WC7R) */double cqd1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{qd1} \right]@f$ (class 8-[LL][RR], WC7I) */ double cqd1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{qd8} \right]@f$ (class 8-[LL][RR], WC7R) */double cqd8R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{qd8} \right]@f$ (class 8-[LL][RR], WC7I) */ double cqd8I[3 * 3 * 3 * 3] = {0.};




    /** @brief @f$ \Re \left[ C_{ledq} \right]@f$ (class 8-[LR][RL], WC5) */ double cledqR[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{ledq} \right]@f$ (class 8-[LR][RL], WC5) */double cledqI[3 * 3 * 3 * 3] = {0.};



    /** @brief @f$ \Re \left[ C_{lequ1} \right]@f$ (class 8-[LR][LR], WC5) */double clequ1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{lequ1} \right]@f$ (class 8-[LR][LR], WC5) */double clequ1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{lequ3} \right]@f$ (class 8-[LR][LR], WC5) */double clequ3R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{lequ3} \right]@f$ (class 8-[LR][LR], WC5) */double clequ3I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{quqd1} \right]@f$ (class 8-[LR][LR], WC5) */double cquqd1R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{quqd1} \right]@f$ (class 8-[LR][LR], WC5) */double cquqd1I[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Re \left[ C_{quqd8} \right]@f$ (class 8-[LR][LR], WC5) */double cquqd8R[3 * 3 * 3 * 3] = {0.};
    /** @brief @f$ \Im \left[ C_{quqd8} \right]@f$ (class 8-[LR][LR], WC5) */double cquqd8I[3 * 3 * 3 * 3] = {0.};
    ///@}

    //-----------------------------------------------------------------------------






    /** @name Maps for I/O
     * Maps that connects the coefficients with the appropriate getter/setter functions
     * (based on their symmetry properties)
     * 
     */

    ///@{
    std::unordered_map<std::string, double*> FermionMasses;
    std::unordered_map<std::string, double*> CKMAngles;
    std::unordered_map<std::string, double*> Operators0F;
    std::unordered_map<std::string, boost::function<void(int, int, double) >> Setter2F;
    std::unordered_map<std::string, boost::function<double(int, int) >> Getter2F;
    std::unordered_map<std::string, boost::function<void(int, int, int, int, double) >> Setter4F;
    std::unordered_map<std::string, boost::function<double(int, int, int, int) >> Getter4F;
    ///@}


    /** @name Output file
     * @brief Stuff for output on file
     */
    ///@{
    /**
     * @brief Prints on file the coefficient @p c
     * 
     * @details This function is used only in the function
     * <tt>SaveOutputFile</tt>. 
     * Currently only "SLHA" printing for WC1,WC2R/I,
     * WC5,WC6R/I,WC7R/I,WC8R/I is implemented. 
     * 
     * @param c coefficient
     * @param name printed name 
     * @param sym symmetry category of the operator
     * @param format chosen format 
     * @param f pointer to file
     */

    static inline void Print(double* c, std::string name,
            std::string sym, std::string format,
            std::ofstream& f);
    ///@}



    /** @name Setters and Getters 
     * Setter and getter function for each symmetry class. 
     * 
     */


    ///@{
    static inline void Yukawa_set(double y[3][3], int i, int j, double val);
    static inline double Yukawa(double y[3][3], int i, int j);

    static inline void WC1_set(double * c, int i, int j, double val);
    static inline double WC1(double * c, int i, int j);

    static inline double WC2R(double * c, int i, int j);
    static inline void WC2R_set(double * c, int i, int j, double val);
    static inline double WC2I(double * c, int i, int j);
    static inline void WC2I_set(double * c, int i, int j, double val);

    static inline double WC3(double * c, int i, int j);
    static inline void WC3_set(double * c, int i, int j, double val);

    static inline void WC5_set(double * c, int i, int j, int k, int l, double val);
    static inline double WC5(double * c, int i, int j, int k, int l);

    static inline double WC6R(double * c, int i, int j, int k, int l);
    static inline void WC6R_set(double * c, int i, int j, int k, int l, double val);
    static inline double WC6I(double * c, int i, int j, int k, int l);
    static inline void WC6I_set(double * c, int i, int j, int k, int l, double val);

    static inline void WC7R_set(double * c, int i, int j, int k, int l, double val);
    static inline double WC7R(double * c, int i, int j, int k, int l);
    static inline double WC7I(double * c, int i, int j, int k, int l);
    static inline void WC7I_set(double * c, int i, int j, int k, int l, double val);

    static inline double WC8R(double * c, int i, int j, int k, int l);
    static inline void WC8R_set(double * c, int i, int j, int k, int l, double val);
    static inline double WC8I(double * c, int i, int j, int k, int l);
    static inline void WC8I_set(double * c, int i, int j, int k, int l, double val);
    ///@}

};


#endif
