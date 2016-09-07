// All functions and structures are initialized in this file

#include <stdio.h>				
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <complex.h>
#include <string.h>
#include <strings.h>

/*--------------------------------------------------------------------*/

/*#define DEBUG*/

/*--------------------------------------------------------------------*/

#define pi    3.1415926535897932385
#define zeta3 1.2020569031595942855
#define hbar  6.58211889e-25 /* in GeV.s */
#define Gn    6.67428e-8  /* in m^3.g^-1.s^-2 */

/*--------------------------------------------------------------------*/

typedef struct relicparam
/* structure containing the cosmological model parameters */
{
	int dummy;
	int entropy_model;
    double dd0,ndd,Tdend;		// Related to Dark Energy
    double sd0,nsd,Tsend;		// Related to Dark entropy
	double Sigmad0,nSigmad,TSigmaend;
	double nt0,nnt,Tnend;
	double mgravitino;
    double table_eff[276][3];   // Reads values from the SgStar files
    double eta0;                // Initial Baryon to photon ratio
    double Nnu;                 // Number of Neutrinos (e+- included)
    double dNnu;                // Number of extra neutrinos (delta N_nu)
    double life_neutron;		// neutron lifetime
    double xinu1,xinu2,xinu3;	// [e-,neutrino], [muon,neutrino], [tau,neutrino] respectively (degeneracy parameters)
    double mass_wimp, g_chi, g_chi_tilde;    // mass and internal degrees of freedom of wimp
    int wimp_added;             // Switch to activate WIMP calculations
    int beta;                   // type of wimp, 1/0 for fermion/boson
    int SM_coupling_wimp;       // wimp coupling to SM particles. 1 for EM, 2 for neutrino, 3 for neut. and eq. neut.
    int selfConjugate;          // 1/0 for self-conjugate/non-self-conjugate WIMP
    double phiW;                // Parameterized WIMP chemical potential
    double Tinit;               // Initial temperature
    /* TO DO: Include a way to vary the chemical potential! */
    int vary_phiW;

}
relicparam;

/*--------------------------------------------------------------------*/
/* Prototypes */

#define NNUCREAC 100	// Expanded from the original 88 reactions
#define NNUC 26

/* general.c */
double max(double x, double y);
double min(double x, double y);
double complex polylog(int n, int m, double x);
double Li2(double x);
double Li3(double x);
double Li4(double x);
double complex CLi2(double complex x);
double complex CLi3(double complex x);
double Cl2(double x);
double Cl3(double x);
double I0(double x);
double I1(double x);
double K0(double x);
double K1(double x);
double K2(double x);
double K3(double x);
double K4(double x);
double Lbessel(double x);
double Mbessel(double x);
double Nbessel(double x);
double K0exp(double x,double z);
double K1exp(double x,double z);
double K2exp(double x,double z);
double expcor(double x);
double kron(int x, int y);
double rand_gauss(void);
int test_integer(char name[]);
int test_file(char *name);

/* omega.c */
void Init_modeleff(int model_eff, struct relicparam* paramrelic);
double heff(double Temp, struct relicparam paramrelic);
double sgStar(double Temp, struct relicparam paramrelic);
double geff(double Temp, struct relicparam paramrelic);
void Init_cosmomodel(struct relicparam* paramrelic);
void Init_cosmomodel_param(double Tinit, double eta, double Nnu, double dNnu, double life_neutron, double xinu1,
                           double xinu2, double xinu3, struct relicparam* paramrelic);
void Init_dark_density(double dd0, double ndd, double T_end, struct relicparam* paramrelic);
void Init_dark_entropy(double sd0, double nsd, double T_end, struct relicparam* paramrelic);
void Init_dark_entropySigmaD(double Sigmad0, double nSigmad, double T_end, struct relicparam* paramrelic);
void Init_nonthermal(double nt0, double nnt, double T_end, struct relicparam* paramrelic);
void Init_wimp(double mass_wimp, double g_chi, double g_chi_tilde, int beta, int SM_coupling_wimp, double phiW,
               int vary_phiW, int selfConjugate, struct relicparam* paramrelic);
double dark_density(double T, struct relicparam paramrelic);
double dark_entropy(double T, struct relicparam paramrelic);
double dark_entropy_derivative(double T, struct relicparam paramrelic);
double dark_entropy_Sigmad(double T, struct relicparam paramrelic);
double nonthermal(double T, struct relicparam paramrelic);
double neutdens(double Tnu, struct relicparam paramrelic);
double neutdens_deriv(double Tnu, struct relicparam paramrelic);

/* bbnrate.c */
void rate_weak(int err, double f[]);
void rate_pn(int err, struct relicparam paramrelic, double f[], double r[], double T9, double Tnu);
void rate_all(int err, double f[], double T9, relicparam paramrelic);

/* bbn.c */
int linearize(double T9, double reacparam[][10], double f[], double r[], int loop, int inc, int ip, double dt, double y0[], double y[], double dydt[], double H, double rhob);
int nucl(int err, struct relicparam paramrelic, double ratioH[]);
int nucl_failsafe(int err, struct relicparam paramrelic, double ratioH[]);
int nucl_witherrors(int err, struct relicparam paramrelic, double ratioH[], double sigma_ratioH[]);
int bbn_excluded(int err, struct relicparam paramrelic);
