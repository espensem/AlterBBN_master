#include "src/include.h"
#include "./ini_files/ini.h"	//For reading .ini file
#include "./ini_files/ini.c"

//------------ For reading input parameters from .ini file--------
typedef struct
{
    double eta;
    double Nnu, dNnu;
    double tau;
    double xinu1, xinu2, xinu3;
    double dd0, ndd;
    double sd0, nsd, Tdend, Tsend;
    double Sigmad0, nSigmad, TSigmaend;
    double mass_wimp;
    int coupling, type_wimp;
    double phiW;
    int vary_phiW;
    double Tinit;
} configuration;


static int handler(void* user, const char* section, const char* name,
                   const char* value)
{
    configuration* pconfig = (configuration*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
      if (MATCH("parameter", "eta")) {
          pconfig->eta = atof(value);
    } else if (MATCH("parameter", "nnu")) {
          pconfig->Nnu = atof(value);
    } else if (MATCH("parameter", "dnnu")) {
          pconfig->dNnu = atof(value);
    } else if (MATCH("parameter", "tau")) {
          pconfig->tau = atof(value);
    } else if (MATCH("parameter", "xi_1")) {
          pconfig->xinu1 = atof(value);
    } else if (MATCH("parameter", "xi_2")) {
          pconfig->xinu2 = atof(value);
    } else if (MATCH("parameter", "xi_3")) {
          pconfig->xinu3 = atof(value);
    } else if (MATCH("parameter", "dd0")) {
          pconfig->dd0 = atof(value);
    } else if (MATCH("parameter", "ndd")) {
          pconfig->ndd = atof(value);
    } else if (MATCH("parameter", "sd0")) {
          pconfig->sd0 = atof(value);
    } else if (MATCH("parameter", "nsd")) {
          pconfig->nsd = atof(value);
    } else if (MATCH("parameter", "Tdend")) {
          pconfig->Tdend = atof(value);
    } else if (MATCH("parameter", "Tsend")) {
          pconfig->Tsend = atof(value);
    } else if (MATCH("parameter", "Sigmad0")) {
          pconfig->Sigmad0 = atof(value);
    } else if (MATCH("parameter", "nSigmad")) {
          pconfig->nSigmad = atof(value);
    } else if (MATCH("parameter", "TSigmaend")) {
          pconfig->TSigmaend = atof(value);
    } else if (MATCH("parameter", "mass_wimp")) {
          pconfig->mass_wimp = atof(value);
    } else if (MATCH("parameter", "type_wimp")) {
          pconfig->type_wimp = atoi(value);
    } else if (MATCH("parameter", "coupling")) {
          pconfig->coupling = atoi(value);
    } else if (MATCH("parameter", "phiW")) {
          pconfig->phiW = atof(value);
    } else if (MATCH("parameter", "vary_phiW")) {
          pconfig->vary_phiW = atoi(value);
    } else if (MATCH("parameter", "T9i")) {
          pconfig->Tinit = atof(value);
    } else {
        return 0;  /* unknown section/name, error */
    }
    return 1;
}
//----------------------------------------------------------------

int main(int argc,char** argv)
{ 
    struct relicparam paramrelic;	// Central Structure in code which contains all values
	double ratioH[NNUC+1],sigma_ratioH[NNUC+1];
	double H2_H,He3_H,Yp,Li7_H,Li6_H,Be7_H;
    double sigma_H2_H,sigma_He3_H,sigma_Yp,sigma_Li7_H,sigma_Li6_H,sigma_Be7_H;
    char cosmoType;
    int cosmo;
    const char *const standard = "standard";
    const char *const darkdens = "darkdens";
    const char *const reheating = "reheating";
    const char *const wimp = "wimp";
    //const char *const combine = "combine";

    /* The type of run may be stated as an input argument. An empty argument will run
     * the default parameter-free SBBN with eta_10=6.10, tau=880.3 and Nnu=3.046. */
    cosmo = 0;          // Default run
    if (argc==2)
    {
        sscanf(argv[1],"%s",&cosmoType);
        if (strcmp(standard, &cosmoType) == 0) cosmo = 1;
        else if (strcmp(darkdens, &cosmoType) == 0) cosmo = 2;
        else if (strcmp(reheating, &cosmoType) == 0) cosmo = 3;
        else if (strcmp(wimp, &cosmoType) == 0) cosmo = 4;
        //else if (strcmp(combine, &cosmoType) == 0) cosmo = 5;
        else
        {
            printf("\t [ERROR]  Wrong argument. Must be 'standard', 'darkdens', 'reheating' or 'wimp'\n");
            exit(1);
        }
    }
    else if (argc>2)
    {
        printf("\t This program takes maximum 1 input parameter: type of cosmology.\n"
               "\t Must be 'standard', 'darkdens', 'reheating' or 'wimp'\n");
        exit(1);
    }
    /*
    char cwd[1024];
       if (getcwd(cwd, sizeof(cwd)) != NULL)
           fprintf(stdout, "Current working dir: %s\n", cwd);
       else
           perror("getcwd() error");
       return 0;
    exit(1);
    */
    // Parsing the input file, storing it in 'config' structure
    configuration config;

    if (ini_parse("./input.ini", handler, &config) < 0) {
        printf("Can't load 'input.ini'\n");
        return 1;
    }

    // For initialization of all variables
    Init_cosmomodel(&paramrelic);

    // If other than default SBBN scenario is chosen, initialize parameters as given through "input.ini".
    if (cosmo == 1)            // standard cosmology
    {
        Init_cosmomodel_param(config.Tinit,config.eta,config.Nnu,config.dNnu,config.tau,config.xinu1,config.xinu2,
                              config.xinu3,&paramrelic);
    }
    else if (cosmo == 2)       // dark density included
    {
        Init_cosmomodel_param(config.Tinit,config.eta,config.Nnu,config.dNnu,config.tau,config.xinu1,config.xinu2,
                              config.xinu3,&paramrelic);
        Init_dark_density(config.dd0,config.ndd,config.Tdend,&paramrelic);
        Init_dark_entropy(config.sd0,config.nsd,config.Tsend,&paramrelic);
    }
    else if (cosmo == 3)      // reheating included
    {
        Init_cosmomodel_param(config.Tinit,config.eta,config.Nnu,config.dNnu,config.tau,config.xinu1,config.xinu2,
                              config.xinu3,&paramrelic);
        Init_dark_density(config.dd0,config.ndd,config.TSigmaend,&paramrelic);
        Init_dark_entropySigmaD(config.Sigmad0,config.nSigmad,config.TSigmaend,&paramrelic);
    }
    else if (cosmo == 4)          // WIMP included
    {
        double gchi, gchi_t;
        int fermion, selfConjugate; // 1/0 for fermion/boson, 1/0 for self-conjugate/non-self-conjugate
        if (config.type_wimp == 1) // Real scalar
        {
            gchi = 1.;
            gchi_t = 1.;
            fermion = 0;
            selfConjugate = 1;
        }
        else if (config.type_wimp == 2) // Complex scalar
        {
            gchi = 2.;
            gchi_t = 2.;
            fermion = 0;
            selfConjugate = 0;
        }
        else if (config.type_wimp == 3) // Majorana fermion
        {
            gchi = 2.;
            gchi_t = 7./4.;
            fermion = 1;
            selfConjugate = 1;
        }
        else if (config.type_wimp == 4) // Dirac fermion
        {
            gchi = 4.;
            gchi_t = 7./2.;
            fermion = 1;
            selfConjugate = 0;
        }
        else
        {
            printf("\t [ERROR] Incorrect input value for parameter 'type_wimp'.");
            exit(1);
        }
        Init_cosmomodel_param(config.Tinit,config.eta,config.Nnu,config.dNnu,config.tau,config.xinu1,config.xinu2,
                              config.xinu3,&paramrelic);
        Init_wimp(config.mass_wimp,gchi,gchi_t,fermion,config.coupling,config.phiW,config.vary_phiW,selfConjugate,
                  &paramrelic);
        /*
        if (cosmo == 5) {
            Init_dark_density(config.dd0,config.ndd,config.Tdend,&paramrelic);
            Init_dark_entropy(config.sd0,config.nsd,config.Tsend,&paramrelic);
            Init_dark_entropySigmaD(config.Sigmad0,config.nSigmad,config.TSigmaend,&paramrelic);
        }
        */
    }
    else {
        // Cosmo = 0. Parameter-free SBBN. Do nothing.
    }

    //OUTPUT
    // Central values
    printf("\n\tYp\t\tH2/H\t\tHe3/H\t\tLi7/H\t\tLi6/H\t\tBe7/H\n");
    nucl(0,paramrelic,ratioH);
    H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
    printf("cent:\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

    // High values
    nucl(1,paramrelic,ratioH);
    H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
    printf("high:\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

    // Low values
    nucl(2,paramrelic,ratioH);
    H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
    printf("low :\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",Yp,H2_H,He3_H,Li7_H,Li6_H,Be7_H);

    // Uncertainties
    if(nucl_witherrors(3,paramrelic,ratioH,sigma_ratioH))
    {
        H2_H=ratioH[3];Yp=ratioH[6];Li7_H=ratioH[8];Be7_H=ratioH[9];He3_H=ratioH[5];Li6_H=ratioH[7];
        sigma_H2_H=sigma_ratioH[3];sigma_Yp=sigma_ratioH[6];sigma_Li7_H=sigma_ratioH[8];sigma_Be7_H=sigma_ratioH[9];
        sigma_He3_H=sigma_ratioH[5];sigma_Li6_H=sigma_ratioH[7];
        printf("+/- :\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",sigma_Yp,sigma_H2_H,sigma_He3_H,sigma_Li7_H,
               sigma_Li6_H,sigma_Be7_H);
        printf("err%%:\t    %.3f\t    %.3f\t    %.3f\t    %.3f\t   %.3f\t    %.3f\n",100*sigma_Yp/Yp,
               100*sigma_H2_H/H2_H,100*sigma_He3_H/He3_H,100*sigma_Li7_H/Li7_H,100*sigma_Li6_H/Li6_H,
               100*sigma_Be7_H/Be7_H);

    }
    else printf("\t [ERROR]  Uncertainty calculation failed\n\n");

    printf("Params: eta=%.2e\tNeff=%.3f\ttau=%.1f\txi1=%.2e\txi2=%.2e\txi3=%.2e\n\n",paramrelic.eta0,
           paramrelic.Nnu+paramrelic.dNnu,paramrelic.life_neutron,paramrelic.xinu1,paramrelic.xinu2,paramrelic.xinu3);

    int compat=bbn_excluded(0,paramrelic);

    if(compat==1) printf("\t [RESULT] Excluded by BBN constraints\n");
    else if(compat==0) printf("\t [RESULT] Compatible with BBN constraints\n");
    else printf("\t [ERROR]  Computation failed\n");

    return 1;
}

