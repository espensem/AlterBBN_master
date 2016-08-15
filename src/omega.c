#include "include.h"			//J Related to Sgstar

/*--------------------------------------------------------------*/

void Init_modeleff(int model_eff, struct relicparam* paramrelic)
/* modifies the model of the QCD equations of state */
{
	int ie,je;
	
	if(model_eff==1)
	{
		const double tableA[276][3]=
		{
#include "sgStar_heff/sgStar_heff_A.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableA[ie][je];
	}
	else if(model_eff==2)
	{
		const double tableB[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableB[ie][je];
	}
	else if(model_eff==3)
	{
	const double tableB2[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B2.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableB2[ie][je];
		}
	else if(model_eff==4)
	{
		const double tableB3[276][3]=
		{
#include "sgStar_heff/sgStar_heff_B3.tab"
  		};	
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableB3[ie][je];
	}
	else if(model_eff==5)
	{
		const double tableC[276][3]=
		{
#include "sgStar_heff/sgStar_heff_C.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableC[ie][je];
	}
	else
	{
		const double tableold[276][3]=
		{
#include "sgStar_heff/sgStar_heff_old.tab"
  		};
		for(ie=0;ie<=275;ie++) for(je=0;je<=2;je++) paramrelic->table_eff[ie][je]=tableold[ie][je];
	}
	return;
}

/*--------------------------------------------------------------*/

double heff(double Temp, struct relicparam paramrelic)
/* computes heff at the temperature Temp */
{
	int ie;
	
	if(Temp>= paramrelic.table_eff[0][0]) return paramrelic.table_eff[0][2];
	
	if(Temp<=0.) return paramrelic.table_eff[275][2];

	ie=1;
	while(Temp<paramrelic.table_eff[ie][0]) ie++;
	
	double heff1,heff2,T1,T2;
	heff1=paramrelic.table_eff[ie][2];
	heff2=paramrelic.table_eff[ie-1][2];
	T1=paramrelic.table_eff[ie][0];
	T2=paramrelic.table_eff[ie-1][0];
			
	return (heff2-heff1)/(T2-T1)*(Temp-T1)+heff1;
}

/*--------------------------------------------------------------*/

double sgStar(double Temp, struct relicparam paramrelic)
/* computes sgStar at the temperature Temp */
{
	int ie;
	
	if(Temp>= paramrelic.table_eff[0][0]) return paramrelic.table_eff[0][1];
	
	if(Temp<=0.) return paramrelic.table_eff[275][1];

	ie=1;
	while(Temp<paramrelic.table_eff[ie][0]) ie++;

	double sgStar1,sgStar2,T1,T2;
	sgStar1=paramrelic.table_eff[ie][1];
	sgStar2=paramrelic.table_eff[ie-1][1];
	T1=paramrelic.table_eff[ie][0];
	T2=paramrelic.table_eff[ie-1][0];
		
	return (sgStar2-sgStar1)/(T2-T1)*(Temp-T1)+sgStar1;
}

/*--------------------------------------------------------------*/

double geff(double Temp, struct relicparam paramrelic)
/* computes geff at the temperature Temp */
{
	double heff0=heff(Temp,paramrelic);
	
	return pow(heff0/sgStar(Temp,paramrelic)*(1.+(heff(Temp*1.001,paramrelic)-heff(Temp*0.999,paramrelic))/0.006/heff0),2.);
}

/*--------------------------------------------------------------*/

void Init_cosmomodel(struct relicparam* paramrelic)
/* initializes the parameters contained in paramrelic */
{
    paramrelic->Tinit=27.;          // Starting at Tnud=2.3 MeV as default
    paramrelic->eta0=6.10e-10;      // Baryon-to-photon ratio (Planck 2015 results XIII)
    paramrelic->Nnu=3.046;          // Number of SM neutrinos, e+- reheating included
    paramrelic->dNnu=0.;            // Number of extra neutrino species (e.g. sterile neutrinos)
    paramrelic->life_neutron=880.3; // Neutron lifetime (PDG2014)
    paramrelic->wimp_added=0;
    paramrelic->vary_phiW=0;
    paramrelic->xinu1=0.;
    paramrelic->xinu2=0.;
    paramrelic->xinu3=0.;
	paramrelic->dd0=paramrelic->ndd=paramrelic->Tdend=0.;
	paramrelic->sd0=paramrelic->nsd=paramrelic->Tsend=0.;
	paramrelic->nt0=paramrelic->nnt=paramrelic->Tnend=0.;
	paramrelic->Sigmad0=paramrelic->nSigmad=paramrelic->TSigmaend=0.;
	Init_modeleff(2,paramrelic);
	return;
}

/*--------------------------------------------------------------*/

void Init_cosmomodel_param(double Tinit, double eta, double Nnu, double dNnu, double life_neutron, double xinu1, double xinu2, double xinu3, struct relicparam* paramrelic)
/* modifies the values of the baryon-to-photon ratio eta, the number of SM neutrinos Nnu, extra neutrino species dNnu
 *  and the neutron lifetime life_neutron */
{
    paramrelic->Tinit=Tinit;
	paramrelic->eta0=eta;
    paramrelic->Nnu=Nnu;
    paramrelic->dNnu=dNnu;
	paramrelic->life_neutron=life_neutron;
    paramrelic->xinu1=xinu1;
    paramrelic->xinu2=xinu2;
    paramrelic->xinu3=xinu3;
	return;
}

/*--------------------------------------------------------------*/

void Init_dark_density(double dd0, double ndd, double T_end, struct relicparam* paramrelic)
/* modifies the parameters of the dark energy density which appears in the Friedmann equation */
{
	paramrelic->dd0=dd0;
	paramrelic->ndd=ndd;
	paramrelic->Tdend=T_end;
	
	return;
}

/*--------------------------------------------------------------*/

void Init_dark_entropy(double sd0, double nsd, double T_end, struct relicparam* paramrelic)
/* modifies the parameters of the dark entropy density which appears in the cosmological equations */
{
	paramrelic->sd0=sd0;
	paramrelic->nsd=nsd;
	paramrelic->Tsend=T_end;
	
	return;
}

/*--------------------------------------------------------------*/

void Init_dark_entropySigmaD(double Sigmad0, double nSigmad, double T_end, struct relicparam* paramrelic)
/* modifies the parameters of the dark entropy generation which appears in the cosmological equations */
{
	paramrelic->Sigmad0=Sigmad0;
	paramrelic->nSigmad=nSigmad;
	paramrelic->TSigmaend=T_end;
	
	return;
}

/*--------------------------------------------------------------*/

void Init_nonthermal(double nt0, double nnt, double T_end, struct relicparam* paramrelic)
/* modifies the parameters of the non-thermal relic particle production which appears in the cosmological equations */
{
	paramrelic->nt0=nt0;
	paramrelic->nnt=nnt;
	paramrelic->Tnend=T_end;
	
	return;
}

/*--------------------------------------------------------------*/

void Init_wimp(double mass_wimp, double g_chi, double g_chi_tilde, int beta, int SM_coupling_wimp, double phiW, int vary_phiW, int selfConjugate, struct relicparam* paramrelic)
/* modifies the parameters of an included light WIMP */
{
    paramrelic->mass_wimp=mass_wimp;
    paramrelic->g_chi=g_chi;
    paramrelic->g_chi_tilde=g_chi_tilde;
    paramrelic->beta=beta;
    paramrelic->SM_coupling_wimp=SM_coupling_wimp;
    paramrelic->wimp_added=1;
    paramrelic->phiW=phiW;
    paramrelic->vary_phiW=vary_phiW;
    paramrelic->selfConjugate=selfConjugate;
    return;
}

/*--------------------------------------------------------------*/

double dark_density(double T, struct relicparam paramrelic)
/* computes the dark energy density at temperature T */
{
	if((paramrelic.dd0==0.)||(T<paramrelic.Tdend)) return 0.;
	
	double rho_rad_1MeV=pi*pi/30.*geff(1.e-3,paramrelic)*1.e-12;
	
	return paramrelic.dd0*rho_rad_1MeV*pow(T/1.e-3,paramrelic.ndd);
}

/*--------------------------------------------------------------*/

double dark_entropy(double T, struct relicparam paramrelic)
/* computes the dark entropy density at temperature T */
{
	if((paramrelic.sd0==0.)&&(paramrelic.Sigmad0==0.)) return 0.;
	
	if((paramrelic.Sigmad0==0.)&&(T<paramrelic.Tsend)) return 0.;
	
	if(paramrelic.Sigmad0==0.)
	{
		double s_rad_1MeV=2.*pi*pi/45.*heff(1.e-3,paramrelic)*1.e-9;
	
		return paramrelic.sd0*s_rad_1MeV*pow(T/1.e-3,paramrelic.nsd);
	}
	else
	{
		double lnT,dlnT,Ttmp;
		int ie,nmax;
		double integ=0.;
		
		nmax=50;
		
		lnT=log(1.e-15);
		
		dlnT=(log(T)-lnT)/nmax;
		
		for(ie=1;ie<nmax;ie++) 
		{
			lnT+=dlnT;
			Ttmp=exp(lnT);
            integ+=sgStar(Ttmp,paramrelic)*dark_entropy_Sigmad(Ttmp,paramrelic)/
                    sqrt(1.+dark_density(Ttmp,paramrelic)/(pi*pi/30.*geff(Ttmp,paramrelic)*pow(Ttmp,4.)))/
                    pow(heff(Ttmp,paramrelic),2.)/pow(Ttmp,5.);
		}
		
        integ+=sgStar(T,paramrelic)*dark_entropy_Sigmad(T,paramrelic)/
                sqrt(1.+dark_density(T,paramrelic)/(pi*pi/30.*geff(T,paramrelic)*pow(T,4.)))/
                pow(heff(T,paramrelic),2.)/pow(T,5.)/2.;
		
		integ*=dlnT;

		double Mplanck=1.2209e19;	

		return 3.*Mplanck*sqrt(5./4./pi/pi/pi)*heff(T,paramrelic)*T*T*T*integ;	
	}
}

/*--------------------------------------------------------------*/

double dark_entropy_derivative(double T, struct relicparam paramrelic)
/* computes the dark energy entropy derivative at temperature T */
{
	if((paramrelic.sd0==0.)&&(paramrelic.Sigmad0==0.)) return 0.;
	
	if((paramrelic.Sigmad0==0.)&&(T<paramrelic.Tsend)) return 0.;
	
	if(paramrelic.Sigmad0==0.)
	{
		return (dark_entropy(T*1.001,paramrelic)-dark_entropy(T*0.999,paramrelic))/0.002/T;
	}
	else
	{
	
		double Mplanck=1.2209e19;	
        return 3.*sgStar(T,paramrelic)/T/heff(T,paramrelic)*(sqrt(geff(T,paramrelic))*
                                                             dark_entropy(T,paramrelic)
                                                             -sqrt(5.*Mplanck/4./pi/pi/pi)/T/T*
                                                             dark_entropy_Sigmad(T,paramrelic)/
                                                             sqrt(1.+dark_density(T,paramrelic)/
                                                                  (pi*pi/30.*geff(T,paramrelic)*pow(T,4.))));
	}
}

/*--------------------------------------------------------------*/

double dark_entropy_Sigmad(double T, struct relicparam paramrelic)
/* computes the dark entropy production at temperature T */
{
	if((paramrelic.sd0==0.)&&(paramrelic.Sigmad0==0.)) return 0.;
	
	if((paramrelic.Sigmad0==0.)&&(T<paramrelic.Tsend)) return 0.;
	
	double Mplanck=1.2209e19;
	
	if(paramrelic.Sigmad0==0.)
	{	
        return 1./Mplanck*(sqrt(24.*pi*(pi*pi/30.*geff(T,paramrelic)*pow(T,4.)
                                        +dark_density(T,paramrelic)))*dark_entropy(T,paramrelic)
                           -sqrt(4.*pi*pi*pi/45.)*heff(T,paramrelic)/
                           sgStar(T,paramrelic)*pow(T,3)*sqrt(1.+dark_density(T,paramrelic)/
                                                              (pi*pi/30.*geff(T,paramrelic)*pow(T,4.)))*
                           dark_entropy_derivative(T,paramrelic));
	}
	else
	{
		if(T<paramrelic.TSigmaend) return 0.;

		double s_rad_1MeV=2.*pi*pi/45.*heff(1.e-3,paramrelic)*1.e-9;
	
		double Sigma_rad_1MeV= 1./Mplanck*sqrt(4.*pi*pi*pi/5.*geff(1.e-3,paramrelic))*(1.e-6)*s_rad_1MeV;
	
		return paramrelic.Sigmad0*Sigma_rad_1MeV*pow(T/1.e-3,paramrelic.nSigmad);
	}
}

/*--------------------------------------------------------------*/

double nonthermal(double T, struct relicparam paramrelic)
/* computes the non-thermally produced relic particle number density at temperature T */
{
    if((paramrelic.nt0==0.)||(T<paramrelic.Tnend)) return 0.;
	
    return paramrelic.nt0*1.e-50*pow(T/1.e-3,paramrelic.nnt);
}

/*--------------------------------------------------------------*/

double neutdens(double Tnu, struct relicparam paramrelic)
/* Computes the neutrino density, including any effects from a neutrino degeneracy */
{
    if((paramrelic.xinu1==0.)&&(paramrelic.xinu2==0.)&&(paramrelic.xinu3==0.))
    {
        /* No degeneracy, relativistic approximation */
        return 2.*pi*pi/30.*7./8.*paramrelic.Nnu*pow(Tnu,4.);
    }

    int ie,je,n;
    double rho=0.;
    double xinu[4];
    double max1,max2,int1,int2;
    double x;

    xinu[1]=paramrelic.xinu1;
    xinu[2]=paramrelic.xinu2;
    xinu[3]=paramrelic.xinu3;

    /* SM neutrinos */
    for(ie=1;ie<=3;ie++)
    {
        /* The factor (paramrelic.Nnu/3.) includes extra DOF from non-rel. e+- and non-inst. nu decoupling */
        if(fabs(xinu[ie])<=0.03)
        {
            rho+=(paramrelic.Nnu/3.)*2.*pi*pi/30.*pow(Tnu,4.)*(7./8.+(15./(4*pi*pi))*xinu[ie]*xinu[ie]+(15./(8.*pow(pi,4.)))*
                                                               pow(xinu[ie],4.));
        }
        else if(fabs(xinu[ie])>=30.)
        {
            rho+=(paramrelic.Nnu/3.)*pow(Tnu,4.)/(8.*pi*pi)*pow(xinu[ie],4.)*(1.+12.*1.645/xinu[ie]/xinu[ie]);
        }
        else
        {
            /* Neutrinos */
            max1=(88.029+xinu[ie])*Tnu;
            int1=0.;
            n=50;
            for(je=1;je<=n-1;je++)
            {
                x=(double)je/(double)n*max1;
                int1+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu-xinu[ie]));
            }
            int1+=0.5*1./(2.*pi*pi)*pow(max1,3.)/(1.+exp(max1/Tnu-xinu[ie]));
            int1*=(paramrelic.Nnu/3.)*max1/(double)n;
            rho+=int1;

            /* Anti-neutrinos */
            max2=(88.029-xinu[ie])*Tnu;
            if(max2>0.)
            {
                int2=0.;
                n=50;
                for(je=1;je<=n-1;je++)
                {
                    x=(double)je/(double)n*max2;
                    int2+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu+xinu[ie]));
                }
                int2+=0.5/(2.*pi*pi)*pow(max2,3.)/(1.+exp(max2/Tnu+xinu[ie]));
                int2*=(paramrelic.Nnu/3.)*max2/(double)n;
                rho+=int2;
            }
        }
    }
    return rho;
}


/*--------------------------------------------------------------*/

double neutdens_deriv(double Tnu, struct relicparam paramrelic)
/* Computes the temperature (Tnu) derivative of the neutrino energy density */
{
    if((paramrelic.xinu1==0.)&&(paramrelic.xinu2==0.)&&(paramrelic.xinu3==0.))
    {
        return 7.*pi*pi/30.*paramrelic.Nnu*pow(Tnu,3.);
    }

    int ie,je,n;
    double drho=0.;
    double xinu[4];
    double max1,max2,int1,int2;
    double x;

    xinu[1]=paramrelic.xinu1;
    xinu[2]=paramrelic.xinu2;
    xinu[3]=paramrelic.xinu3;

    /* SM neutrinos */
    for(ie=1;ie<=3;ie++)
    {
        if(fabs(xinu[ie])<=0.03)
        {
            drho+=(paramrelic.Nnu/3.)*4.*pi*pi/15.*pow(Tnu,3.)*(7./8.+(15./(4*pi*pi))*xinu[ie]*xinu[ie]
                                                                +(15./(8.*pow(pi,4.)))*pow(xinu[ie],4.));
        }
        else if(fabs(xinu[ie])>=30.)
        {
            drho+=(paramrelic.Nnu/3.)*pow(Tnu,3.)/(2.*pi*pi)*pow(xinu[ie],4.)*(1.+12.*1.645/xinu[ie]/xinu[ie]);
        }
        else
        {
            max1=(88.029+xinu[ie])*Tnu;
            int1=0.;
            n=50;
            for(je=1;je<=n-1;je++)
            {
                x=(double)je/(double)n*max1;
                int1+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu-xinu[ie]));
            }
            int1+=0.5*1./(2.*pi*pi)*pow(max1,3.)/(1.+exp(max1/Tnu-xinu[ie]));
            int1*=(paramrelic.Nnu/3.)*4.*max1/Tnu/(double)n;
            drho+=int1;

            max2=(88.029-xinu[ie])*Tnu;
            if(max2>0.)
            {
                int2=0.;
                n=50;
                for(je=1;je<=n-1;je++)
                {
                    x=(double)je/(double)n*max2;
                    int2+=1./(2.*pi*pi)*pow(x,3.)/(1.+exp(x/Tnu+xinu[ie]));
                }
                int2+=0.5/(2.*pi*pi)*pow(max2,3.)/(1.+exp(max2/Tnu+xinu[ie]));
                int2*=(paramrelic.Nnu/3.)*4.*max2/Tnu/(double)n;
                drho+=int2;
            }
        }
    }
    return drho;
}




