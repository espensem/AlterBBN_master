'''
* ARES - A PROGRAM FOR COMPILING AND EXECUTING ALTERBBN.
*
* The program consists of the following classes:
* Ares() - Main class. Calls on the helper classes Element() and MakePlots(). 
*          Located in main.py.
* Element() - Helper class used to create objects of the elements for which AlterBBN
*             computes the resulting abundances. Holds four arrays: low, high, central and err,
*             which containts the AlterBBN results.
*             Located in element.py.
* MakePlots() - Helper class which holds methods for plotting. The user may easily create
*               new methods for customized plotting, which need to be called within the
*               main class Ares().
*               Located in plotting.py.
*
* This file processes the arguments needed to run the program, instantiates the main class
* Ares() and calls on its methods set_xvals() (for a multirun) and execute().
* The arguments may be set within this file or given as command line inputs upon the execution of 
* this file. Terminal inputs will override the arguments set within this file.
*
* USAGE:
*
* Inforun (Print a short program manual to terminal): 
* python ares.py info
*
* Singlerun (default -> cosmoType=paramfree, direct_output=False):
* python ares.py
*
* Singlerun (<optional argument> -> will override arguments set within this file):
* python ares.py single <cosmoType> <direct_output>
*
* Multirun (<optional argument> -> will override arguments set within this file):
* python ares.py multi <cosmoType> <varyParam> <lowVal> <highVal> <n_points> <spacingType> 
*
* The user may input as many optional arguments as desired, the command line arguments will
* then override the ones set here. HOWEVER, the order of the command line arguments must
* follow the order given above!
'''


'''
* General imports
'''
import sys
from main import Ares


'''
* User defined arguments and settings.
* Providing command line arguments upon the execution of this file will override
* the arguments and set here.
* Automatic plotting options must be set within the AlterBBN inputfile. 
''' 
#----- Allowed arguments and settings (DO NOT CHANGE THIS!) -----
runTypes = ["single", "multi"]
cosmoTypes = ["paramfree", "standard", "darkdens", "reheating", "wimp"]
varyParams = ["eta", "dnnu", "tau", "xinu1", "xinu2", "xinu3", "mass_wimp", "phiW"]
spacingTypes = ["log", "normal"]
#program = any (may be changed from the default="primary.c" if another main in AlterBBN is defined)
#direct_output = True or False
#lowVal = any float
#highVal = any float
#n_points = any int

###################### THIS MAY BE CHANGED BY THE USER ###########################
#----- General arguments and settings -----
program = "primary.c"	# The program to be run and compiled. Not a command line input option.
runType = 'single'	# Execute AlterBBN once (single) or multiple times (multi).
cosmoType = 'standard'	# Type of cosmology to be analyzed by AlterBBN.
direct_output = False	# Direct the results to 'alterdata.txt' (True),
			# or print them to the terminal (False) for singlerun.
			# Results are automatically directed to 'alterdata.txt' for multirun.
plot_none = False	# Switch to override all plot settings given in AlterBBN inputfile
			# and unable all plots if set to True.

#----- Free parameter arguments -----
varyParam = "eta"
lowVal = 1e-10
highVal = 1e-9
n_points = 3
spacingType = 'log'
##################################################################################


'''
* Method to validate command line input.
*
* @param arg: the input argument to be vaidated.
* @param arg_type: what type the argument should be, given as a string.
* Must be "int", "float", "bool" or "string".
'''
def check_input(arg, arg_type):
    try:
        if arg_type == "int":
            parameter = float(arg)	# In case n_points is inputted as a float
            parameter = int(parameter)
        elif arg_type == "float":
            parameter = float(arg)
        elif arg_type == "bool":
            # Allow for non-capitalized first letter. 
            if (arg == "true") or (arg == "True") or (arg == "TRUE"):
                parameter = True
            elif (arg == "false") or (arg == "False") or (arg == "FALSE"):
                parameter = False
            else:
                print "\t [ERROR] Wrong value for boolean argument.\n"+\
                    "\t         Try using only non-capitalized letters."
                sys.exit()
        elif arg_type == "string":
            # Type-safe the parameter input for the most common ways to write the keywords
            if (arg == "single") or (arg == "Single") or (arg == "SINGLE"):
                parameter = "single"
            elif (arg == "multi") or (arg == "Multi") or (arg == "MULTI"):
                parameter = "multi"
            elif (arg == "paramfree") or (arg == "Paramfree") or (arg == "PARAMFREE"):
                parameter = "paramfree"
            elif (arg == "standard") or (arg == "Standard") or (arg == "STANDARD"):
                parameter = "standard"
            elif (arg == "darkdens") or (arg == "Darkdens") or (arg == "DARKDENS") \
                 or (arg == "DarkDens") or (arg == "dark_dens") or (arg == "DARK_DENS") \
                 or (arg == "DARK_dens") or (arg == "dark_DENS") or (arg == "Dark_Dens") \
                 or (arg == "Dark_dens") or (arg == "dark_Dens"):
                parameter = "darkdens"
            elif (arg == "reheating") or (arg == "Reheating") or (arg == "REHEATING") \
                 or (arg == "reHeating") or (arg == "ReHeating") or (arg == "re_heating") \
                 or (arg == "re_Heating") or (arg == "Re_Heating") or (arg == "RE_HEATING"):
                parameter = "reheating"
            elif (arg == "wimp") or (arg == "Wimp") or (arg == "WIMP"):
                parameter = "wimp"
            elif (arg == "eta") or (arg == "Eta") or (arg == "ETA"):
                parameter = "eta"
            elif (arg == "tau") or (arg == "Tau") or (arg == "TAU"):
                parameter = "tau"
            elif (arg == "dnnu") or (arg == "dNnu") or (arg == "DNNU") or (arg == "Dnnu") \
                 or (arg == "DNnu") or (arg == "dnNU"):
                parameter = "dnnu"
            elif (arg == "xinu1") or (arg == "xinu_1") or (arg == "XINU1") or (arg == "XINU_1") \
                 or (arg == "Xinu1") or (arg == "Xinu_1"):
                parameter = "xinu1"
            elif (arg == "xinu2") or (arg == "xinu_2") or (arg == "XINU2") or (arg == "XINU_2") \
                 or (arg == "Xinu2") or (arg == "Xinu_2"):
                parameter = "xinu2"
            elif (arg == "xinu3") or (arg == "xinu_3") or (arg == "XINU3") or (arg == "XINU_3") \
                 or (arg == "Xinu3") or (arg == "Xinu_3"):
                parameter = "xinu3"
            elif (arg == "mass_wimp") or (arg == "mass_WIMP") or (arg == "MASS_WIMP") \
                 or (arg == "Mass_WIMP") or (arg == "Mass_Wimp") or (arg == "MASS_wimp") \
                 or (arg == "mass_Wimp") or (arg == "Mass_wimp") or (arg == "massWIMP") \
                 or (arg == "MASSwimp") or (arg == "MassWimp") or (arg == "MassWIMP") \
                 or (arg == "massWimp"):
                parameter = "mass_wimp"
            else:
                print "\t [ERROR] Wrong value for string argument.\n"+\
                    "\t         Try using only non-capitalized letters."
                sys.exit()
        else:
            print "\t [ERROR] Wrong argument in method 'check_input' in 'ares.py'.\n"+\
            "Must be ['int', 'float', 'bool', 'string']."
            sys.exit()
    except (ValueError, TypeError):
        print "\t [ERROR] Input argument could not be converted to the desired type.\n"

    return parameter
        
        

'''
* Process command line input
'''
inforun = False
if len(sys.argv) < 2:
    # Make default if no arguments are given
    runType = "single"
    cosmoType = "paramfree"
    direct_output = False
elif len(sys.argv) == 2:
    if sys.argv[1] == 'info':
        inforun = True
    else:
        runType = check_input(sys.argv[1], "string")
    if sys.argv[1] != "info":
        if sys.argv[1] not in runTypes:
            print "\t [ERROR] Wrong argument '%s'. Run with argument 'info' to read documentation." % \
                sys.argv[1]
            sys.exit()

elif len(sys.argv) == 3:
    runType = check_input(sys.argv[1], "string")
    cosmoType = check_input(sys.argv[2], "string")
    if runType not in runTypes:
        print "\t [ERROR] Wrong argument '%s'. Run with argument 'info' to read documentation." % sys.argv[1]
        sys.exit()
    elif cosmoType not in cosmoTypes:
        print "\t [ERROR] Wrong argument '%s'. Run with argument 'info' to read documentation." % sys.argv[2]
        sys.exit()

elif len(sys.argv) > 3:
    runType = check_input(sys.argv[1], "string")
    cosmoType = check_input(sys.argv[2], "string")

    # Maximum number of input arguments for single run is 3:
    if (len(sys.argv) > 4) and (sys.argv[1] == "single"):
        print "\t [ERROR] Too many input arguments for runtype 'single'."
        sys.exit()

    if len(sys.argv) == 4:
        if sys.argv[1] == "single":
            direct_output = check_input(sys.argv[3], "bool")
            print direct_output
            print type(direct_output)
        else:
            varyParam = check_input(sys.argv[3], "string")
    elif len(sys.argv) == 5:
        varyParam = check_input(sys.argv[3], "string")
        lowVal = check_input(sys.argv[4], "float")
    elif len(sys.argv) == 6:
        varyParam = check_input(sys.argv[3], "string")
        lowVal = check_input(sys.argv[4], "float")
        highVal = check_input(sys.argv[5], "float")
    elif len(sys.argv) == 7:
        varyParam = check_input(sys.argv[3], "string")
        lowVal = check_input(sys.argv[4], "float")
        highVal = check_input(sys.argv[5], "float")
        n_points = check_input(sys.argv[6], "int")
    elif len(sys.argv) == 8:
        varyParam = check_input(sys.argv[3], "string")
        lowVal = check_input(sys.argv[4], "float")
        highVal = check_input(sys.argv[5], "float")
        n_points = check_input(sys.argv[6], "int")
        spacingType = check_input(sys.argv[7], "string")
    else:
        print "\t [ERROR] Too many input arguments given."
        sys.exit()
            
    if runType not in runTypes:
        print "\t [ERROR] Wrong argument '%s'. Run without arguments to read documentation." % runType
        sys.exit()
    elif cosmoType not in cosmoTypes:
        print "\t [ERROR] Wrong argument '%s'. Run without arguments to read documentation." % cosmoType
        sys.exit()
    elif varyParam not in varyParams:
        print "\t [ERROR] Wrong argument '%s'. Run without arguments to read documentation." % varyParam
        sys.exit()
    elif spacingType not in spacingTypes:
        print "\t [ERROR] Wrong argument '%s'. Run without arguments to read documentation." % spacingType
        sys.exit()



'''
* Run the program, either by invoking Ares() or by printing short manual to terminal.
'''
if not inforun:
    print ''
    print "\t [INFO]  ARES is a tool to comile and run the C-written program AlterBBN.\n"+\
        "\t A user manual can be found by running the program with the argument 'info'."
    print ''
    if runType == "single":
        run = Ares(program=program, singlerun=True, cosmoType=cosmoType,
                   direct_output=direct_output, plot_none=plot_none)
        run.execute()
    elif runType == "multi":
        if cosmoType == "paramfree":
            print "\t [WARNING] Cosmology type 'PARAMFREE' is incompatible with run type \n"+\
                "\t 'MULTI'. Type of run changed to 'SINGLE'."
            run = Ares(program=program, singlerun=True, cosmoType=cosmoType,
                       direct_output=direct_output, plot_none=plot_none)
            run.execute()
        else:
            run = Ares(param=varyParam, program=program, singlerun=False, cosmoType=cosmoType,
                       direct_output=direct_output, plot_none=plot_none)
            run.set_xvals(lowVal, highVal, n_points, spacingType)
            run.execute()

else:
    print ""
    print "\t ARES - A PROGRAM FOR COMPILING AND EXECUTING ALTERBBN.\n\n"+\
        "\t The program consists of the following classes:\n"+\
	"\t Ares() - Main class. Calls on the helper classes Element() and MakePlots().\n"+\
    	"\t          Located in main.py.\n"+\
	"\t Element() - Helper class used to create objects of the elements for which\n"+\
	"\t             AlterBBN computes the resulting abundances. Holds four arrays:\n"+\
        "\t             low, high, central and err, which containts the AlterBBN results.\n"+\
	"\t             Located in element.py.\n"+\
	"\t MakePlots() - Helper class which holds methods for plotting. The user may easily\n"+\
        "\t               create new methods for customized plotting, which need to be called\n"+\
        "\t               within the main class Ares().\n"+\
	"\t               Located in plotting.py.\n\n"+\
	"\t ares.py processes the arguments needed to run the program, instantiates the main\n"+\
	"\t Ares() and calls on its methods set_xvals() (for a multirun) and execute().\n"+\
	"\t The arguments may be set within ares.py or given as command line inputs upon\n"+\
        "\t the execution of ares.py. Terminal inputs will override the arguments set within\n"+\
        "\t ares.py.\n\n"+\
        "\t USAGE:\n\n"+\
	"\t Inforun (Print a short program manual to terminal):\n"+\
	"\t python ares.py info\n\n"+\
	"\t Singlerun (default -> cosmoType=paramfree, direct_output=False):\n"+\
	"\t python ares.py\n\n"+\
	"\t Singlerun (<optional argument> -> will override arguments set within this file):\n"+\
	"\t python ares.py single <cosmoType> <direct_output>\n\n"+\
	"\t Multirun (<optional argument> -> will override arguments set within this file):\n"+\
	"\t python ares.py multi <cosmoType> <varyParam> <lowVal> <highVal> <n_points> "+\
        "<spacingType>\n\n"+\
	"\t The user may input as many optional arguments as desired, the command line\n"+\
        "\t arguments will then override the ones set in ares.py. HOWEVER, the order of the\n"+\
        "\t command line arguments must follow the ordering given above!\n\n"+\
        "\t Use the AlterBBN inputfile ('input.ini') for initialization of the cosmological\n"+\
        "\t parameters.\n\n"+\
        "\t LIST OF ALLOWED ARGUMENTS:\n\n"+\
        "\t 1) Type of run (runType).\n"+\
        "\t\t - 'SINGLE' for running AlterBBN once. Can be combined with arguments\n"+\
        "\t\t   <cosmoType> and <direct_output>.\n"+\
        "\t\t - 'MULTI' for running AlterBBN multiple times, with one of\n"+\
        "\t\t   [eta,dnnu,tau,xinu1,xinu2,xinu3,mass_wimp,phiW] set to vary.\n"+\
        "\t\t   Can be combined with arguments <cosmoType>, <varyparam>, <lowVal>, <highVal>,\n"+\
        "\t\t   <n_points> and <spacingType>.\n\n"+\
        "\t 2) Type of cosmology (<cosmoType>).\n"+\
        "\t\t - 'PARAMFREE' for running the SBBN scenario with fixed cosmological\n"+\
        "\t\t   parameters set within the method 'Init_cosmomodel' in AlterBBN's 'omega.c'.\n"+\
        "\t\t   Can not be combined with a multirun.\n"+\
        "\t\t - 'STANDARD' for running the SBBN scenario plus any optional neutrino\n"+\
        "\t\t   degeneracy, with cosmological parameters set in the AlterBBN inputfile.\n"+\
        "\t\t - 'DARKDENS' for including generalized physical effects that may change the\n"+\
        "\t\t   universal expansion rate, here parameterized through an additional\n"+\
        "\t\t   'dark energy' contribution in the first Friedmann equation, adding to the\n"+\
        "\t\t   total energy density of the universe.\n"+\
        "\t\t - 'REHEATING' for including generalized effects of increased entropy due to\n"+\
        "\t\t   reheating prior to BBN. This is parameterized through a 'dark energy\n"+\
        "\t\t   production'.\n"+\
        "\t\t - 'WIMP' for including scenarios with light WIMPs.\n\n"+\
        "\t If runType='SINGLE' (output always directed to outputfile for multirun):\n"+\
        "\t 3) Handling of result output (<direct_output>).\n"+\
        "\t\t - 'TRUE' for directing the AlterBBN result output to AlterBBN outputfile\n"+\
        "\t\t   ('alterdata.txt').\n"+\
        "\t\t - 'FALSE' for directing the AlterBBN result output to the terminal.\n\n"+\
        "\t If runType='MULTI':\n"+\
        "\t 3) Free parameter (<varyParam>).\n"+\
        "\t\t - 'ETA' for the baryon-to-photon ratio.\n"+\
        "\t\t - 'TAU' for the mean neutron lifetime.\n"+\
        "\t\t - 'DNNU' for the number of extra neutrino species.\n"+\
        "\t\t - 'XINU1' for the electron neutrino degeneracy parameter.\n"+\
        "\t\t - 'XINU2' for the muon neutrino degeneracy parameter.\n"+\
        "\t\t - 'XINU3' for the tau neutrino degeneracy parameter.\n"+\
        "\t\t - 'MASS_WIMP' for the WIMP mass.\n"+\
        "\t\t - 'PHIW' for the WIMP degeneracy parameter.\n\n"+\
        "\t 4) Lower limit on the free parameter (<lowVal>).\n"+\
        "\t\t - Any float.\n\n"+\
        "\t 5) Upper limit on the free parameter (<highVal>)\n"+\
        "\t\t - Any float.\n\n"+\
        "\t 6) Number of evaluation points for the free parameter (<n_points>).\n"+\
        "\t\t - Any integer.\n\n"+\
        "\t 7) Type of spacing of the x-values (<spacingType>).\n"+\
        "\t\t - 'NORMAL' for linearly spaced x-values.\n"+\
        "\t\t - 'LOG' for logarithmically spaced x-values."
    print ""
