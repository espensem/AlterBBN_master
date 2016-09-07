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
* Imports
'''
import sys, os
from main import Ares
from check_settings import GetParamDict


'''
* Define allowed arguments for run, to validate user input.
'''
runTypes = ["single", "multi"]
cosmoTypes = ["paramfree", "standard", "darkdens", "reheating", "wimp"]
varyParams = ["eta", "dnnu", "tau", "xinu1", "xinu2", "xinu3", "mass_wimp", "phiW"]
spacingTypes = ["log", "normal"]



'''
* Collect user input from settings.txt
'''
def CollectData(check_params):
    # Initiate parameter dictionary
    parameters = {}

    filePath = os.path.join(os.path.dirname(__file__), "settings.txt")
    # Check if 'settings.txt' exists. If not, create a new file with default parameters,
    # and run with those.
    if os.path.isfile(filePath):
        # File exists, proceeding
        pass
    else:
        pass
        # create new file
        # filePath = new file
    
    # Read the user defined values and put into dictionary
    infile = open(filePath, "r")
    for line in infile:
        if not line.strip() or line.startswith("#"):
            # Empty line, or commentline
            pass
        elif line.split()[1] == "program":
            parameters["program"] = line.split("=")[1].strip()
            check_params["program"] = True
        elif line.split()[1] == "runType":
            parameters["runType"] = CheckInput(line.split("=")[1].strip(), "string")
            check_params["runType"] = True
        elif line.split()[1] == "cosmoType":
            parameters["cosmoType"] = CheckInput(line.split("=")[1].strip(), "string")
            check_params["cosmoType"] = True
        elif line.split()[1] == "direct_output":
            parameters["direct_output"] = CheckInput(line.split("=")[1].strip(), "bool")
            check_params["direct_output"] = True
        elif line.split()[1] == "plot_none":
            parameters["plot_none"] = CheckInput(line.split("=")[1].strip(), "bool")
            check_params["plot_none"] = True
        elif line.split()[1] == "plot_abundancesVSeta":
            parameters["plot_abundancesVSeta"] = CheckInput(line.split("=")[1].strip(), "bool")
            check_params["plot_abundancesVSeta"] = True
        elif line.split()[1] == "varyParam":
            parameters["varyParam"] = CheckInput(line.split("=")[1].strip(), "string")
            check_params["varyParam"] = True
        elif line.split()[1] == "lowVal":
            parameters["lowVal"] = CheckInput(line.split("=")[1].strip(), "float")
            check_params["lowVal"] = True
        elif line.split()[1] == "highVal":
            parameters["highVal"] = CheckInput(line.split("=")[1].strip(), "float")
            check_params["highVal"] = True
        elif line.split()[1] == "nVals":
            parameters["nVals"] = CheckInput(line.split("=")[1].strip(), "int")
            check_params["nVals"] = True
        elif line.split()[1] == "spacingType":
            parameters["spacingType"] = CheckInput(line.split("=")[1].strip(), "string")
            check_params["spacingType"] = True
        else:
            pass

    infile.close()
    # If not all the parameters are defined, exit with message
    for key in check_params:
        if not check_params[key]:
            print "\t [ERROR] In reading 'settings.txt'. Check that parameter '%s' exists." \
                % key
            sys.exit()
    
    return parameters
        


'''
* Process command line input. Arguments given here will override the corresponding 
* arguments given in settings.txt.
'''
def ProcessInput(parameters):
    inforun = False
    if len(sys.argv) < 2:
        # Make default if no arguments are given
        parameters["runType"] = "single"
        parameters["cosmoType"] = "paramfree"
        parameters["direct_output"] = False
    elif len(sys.argv) == 2:
        if sys.argv[1] == 'info':
            inforun = True
        else:
            if sys.argv[1] not in runTypes:
                print "\t [ERROR] Wrong argument '%s'. Run with argument " % sys.argv[1] +\
                    "'info' to read documentation." 
                sys.exit()
            else:
                parameters["runType"] = CheckInput(sys.argv[1], "string")
        
    elif len(sys.argv) == 3:
        parameters["runType"] = CheckInput(sys.argv[1], "string")
        parameters["cosmoType"] = CheckInput(sys.argv[2], "string")
        if sys.argv[1] not in runTypes:
            print "\t [ERROR] Wrong argument '%s'. Run with argument " % sys.argv[1] +\
                    "'info' to read documentation." 
            sys.exit()
        elif sys.argv[2] not in cosmoTypes:
            print "\t [ERROR] Wrong argument '%s'. Run with argument " % sys.argv[2] +\
                    "'info' to read documentation." 
            sys.exit()

    elif len(sys.argv) > 3:
        parameters["runType"] = CheckInput(sys.argv[1], "string")
        parameters["cosmoType"] = CheckInput(sys.argv[2], "string")

        # Maximum number of input arguments for single run is 3:
        if (len(sys.argv) > 4) and (sys.argv[1] == "single"):
            print "\t [ERROR] Too many input arguments for runtype 'single'."
            sys.exit()

        if len(sys.argv) == 4:
            if sys.argv[1] == "single":
                parameters["direct_output"] = CheckInput(sys.argv[3], "bool")
            else:
                parameters["varyParam"] = CheckInput(sys.argv[3], "string")
        elif len(sys.argv) == 5:
            parameters["varyParam"] = CheckInput(sys.argv[3], "string")
            parameters["lowVal"] = CheckInput(sys.argv[4], "float")
        elif len(sys.argv) == 6:
            parameters["varyParam"] = CheckInput(sys.argv[3], "string")
            parameters["lowVal"] = CheckInput(sys.argv[4], "float")
            parameters["highVal"] = CheckInput(sys.argv[5], "float")
        elif len(sys.argv) == 7:
            parameters["varyParam"] = CheckInput(sys.argv[3], "string")
            parameters["lowVal"] = CheckInput(sys.argv[4], "float")
            parameters["highVal"] = CheckInput(sys.argv[5], "float")
            parameters["nVals"] = CheckInput(sys.argv[6], "int")
        elif len(sys.argv) == 8:
            parameters["varyParam"] = CheckInput(sys.argv[3], "string")
            parameters["lowVal"] = CheckInput(sys.argv[4], "float")
            parameters["highVal"] = CheckInput(sys.argv[5], "float")
            parameters["nVals"] = CheckInput(sys.argv[6], "int")
            parameters["spacingType"] = CheckInput(sys.argv[7], "string")
        else:
            print "\t [ERROR] Too many input arguments given."
            sys.exit()

        if runType not in runTypes:
            print "\t [ERROR] Wrong argument '%s'. Run without arguments to read documentation." \
                % parameters["runType"]
            sys.exit()
        elif cosmoType not in cosmoTypes:
            print "\t [ERROR] Wrong argument '%s'. Run without arguments to read documentation." \
                % parameters["cosmoType"]
            sys.exit()
        elif varyParam not in varyParams:
            print "\t [ERROR] Wrong argument '%s'. Run without arguments to read documentation." \
                % parameters["varyParam"]
            sys.exit()
        elif spacingType not in spacingTypes:
            print "\t [ERROR] Wrong argument '%s'. Run without arguments to read documentation." \
                % parameters["spacingType"]
            sys.exit()

    return inforun, parameters


'''
* Method to validate command line input.
*
* @param arg: the input argument to be vaidated.
* @param arg_type: what type the argument should be, given as a string.
* Must be "int", "float", "bool" or "string".
'''
def CheckInput(arg, arg_type):
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
            elif (arg == "normal") or (arg == "Normal") or (arg == "NORMAL"):
                parameter = "normal"
            elif (arg == "log") or (arg == "Log") or (arg == "LOG"):
                parameter = "log"
            else:
                print "\t [ERROR] Wrong value for string argument.\n"+\
                    "\t         Try using only non-capitalized letters."
                sys.exit()
        else:
            print "\t [ERROR] Wrong argument in method 'CheckInput' in 'ares.py'.\n"+\
            "Must be ['int', 'float', 'bool', 'string']."
            sys.exit()
    except (ValueError, TypeError):
        print "\t [ERROR] Input argument could not be converted to the desired type.\n"

    return parameter



if __name__ == "__main__":

    # Get dictionary of valid input parameter names, to check for valid input in settings.txt
    check_params = GetParamDict()
    # Collect data from settings.txt, returned as a dictionary
    parameters = CollectData(check_params)
    # Update the dictionary if input arguments are given. Also get inforun switch.
    inforun, parameters = ProcessInput(parameters)
    '''
    * Run the program, either by invoking Ares() or by printing short manual to terminal.
    '''
    if not inforun:
        print ''
        print "\t [INFO]  ARES is a tool to comile and run the C-written program AlterBBN.\n"+\
            "\t A user manual can be found by running the program with the argument 'info'."
        print ''
        if parameters["runType"] == "single":
            run = Ares(program=parameters["program"], singlerun=True,
                       cosmoType=parameters["cosmoType"],
                       direct_output=parameters["direct_output"])
            run.execute()
        elif parameters["runType"] == "multi":
            if parameters["cosmoType"] == "paramfree":
                print "\t [WARNING] Cosmology type 'PARAMFREE' is incompatible with run type \n"+\
                    "\t 'MULTI'. Type of run changed to 'SINGLE'."
                run = Ares(program=parameters["program"], singlerun=True,
                       cosmoType=parameters["cosmoType"],
                       direct_output=parameters["direct_output"])
                run.execute()
            else:
                run = Ares(param=parameters["varyParam"],
                           program=parameters["program"], singlerun=False,
                           cosmoType=parameters["cosmoType"],
                           direct_output=parameters["direct_output"],
                           plot_abundancesVSeta=parameters["plot_abundancesVSeta"],
                           plot_none=parameters["plot_none"])
                run.set_xvals(parameters["lowVal"], parameters["highVal"],
                              parameters["nVals"], parameters["spacingType"])
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
            "\t 6) Number of evaluation points for the free parameter (<nVals>).\n"+\
            "\t\t - Any integer.\n\n"+\
            "\t 7) Type of spacing of the x-values (<spacingType>).\n"+\
            "\t\t - 'NORMAL' for linearly spaced x-values.\n"+\
            "\t\t - 'LOG' for logarithmically spaced x-values."
        print ""
