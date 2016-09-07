'''
* General imports.
'''
import os, sys
import numpy as np
import subprocess
import fileinput
import time
from element import Element
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from inspect import currentframe, getframeinfo
'''
* Import of user defined classes.
'''
from plotting import MakePlots


'''
* Main class for compiling and running AlterBBN.
*
* @author Espen Sem Jenssen
'''
class Ares():

    '''
    * Initializer with default arguments added.
    *
    * @param eta: Parameter to vary in multirun. Allowed arguments are
    * "eta", "nnu", "dnnu", "tau", "xinu1", "xinu2", "xinu3", "mass_wimp", "phiW".
    * @param program: AlterBBN main program to be compiled and run.
    * @param singlerun: True for a single compile and run, False for one compile
    * with multiple runs with param varied for each run.
    * @param cosmoType: Type of cosmology to execute. Allowed arguments are
    * "paramfree", "standard", "darkdens", "reheating", "wimp". "paramfree" is not
    * compatible with multirun.
    * @param direct_output: Only relevant in singlerun. True if output should be directed
    * to outputfile "alterdata.txt", False if output should be printed to the terminal.
    '''
    def __init__(self, param="eta", program="primary.c", cosmoType="paramfree", singlerun=True,
                 direct_output=False, plot_abundancesVSeta=False, plot_none=False):

        # Cross-platform path directories
        self.aresPath = os.path.dirname(__file__)      			# Current directory
        self.alterPath = os.path.split(self.aresPath)[0]		# AlterBBN directory
        self.plotPath = os.path.join(self.aresPath, "Plots")		# Plots/ folder
        self.program = program  					# C program
        self.executable = self.program.split(".")[0]+".x"		# C executable
        self.outfile = os.path.join(self.aresPath, "alterdata.txt")	# AlterBBN output-file
        self.inputfile = os.path.join(self.alterPath, "input.ini")	# AlterBBN input-file
        # General initializations
        self.param = param
        self.validParams1 = ["eta","xinu1","xinu2","xinu3","phiW"]
        self.validParams2 = ["dnnu","tau","mass_wimp"]
        self.cosmoType = cosmoType
        self.singlerun = singlerun
        self.direct_output = direct_output
        self.spacingType = ""
        self.x = 0
        # The following are declared as Element objects in method extract_data()
        self.He4, self.H2, self.He3, self.Li7, self.Li6, self.Be7 = 0, 0, 0, 0, 0, 0
        # Dictionary holding all element objects
        self.elements = {}
        # Plotting switches, defaults to False for all
        self.plot_abundancesVSeta = plot_abundancesVSeta
        self.plot_abundancesVStime = False
        self.plot_NeffVSmassWIMP = False
        # Switch to override all plot settings given in AlterBBN inputfile, and unable all.
        self.plot_none = plot_none

        '''
        print self.aresPath
        print self.alterPath
        print self.program
        print self.executable
        sys.exit()
        '''



    '''
    * Method for setting the x values in multirun.
    *
    * @param low: Minimum value of the parameter to vary
    * @param high: Maximum value of the parameter to vary
    * @param N: Number of values in the interval [low, high] to analyze
    * @param spacing: The type of spacing between the values in the interval [low, high].
    * May be 'normal' for normally spaced values and 'log' for logarithmically spaced values.
    *
    * @return void: sets self.x
    '''
    def set_xvals(self, low, high, N, spacing):
        self.spacingType = spacing
        if spacing == 'normal':
            self.x = np.linspace(low, high, N)
        elif spacing == 'log':
            self.x = np.logspace(np.log10(low), np.log10(high), num=N)



    '''
    * Method for compiling and executing AlterBBN, given the default/user defined settings
    * and arguments.
    *
    * @return void: compiles self.program and executes self.executable either
    * once for singlerun or N times for multirun.
    '''
    def execute(self):
        # Print run info to terminal
        if self.singlerun:
            print "\t [INFO]   Running Ares with the settings:\n"+\
                "\t          runType=single, cosmoType=%s, direct_output=%s\n" \
                % (self.cosmoType, self.direct_output)
        else:
            print "\t [INFO]   Running Ares with the settings:\n"+\
                "\t          runType=multi, cosmoType=%s, varyParam=%s," \
                % (self.cosmoType, self.param)
            print "\t          lowVal=%.3e, highVal=%.3e, nVals=%d, spacingType=%s.\n" \
                % (self.x[0], self.x[-1], len(self.x), self.spacingType)

        t0 = time.time()
        # Change working directory to AlterBBN directory
        os.chdir(self.alterPath)
        # Compile .c-file
        try:
            frameinfo = getframeinfo(currentframe())
            subprocess.check_call(["make", self.program])
            print ""
            print "\t [INFO]   '%s'" % (os.path.join(os.getcwd(), self.program))
            print "\t          compiled successfully." 
        except subprocess.CalledProcessError as e:
            print frameinfo.filename, frameinfo.lineno+1
            print "\t [ERROR]  Could not compile '%s'." % (os.path.join(os.getcwd(),
                                                                        self.program))
            print "\t          Check if file exists. Check for errors in the makefile."
            sys.exit()

        # SINGLERUN
        if self.singlerun:
            print "\t [RESULT]"
            # Output directed to self.outfile
            if self.direct_output:
                # Clear outfile, and prepare for AlterBBN output
                open(self.outfile, 'w').close()
                outapp = open(self.outfile, 'a')
                with outapp as out:
                    # The SBBN scenario, using default parameters, requires no input arguments
                    # for AlterBBN executable
                    if self.cosmoType == "paramfree":
                        try:
                            frameinfo = getframeinfo(currentframe())
                            subprocess.call(["./"+self.executable], stdout=out)
                        except Exception as e:
                            print frameinfo.filename, frameinfo.lineno+1
                            print "\t [ERROR]  Could not run '%s'." % (os.path.
                                                                       join(os.getcwd(),
                                                                            self.executable))
                            print "\t          Check if executable exists."
                            print "\t          Check for errors in the makefile."
                            print "\t          Check if outfile %s exists." % self.outfile
                            sys.exit()
                    # Other scenarios must be given as input argument for AlterBBN executable
                    else:
                        try:
                            frameinfo = getframeinfo(currentframe())
                            subprocess.call(["./"+self.executable, self.cosmoType], stdout=out)
                        except Exception as e:
                            print frameinfo.filename, frameinfo.lineno+1
                            print "\t [ERROR]  Could not run '%s'." % (os.path.
                                                                       join(os.getcwd(),
                                                                            self.executable))
                            print "\t          Check if executable exists."
                            print "\t          Check for errors in the makefile."
                            print "\t          Check if outfile %s exists." % self.outfile
                            sys.exit()
                    print ''
                outapp.close()
            # Output to terminal
            else:
                if self.cosmoType == "paramfree":
                    try:
                        frameinfo = getframeinfo(currentframe())
                        subprocess.call(["./"+self.executable])
                    except Exception as e:
                        print frameinfo.filename, frameinfo.lineno+1
                        print "\t [ERROR]  Could not run '%s'." % (os.path.
                                                                   join(os.getcwd(),
                                                                        self.executable))
                        print "\t          Check if executable exists."
                        print "\t          Check for errors in the makefile."
                        sys.exit()
                else:
                    try:
                        frameinfo = getframeinfo(currentframe())
                        subprocess.call(["./"+self.executable, self.cosmoType])
                    except Exception as e:
                        print frameinfo.filename, frameinfo.lineno+1
                        print "\t [ERROR]  Could not run '%s'." % (os.path.
                                                                   join(os.getcwd(),
                                                                        self.executable))
                        print "\t          Check if executable exists."
                        print "\t          Check for errors in the makefile."
                        sys.exit()
                print ''
            os.chdir(self.aresPath)

        # MULTIRUN
        else:
            # Clear outfile, and prepare for AlterBBN output
            open(self.outfile, 'w').close()
            # Iterate through the self.param values that are to be analyzed 
            for i in range(len(self.x)):
                # Adjust terminal info to fit self.param
                if self.param in self.validParams1:
                    print "\t [INFO]   Iteration %d: %s=%.3e" % (i+1, self.param, self.x[i])
                else:
                    print "\t [INFO]   Iteration %d: %s=%.3f" % (i+1, self.param, self.x[i])
                # Search through self.inputfile to find the parameter to vary, and change it
                for line in fileinput.input(self.inputfile, inplace=True):
                    if '=' in line and line.startswith(self.param):
                        # Overwrite the value of self.param in self.inputfile, adjusting the
                        # representation form depending on what self.param is
                        # BE CAREFUL NOT TO ERASE IMPORTANT INFORMATION FROM SELF.INPUTFILE
                        # IF MODIFYING THIS!
                        if self.param in self.validParams1:
                            print self.param+' = %.3e' % self.x[i]
                        else:
                            print self.param+' = %.3f' % self.x[i]
                        continue
                    else:
                        # Do nothing if self.param is not in line
                        print line.strip()
                # Execute the program, direct output to self.outfile
                outapp = open(self.outfile, 'a')
                with outapp as out:
                    subprocess.call(["./"+self.executable, self.cosmoType], stdout=out)
            outapp.close()
            # Change working directory back to Ares directory
            os.chdir(self.aresPath)

        self.t1 = time.time() - t0
        print "\t [INFO]   Time elapsed reading, compiling and running:\n\t\t  %.4f seconds" % self.t1
        # EXECUTE PLOTTING METHODS
        if not self.singlerun:
            # Check self.inputfile if the user wants to make any plots
            for line in fileinput.input(self.inputfile):
                if "=" in line and line.startswith("abundancesVSeta"):
                    if line.split("=")[1].strip() == "true" or \
                       line.split("=")[1].strip() == "True" or \
                       line.split("=")[1].strip() == "TRUE":
                        self.plot_abundancesVSeta = True
                    elif line.split("=")[1].strip() == "false" or \
                         line.split("=")[1].strip() == "False" or \
                         line.split("=")[1].strip() == "FALSE":
                        self.plot_abundancesVSeta = False
                    else:
                        print "\t [WARNING] Invalid value '%s' for parameter '%s' " \
                            % (line.split("=")[1].strip(), line.split("=")[0].strip()) +\
                            "in AlterBBN inputfile:" 
                        print "\t           %s" % self.inputfile
                        print "\t           Default value 'False' will be used:"
                        self.plot_abundancesVSeta = False
                    
            print "\t [INFO]   Initializing plotting"
            t0 = time.time()
            # Extract data from file
            self.extract_data()
            # Initialize different plots
            makePlots = MakePlots()
            if self.plot_abundancesVSeta and not self.plot_none:
                makePlots.abundancesVSeta("eta", self.x, self.elements)

            # ADD USER DEFINED PLOTTING METHODS HERE 

            print "\t [INFO]   Saved plot 'abundancesVSeta.eps' to folder Plots/"

            self.t2 = time.time() - t0
            print "\t [INFO]   Plotting finished"
            print "\t [INFO]   Time elapsed plotting:\n\t\t  %.4f seconds" % self.t2
            print "\t [INFO]   Total elapsed time:\n\t\t  %.4f seconds" % (self.t1+self.t2)




    '''
    * Method for extracting AlterBBN output.
    *
    * @return void: Collects data from self.outfile and instantiates an Element object for
    * each element found in self.outfile. Each object is added to the self.elements dictionary.
    '''
    def extract_data(self):
        # Check if self.outfile exists
        if not os.path.exists(self.outfile):
            print "\t [ERROR]  The file '%s' does not exist. Exiting program." %self.outfile
            print ""
            sys.exit()
        # Find the length of the arrays
        infile = open(self.outfile, 'r')
        infile.readline()
        N = 0
        for line in infile:
            if 'Yp' in line: N += 1
        infile.close()
        # Find the ordering of the elements in self.outfile
        colnrs = {'He4':0, 'H2':0, 'He3':0, 'Li7':0, 'Li6':0, 'Be7':0}
        infile = open(self.outfile, 'r')
        for line in infile:
            if not line.strip(): continue	# If there are blank lines, skip them
            words = line.split()
            for i in range(len(words)):
                if words[i].startswith('Yp'): colnrs['He4'] = i
                elif words[i].startswith('H2'): colnrs['H2'] = i
                elif words[i].startswith('He3'): colnrs['He3'] = i
                elif words[i].startswith('Li7'): colnrs['Li7'] = i
                elif words[i].startswith('Li6'): colnrs['Li6'] = i
                elif words[i].startswith('Be7'): colnrs['Be7'] = i
            break		# Need only consider first line
        infile.close()
        # Create Element objects 
        self.He4 = Element('He4', N, colnrs['He4'])
        self.H2 = Element('H2', N, colnrs['H2'])
        self.He3 = Element('He3', N, colnrs['He3'])
        self.Li7 = Element('Li7', N, colnrs['Li7'])
        self.Li6 = Element('Li6', N, colnrs['Li6'])
        self.Be7 = Element('Be7', N, colnrs['Be7'])
        # Fill Element object arrays with data from self.outfile
        infile = open(self.outfile, 'r')
        i = 0
        for line in infile:
            if not line.strip(): continue	# Skip blank lines
            if '[RESULT]' in line:
                 i += 1
            elif 'cent' in line:
                colonsplit = line.split(':')
                values = colonsplit[1].split()
                self.He4.elemAppend(i, centVal=float(values[colnrs['He4']]))
                self.H2.elemAppend(i, centVal=float(values[colnrs['H2']]))
                self.He3.elemAppend(i, centVal=float(values[colnrs['He3']]))
                self.Li7.elemAppend(i, centVal=float(values[colnrs['Li7']]))
                self.Li6.elemAppend(i, centVal=float(values[colnrs['Li6']]))
                self.Be7.elemAppend(i, centVal=float(values[colnrs['Be7']]))
            elif '+/-' in line:
                colonsplit = line.split(':')
                values = colonsplit[1].split()
                self.He4.elemAppend(i, errVal=float(values[colnrs['He4']]))
                self.H2.elemAppend(i, errVal=float(values[colnrs['H2']]))
                self.He3.elemAppend(i, errVal=float(values[colnrs['He3']]))
                self.Li7.elemAppend(i, errVal=float(values[colnrs['Li7']]))
                self.Li6.elemAppend(i, errVal=float(values[colnrs['Li6']]))
                self.Be7.elemAppend(i, errVal=float(values[colnrs['Be7']]))
            elif 'low' in line:
                colonsplit = line.split(':')
                values = colonsplit[1].split()
                self.He4.elemAppend(i, lowVal=float(values[colnrs['He4']]))
                self.H2.elemAppend(i, lowVal=float(values[colnrs['H2']]))
                self.He3.elemAppend(i, lowVal=float(values[colnrs['He3']]))
                self.Li7.elemAppend(i, lowVal=float(values[colnrs['Li7']]))
                self.Li6.elemAppend(i, lowVal=float(values[colnrs['Li6']]))
                self.Be7.elemAppend(i, lowVal=float(values[colnrs['Be7']]))
            elif 'high' in line:
                colonsplit = line.split(':')
                values = colonsplit[1].split()
                self.He4.elemAppend(i, highVal=float(values[colnrs['He4']]))
                self.H2.elemAppend(i, highVal=float(values[colnrs['H2']]))
                self.He3.elemAppend(i, highVal=float(values[colnrs['He3']]))
                self.Li7.elemAppend(i, highVal=float(values[colnrs['Li7']]))
                self.Li6.elemAppend(i, highVal=float(values[colnrs['Li6']]))
                self.Be7.elemAppend(i, highVal=float(values[colnrs['Be7']]))
            
        infile.close()
        # Add Element objects to self.elements dictionary
        self.elements["He4"] = self.He4
        self.elements["H2"] = self.H2
        self.elements["He3"] = self.He3
        self.elements["Li7"] = self.Li7
        self.elements["Li6"] = self.Li6
        self.elements["Be7"] = self.Be7
        print "\t [INFO]   Data successfully extracted from:\n"+\
            "\t          '%s'." % self.outfile
        if self.param in self.validParams1:
            print """\t [INFO]   Varying the parameter '%s' with %d %s spaced values between %.3e and\n\t\t  %.2e""" % (self.param, len(self.x), self.spacingType, min(self.x), max(self.x))
        else:
            print """\t [INFO]   Varying the parameter '%s' with %d %s spaced values between %.3f and\n\t\t  %.2f""" % (self.param, len(self.x), self.spacingType, min(self.x), max(self.x))



            
            


