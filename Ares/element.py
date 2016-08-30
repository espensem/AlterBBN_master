'''
* General imports.
'''
import numpy as np
import sys


'''
* Class to create objects of elements, holding information about each element in AlterBBN
* output in the form of arrays for low, high and central abundance values as well as the
* total error from the AlterBBN calculation. 
*
* @author Espen Sem Jenssen
'''
class Element():

    '''
    * Initializer.
    *
    * @param elemName: name of the element (e.g. 'He4', 'H2' etc...)
    * @param N: the number of points where the abundancies have been analyzed.
    * @param cnr: the column number in the AlterBBN outputfile where the results from
    * the element are found.
    '''
    def __init__(self, elemName, N, cnr):
        self.elemName = elemName
        self.N = N
        # Initialize arrays.
        self.low, self.cent, self.high, self.err = np.zeros(self.N), np.zeros(self.N), \
                                                   np.zeros(self.N), np.zeros(self.N)
        self.cnr = cnr



    '''
    * Method for appending values to the four arrays. Will append one of: low, high, cent or err.
    *
    * @param i: the evaluation point index (corresponds to the index in the array x in Ares()
    * where the corresponding varyParam value is found).
    * @param lowVal: set to the relevant value if the low value is to be appended.
    * @param centVal: set to the relevant value if the central value is to be appended.
    * @param highVal: set to the relevant value if the high value is to be appended.
    * @param errVal: set to the relevant value if the error is to be appended.
    *
    * @return void: appends a value to an array.
    '''
    def elemAppend(self, i, lowVal=None, centVal=None, highVal=None, errVal=None):
        '''
        Append either low, central, high or error value.
        '''
        if lowVal is not None:
            self.low[i] = lowVal
        elif centVal is not None:
            self.cent[i] = centVal
        elif highVal is not None:
            self.high[i] = highVal
        elif errVal is not None:
            self.err[i] = errVal
        else:
            print """\t[ERROR]  Need at least two arguments in function 'elemAppend'.\n\t\tElement.elemAppend(self, i, lowVal=None, centVal=None, highVal=None, errVal=None).""" 
            sys.exit()

