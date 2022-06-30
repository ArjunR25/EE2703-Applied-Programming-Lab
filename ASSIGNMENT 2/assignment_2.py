'''
SPICE - PART 2
-Reads input netlist file to solve the circuit, works for dc and single frequency ac circuits

-Uses Modified Nodal Analysis(MNA) to report nodal voltages and currents through voltage sources
-Reports the same quantities in phasors in case of ac circuit

Feb. 09, 2022

By Arjun R
EE20B016
'''

import numpy as np
from sys import argv, exit

"""
It's recommended to use constant variables than hard-coding them everywhere.
"""
CIRCUIT = '.circuit'
END = '.end'
AC = '.ac'

"""
Check if the user has given required and only the required inputs
Otherwise, show them the expected usage.
"""
if len(argv) != 2:
    print('\nUsage: %s <inputfile>' % argv[0])
    exit()

class RLC:
    """
    Creates an object for passive elements R, L, C
    Methods: applyToMatrix  
    """
    def __init__(self, elementName, inNode, outNode, value):
        self.eleName = elementName
        self.inNode = inNode
        self.outNode = outNode
        self.value = value

    def applyToMatrix(self, nodeMatrix):
        """method to apply the corresponding stamp for the given element 
        to the incidence matrix"""
        if self.eleName[0] == 'R':
            G = 1/self.value
            nodeMatrix[Table[self.inNode], Table[self.inNode]] += G
            nodeMatrix[Table[self.inNode], Table[self.outNode]] -= G
            nodeMatrix[Table[self.outNode], Table[self.inNode]] -= G
            nodeMatrix[Table[self.outNode], Table[self.outNode]] += G

        elif self.eleName[0] == 'L':
            YL = 1/complex(0, w*self.value)
            nodeMatrix[Table[self.inNode], Table[self.inNode]] += YL
            nodeMatrix[Table[self.inNode], Table[self.outNode]] -= YL
            nodeMatrix[Table[self.outNode], Table[self.inNode]] -= YL
            nodeMatrix[Table[self.outNode], Table[self.outNode]] += YL

        elif self.eleName[0] == 'C':
            YC = complex(0, w*self.value)
            nodeMatrix[Table[self.inNode], Table[self.inNode]] += YC
            nodeMatrix[Table[self.inNode], Table[self.outNode]] -= YC
            nodeMatrix[Table[self.outNode], Table[self.inNode]] -= YC
            nodeMatrix[Table[self.outNode], Table[self.outNode]] += YC

class IndSrcs:
    """
    Creates an object for independent sources
    Methods: applyToMatrix  
    """
    def __init__(self, elementName, inNode, outNode, value, phase = 0):
        self.eleName = elementName
        self.inNode = inNode
        self.outNode = outNode
        self.value = value
        self.phase = phase

    def applyToMatrix(self, nodeMatrix, bMatrix):
        """method to apply the corresponding stamp for the given element 
        to the incidence and source matrices"""
        if self.eleName[0] == 'V':
            nodeMatrix[Table[self.inNode], Table[self.eleName]] += 1
            nodeMatrix[Table[self.outNode], Table[self.eleName]] -= 1
            nodeMatrix[Table[self.eleName], Table[self.inNode]] += 1
            nodeMatrix[Table[self.eleName], Table[self.outNode]] -= 1
            bMatrix[Table[self.eleName]] += self.value*complex(np.cos(self.phase), np.sin(self.phase)) # phasor form

        elif self.eleName[0] == 'I':
            bMatrix[Table[self.inNode]] -= self.value*complex(np.cos(self.phase), np.sin(self.phase))
            bMatrix[Table[self.outNode]] += self.value*complex(np.cos(self.phase), np.sin(self.phase))

class ICtrlSrcs:
    """
    Creates an object for current-controlled sources
    Methods: applyToMatrix  
    """
    def __init__(self, elementName, inNode, outNode, VnodeIn, VnodeOut, Vname, value):
        self.eleName = elementName
        self.inNode = inNode
        self.outNode = outNode
        self.VIn = VnodeIn # in node for 0V source through which current enters
        self.VOut = VnodeOut # out node for 0V source through which current leaves
        self.Vname = Vname
        self.value = value
    
    def applyToMatrix(self, nodeMatrix):
        """method to apply the corresponding stamp for the given element 
        to the incidence matrix"""
        if self.eleName[0] == 'H': #CCVS
            nodeMatrix[Table[self.inNode], Table[self.eleName]] += 1
            nodeMatrix[Table[self.outNode], Table[self.eleName]] -= 1
            #nodeMatrix[Table[self.VIn], Table[self.Vname]] += 1
            #nodeMatrix[Table[self.VOut], Table[self.Vname]] -= 1
            #nodeMatrix[Table[self.Vname], Table[self.VIn]] += 1
            #nodeMatrix[Table[self.Vname], Table[self.VOut]] -= 1
            nodeMatrix[Table[self.eleName], Table[self.inNode]] += 1
            nodeMatrix[Table[self.eleName], Table[self.outNode]] -= 1
            nodeMatrix[Table[self.eleName], Table[self.Vname]] -= self.value
            
        if self.eleName[0] == 'F': #CCCS
            nodeMatrix[Table[self.inNode], Table[self.Vname]] += self.value
            nodeMatrix[Table[self.outNode], Table[self.Vname]] -= self.value
            #nodeMatrix[Table[self.VIn], Table[self.Vname]] += 1
            #nodeMatrix[Table[self.VOut], Table[self.Vname]] -= 1
            #nodeMatrix[Table[self.Vname], Table[self.VIn]] += 1
            #nodeMatrix[Table[self.Vname], Table[self.VOut]] -= 1

class VCtrlSrcs:
    """
    Creates an object for voltage-controlled sources
    Methods: applyToMatrix  
    """
    def __init__(self, elementName, inNode, outNode, VnodeIn, VnodeOut, value):
        self.eleName = elementName
        self.inNode = inNode
        self.outNode = outNode
        self.VIn = VnodeIn
        self.VOut = VnodeOut
        self.value = value

    def applyToMatrix(self, nodeMatrix):
        """method to apply the corresponding stamp for the given element 
        to the incidence matrix"""
        if self.eleName[0] == 'E': #VCVS
            nodeMatrix[Table[self.inNode], Table[self.eleName]] += 1
            nodeMatrix[Table[self.outNode], Table[self.eleName]] -= 1
            nodeMatrix[Table[self.eleName], Table[self.inNode]] += 1
            nodeMatrix[Table[self.eleName], Table[self.outNode]] -= 1
            nodeMatrix[Table[self.eleName], Table[self.VIn]] -= self.value
            nodeMatrix[Table[self.eleName], Table[self.VOut]] += self.value

        if self.eleName[0] == 'G': #VCCS
            nodeMatrix[Table[self.inNode], Table[self.VIn]] += self.value
            nodeMatrix[Table[self.inNode], Table[self.VOut]] -= self.value
            nodeMatrix[Table[self.outNode], Table[self.VIn]] -= self.value
            nodeMatrix[Table[self.outNode], Table[self.VOut]] += self.value

def get_key(val, dictionary):
    """Gets corrsponding key for provided value from the given dictionary"""
    for key, value in dictionary.items():
         if val == value:
             return key


"""
The user might input a wrong file name by mistake.
In this case, the open function will throw an IOError.
Taking care of it using try-catch:
"""
try:
    with open(argv[1]) as f:
        lines = f.readlines()
        start = -1; end = -2; ac = -1

        for line in lines:
            # extracting circuit definition start and end lines
            if AC == line[:len(AC)]:
                ac = lines.index(line)
            if CIRCUIT == line[:len(CIRCUIT)]:
                start = lines.index(line)
            elif END == line[:len(END)]:
                end = lines.index(line)

        if start >= end:                # validating circuit block
            print('Invalid circuit definition')
            exit(0)

        # ac here denotes the line number where .ac command appears in the netlist
        # in case there is no .ac command, the value of ac remains -1 as initialized
        if ac == -1:
            w = 0 # w is the angular frequency; for dc, its value is zero
        else:
            line_AC = lines[ac].split('#')[0].split() # list of tokens in .ac command
            acSrc = line_AC[1] # name of AC source
            w = 2*np.pi*float(line_AC[2])

        Table = {'GND': 0,} # Table of distinct nodes and voltage sources(to define rows (and columns) of the incidence matrix)
        
        Vlist = list() # List of names of all voltage sources
        elements = list() # Whole block of circuit definition stored as a list
        ele_dict = dict() # Dictionary used to map keys (taken as element names) with class instances
        count = 0 # Stores the row number assigned to the particular node/current branch equation

        for line in lines[start+1:end]:
            eleList = line.split('#')[0].split() 
            elements.append(eleList)
            if eleList[0][0] in ['V', 'H', 'E']:
                Vlist.append(eleList[0])

            # The next two 'if' blocks are appending unique node names and corresponding row numbers in the Table dictionary
            if (eleList[1] not in Table.keys()):
                count += 1
                Table[eleList[1]] = count
            
            if (eleList[2] not in Table.keys()):
                count += 1
                Table[eleList[2]] = count
        
        # Appending voltage source names and corresponding row numbers to Table dictionary
        for V in Vlist:
            count += 1
            Table[V] = count
        
        M = np.zeros((len(Table), len(Table)), dtype=complex) # INCIDENCE MATRIX
        b = np.zeros(len(Table), dtype=complex) # SOURCE MATRIX
        np_elements = np.array([element[0] for element in elements],) # numpy array of element names to check for duplication

        """Check for repetition of element names (Spice does not allow it)"""
        if len(np_elements) != len(np.unique(np_elements)):
            print("Invalid netlist file!\nDuplicate element declaration detected")
            exit()

        for element in elements:
            if element[0][0] in ['R', 'L', 'C']: # Passive elements
                # creating a class instance and mapping the element name (element[0]) to it
                ele_dict[element[0]] = RLC(element[0], element[1], element[2], float(element[3]),)
                RLC.applyToMatrix(ele_dict[element[0]], M) # adding MNA stamp to the matrix

            elif element[0][0] in ['V', 'I']: # Independent sources
                ele_dict[element[0]] = IndSrcs(element[0], element[1], element[2], 
                [float(element[4]) if element[3]=='dc' else float(element[4])/2][0], [float(element[5]) if element[3]=='ac' else 0][0],)
                IndSrcs.applyToMatrix(ele_dict[element[0]], M, b)
            
            elif element[0][0] in ['H', 'F']: # Current-controlled sources
                VnodeIn, VnodeOut = [[ele[1], ele[2]] for ele in elements if ele[0] == element[3]][0]
                ele_dict[element[0]] = ICtrlSrcs(element[0], element[1], element[2], VnodeIn, VnodeOut, 
                element[3], float(element[4]),) 
                ICtrlSrcs.applyToMatrix(ele_dict[element[0]], M)

            elif element[0][0] in ['E', 'G']: # Voltage-controlled sources
                ele_dict[element[0]] = VCtrlSrcs(element[0], element[1], element[2],  
                element[3], element[4], float(element[5]),)
                VCtrlSrcs.applyToMatrix(ele_dict[element[0]], M)
 
        """ Since only n-1 nodal equations are required but I have accepted n, I am 
        removing GND nodal equation so that the solver doesn't encounter a singular matrix"""
        x = np.linalg.solve(M[1:, 1:], b[1:]) 
        x = iter(np.insert(x, 0, 0)) 
        # iter() returns an iterator object and the subsequent elements can be retrieved
        # by calling function next(<iterator>)
        # insert function inserts <arg2> at index <arg1> of array <arg0>
 
        print("\nTable of nodes/branches: ", Table, "\n")
        for row in sorted(Table.values()):
            node_name = get_key(row, Table)
            if w != 0: # in case of AC circuit, report answers in phasors
                if node_name[0] in ['V', 'H', 'E']: # Voltage sources 
                    print("Current through {}: ({c.real:.3e} + {c.imag:.3e}j) A".format(node_name, c = next(x)))
                else:
                    print("Voltage at node {}: ({c.real:.3e} + {c.imag:.3e}j) V".format(node_name, c = next(x)))
            else:
                if node_name[0] in ['V', 'H', 'E']: # Voltage sources 
                    print("Current through {}: {c.real:.3e} A".format(node_name, c = next(x)))
                else:
                    print("Voltage at node {}: {c.real:.3e} V".format(node_name, c = next(x)))

        print("\n")

except IOError:
    # To catch wrong file name error
    print('Invalid file')
    exit()
except:
    # Handles all other exceptions
    print('ERROR: Input file doesn\'t follow standard netlist file format')
    exit()
