'''
SPICE - PART 1
-Reads input netlist file to construct the circuit
-Gives a brief description, along with the tokens in reverse order

Jan 26, 2022

By Arjun R
EE20B016
'''

from sys import argv, exit

"""
Recommended to use constant variables than hard-coding them everywhere.
"""
CIRCUIT = '.circuit'
END = '.end'

"""
It's a good practice to check if the user has given required and only the required inputs.
"""
if len(argv) != 2:
    print('\nERROR: Invalid number of arguments!')
    print('Expected usage: {} <inputfile>'.format(argv[0]))
    exit()

"""
Check whether the input file is a valid netlist file
"""
if (not argv[1].endswith(".netlist")):
    print("ERROR: Wrong file type!")
    exit()

def analyseCkt(elementType, elementList):
    """
    This function analyses each branch and its tokens and prints the details
    of each element present
    
    elementType: element type like any of R, L, C, V, I, or
        current-controlled or voltage-controlled sources
    elementList: List of tokens signifying each branch
    """
    ls = elementList
    if elementType in ['R', 'L', 'C', 'V', 'I']:
        #Type 1 Branch (4 tokens)
        eleName = ls[0]
        node1 = ls[1]; node2 = ls[2]
        val = ls[3]
        if elementType == 'R':
            print('ELEMENT: Resistor')
        elif elementType == 'L':
            print('ELEMENT: Inductor')
        elif elementType == 'C':
            print('ELEMENT: Capacitor')
        elif elementType == 'V':
            print('ELEMENT: Voltage Source')
        elif elementType == 'I':
            print('ELEMENT: Current Source')
        print('Name: {}, Value: {}\nFrom Node: {}, To Node: {}\n'
            .format(eleName, val, node1, node2)) 
        output = ' '.join([val, node2, node1, eleName])
        return output # returns each line to be printed in reverse

    elif elementType in ['E', 'G']: #Voltage-controlled sources
        #Type 2 Branch (6 tokens)
        eleName = ls[0]
        node1 = ls[1]; node2 = ls[2]; node3 = ls[3]; node4 = ls[4]
        val = ls[5]
        if elementType == 'E':
            print('ELEMENT: Voltage Controlled Voltage Source')
        elif elementType == 'G':
            print('ELEMENT: Voltage Controlled Current Source')
        print('Name: {}, Value: {}\nFrom Node: {}, To Node: {}'
            .format(eleName, val, node1, node2))
        print('Controlling voltage nodes: From Node: {}, To Node: {}\n'
            .format(node3, node4)) 
        output = ' '.join([val, node4, node3, node2, node1, eleName])
        return output

    elif elementType in ['H', 'F']: #Current-controlled sources
        #Type 3 Branch (5 tokens)
        eleName = ls[0]
        node1 = ls[1]; node2 = ls[2]
        Vname = ls[3]
        val = ls[4]
        if elementType == 'H':
            print('ELEMENT: Current Controlled Voltage Source')
        elif elementType == 'F':
            print('ELEMENT: Current Controlled Current Source')
        print('Name: {}, Value: {}\nFrom Node: {}, To Node: {}'
            .format(eleName, val, node1, node2))
        print('Controlling current through voltage source {}\n'
            .format(Vname))
        output = ' '.join([val, Vname, node2, node1, eleName])
        return output

"""
The user might input a wrong file name by mistake.
In this case, the open function will throw an IOError.
"""
try:
    with open(argv[1]) as f:
        lines = f.readlines()           # Reads the entire file as a list of lines
        start = -1; end = -2
        for line in lines:              # extracting circuit definition start and end lines
            if CIRCUIT == line[:len(CIRCUIT)]:
                start = lines.index(line)
            elif END == line[:len(END)]:
                end = lines.index(line)
                break
        if start >= end:                # validating circuit block
            print('Invalid circuit definition')
            exit(0)
        final = ''                      # initializing string for final reversed output
        for line in reversed(lines[start+1:end]):
            eleList = line.split('#')[0].split() 
            eleType = eleList[0][0].upper()    # extracts first letter of each line
                                               # also checks for lower case
            final = final + analyseCkt(eleType, eleList) + '\n'

        print('\n'+final)

except IOError:
    # To catch wrong file name error
    print('ERROR: Invalid file')
    exit()

except:
    # Handles all other exceptions
    print('ERROR: Input file doesn\'t follow standard netlist file format')
    exit()

