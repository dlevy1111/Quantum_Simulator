# starting over
# this time, the simulator will do the following:
# parallel gates
# 2-qubit gates
# a visual section so that you know what your circuit looks like

# my plan:
# instead of making the circuit gate by gate, the simulator will make it column bu column
# this will be done by the user, who will specify all parallel gates, as opposed to specifying gates serially

import numpy as np
from sympy import *

class Quantum_Circuit:
    def __init__(self, num_qubits) -> None:
        self.num_qubits = num_qubits # how many qubits does our circuit operate on?
        self.matrix = np.eye(2**num_qubits, dtype=complex)
        self.visual_circuit = "" # the visual part
    
    def make_column(self,column : list) -> None: # column is the list of tuples. each tuple will include the following information: gate type, qubit-that-is-being-acted-on, and if necessary, information that is pertinent to the gate (rotations, controls)
        # example input: column = [("h", 0), ("x",1), ("x",2)]
        # the user won't have to specify qubits that don't have gates operating on them (i.e. identity operations)
        # self.graphic_designer(column) # a method that creates the picture
        
        

        if len(column) != self.num_qubits: # this loop only occurs if the circuit doesn't have gates on all qubits.
            all_qubits = set(range(0,self.num_qubits)) # using set operations
            acting_qubits = set([column[i][1] for i in range(0,len(column))])
            unaltered_qubits = all_qubits - acting_qubits
            for i in range(0,self.num_qubits):
                if i in unaltered_qubits:
                    column.insert(i, ("i", i)) # should be all fixed now
        print(column)
        if 0 not in column[0]:
            column_matrix = self._logic("i")
        else:
            column_matrix = self._logic(column[0][0])
        for i in range(1,self.num_qubits): # skip the first qubit because we need to make the column matrix at some point (before this loop); (if num_qubits == 1, the loop doesn't happen anyway.)
            gate_information = column[i]
            if i not in gate_information: # qubit "i" will not be acted on in this column 
                column_matrix = np.kron(self._logic("i"), column_matrix)
            else:
                column_matrix = np.kron(self._logic(gate_information[0]), column_matrix)
        self.matrix = column_matrix @ self.matrix

    def print_matrix(self):
        init_printing()
        matrix_out = Matrix(self.matrix)
        pprint(matrix_out)

    def _logic(self, strnput) -> np.ndarray: # strnpt is a portmanteau of string input
        all_outputs = { # created to make it easier to return functions, as well as keep the code clean
            "i": self._i(),
            "h": self._h(),
            "x": self._x(),
            "y": self._y(),
            "z": self._z()
            # with more room for additions
        }
        return all_outputs[strnput]

    def _i(self) -> np.ndarray:
        I = np.matrix('1, 0; 0, 1') # identity matrix
        return I
    
    def _h(self) -> np.ndarray:
        H = np.matrix('0.707, 0.707; 0.707, -0.707') # hadamard matrix
        return H
    
    def _x(self) -> np.ndarray: 
        X = np.matrix('0, 1; 1, 0') # pauli x matrix
        return X
    
    def _y(self) -> np.ndarray:
        Y = np.matrix("0, -1j; 1j, 0") # pauli y matrix
        return Y

    def _z(self) -> np.ndarray:
        Z = np.matrix("1, 0; 0, -1") # pauli z matrix
        return Z
    
    def direct_sum(self, a : np.matrix, b : np.matrix) -> np.matrix:
        c = np.kron(a, np.eye(b.shape[0])) + np.kron(np.eye(a.shape[0]), b) # thanks chatgpt
        return c
    

num_qubits = 3

qc = Quantum_Circuit(num_qubits)
# qc.make_column([("h", 0), ("x", 1), ("z", 2)])
qc.make_column([("h", 2)])


n = np.array([1, 0])
statevector = np.array([1,0], dtype=complex)
for i in range(1,num_qubits):
    statevector = np.kron(n, statevector)

qc.print_matrix()
pprint(Matrix(qc.matrix @ statevector))
