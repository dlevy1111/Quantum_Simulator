import numpy as np
from sympy import *


class Quantum_Circuit:
    def __init__(self, num_qubits):
        self.num_qubits = num_qubits # how many qubits does our quantum circuit operate on?
        self.matrix = np.eye(2**self.num_qubits, dtype=complex)
    
    def h(self, qubits): # hadamard operation on a qubit
        if type(qubits) == int: # support for stacked operations
            qubits = [qubits]
        operator = None
        H = np.matrix('0.707, 0.707; 0.707, -0.707') # hadamard matrix
        I = np.matrix('1, 0; 0, 1') # identity matrix
        for i in range(0,self.num_qubits):
            if i == 0:
                if i in qubits: # gate layer operation
                    operator = H
                else:
                    operator = I
            else:
                if i in qubits:
                    operator = np.kron(H, operator)
                else:
                    operator = np.kron(I, operator)
        self.matrix = operator @ self.matrix

    def x(self, qubits): # pauli-x gate
        if type(qubits) == int: # support for stacked operations
            qubits = [qubits]
        operator = None
        X = np.matrix('0, 1; 1, 0') # x matrix
        I = np.matrix('1, 0; 0, 1') # identity matrix
        for i in range(0,self.num_qubits):
            if i == 0:
                if qubits == i:
                    operator = X
                else:
                    operator = I
            else:
                if qubits == i:
                    operator = np.kron(X, operator)
                else:
                    operator = np.kron(I, operator)
        self.matrix = operator @ self.matrix
    
    def y(self, qubits): # y operation
        if type(qubits) == int: # support for stacked operations
            qubits = [qubits]
        operator = None
        Y = np.matrix("0, -1j; 1j, 0") # y matrix
        I = np.matrix("1, 0; 0, 1") # identity matrix
        for i in range(0,self.num_qubits):
            if i == 0:
                if i == qubits:
                    operator = Y
                else:
                    operator = I
            else:
                if i == qubits:
                    operator = np.kron(Y, operator)
                else:
                    operator = np.kron(I, operator)
        self.matrix = operator @ self.matrix
    
    def z(self, qubits): # z operator
        if type(qubits) == int: # support for stacked operations
            qubits = [qubits]
        operator = None
        Z = np.matrix("1, 0; 0, -1") # z matrix
        I = np.matrix("1, 0; 0, 1") # identity matrix
        for i in range(0, self.num_qubits): 
            if i == 0:
                if i == qubits:
                    operator = Z
                else:
                    operator = I # gate layer operation
            else:
                if i == qubits: 
                    operator = np.kron(Z, operator)
                else:
                    operator = np.kron(I, operator)
        self.matrix = operator @ self.matrix
    
    # def s(self, qubit):
        
init_printing()
n_qubits = 2
a = Quantum_Circuit(n_qubits)
# print(a.matrix)
n = np.array([1, 0])
statevector = np.array([1,0], dtype=complex)
for i in range(1,n_qubits):
    statevector = np.kron(n, statevector)
a.h(0)
a.y([0])
# a.h(0)



print(a.matrix)
print(a.matrix @ statevector)
