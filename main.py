import numpy as np



class Quantum_Circuit:
    def __init__(self, num_qubits):
        self.num_qubits = num_qubits # how many qubits does our quantum circuit operate on?
        self.matrix = np.eye(2**self.num_qubits)
    
    def h(self, qubit): # hadamard operation on a qubit
        a = 1

    def x(self, qubit):
        operator = None
        X = np.matrix('0, 1; 1, 0')
        I = np.matrix('1, 0; 0, 1')
        for i in range(0,self.num_qubits):
            if i == 0:
                if qubit == i:
                    operator = X
                else:
                    operator = I
            else:
                if qubit == i:
                    operator = np.kron(X, operator)
                else:
                    operator = np.kron(I, operator)
        self.matrix = self.matrix @ operator
        


a = Quantum_Circuit(2)
# print(a.matrix)
a.x(1) # doesn't work
a.x(0)
print(a.matrix)
# print(a.statevector[0], a.statevector[])
