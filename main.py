# starting over
# this time, the simulator will do the following:
# parallel gates ✔
# 2-qubit gates ✔
# n-qubit gates (should not be that hard honestly) (no longer necessary because the user has access to a complete gate set (clifford + T))
# a visual section so that you know what your circuit looks like

# my plan:
# instead of making the circuit gate by gate, the simulator will make it column bu column
# this will be done by the user, who will specify all parallel gates, as opposed to specifying gates serially

import numpy as np
from sympy import *

class Quantum_Circuit:
    def __init__(self, num_qubits, mode = "none") -> None:
        self.num_qubits = num_qubits # how many qubits does our circuit operate on?
        self.matrix = np.eye(2**num_qubits, dtype=complex)
        self.visual_circuit = "" # the visual part
        self.mode = mode # mode means printing like ibmq (with reversed endian) or not.
        self._u_values = np.matrix('1, 0; 0, 1')
    
    def make_column(self,column : list) -> None: # column is the list of tuples. each tuple will include the following information: gate type, qubit-that-is-being-acted-on, and if necessary, information that is pertinent to the gate (rotations, controls)
        # example input: column = [("h", 0), ("x",1), ("x",2)]
        # the user won't have to specify qubits that don't have gates operating on them (i.e. identity operations)
        # self.graphic_designer(column) # a method that creates the picture
        # second example: column = [("cnot", (0, 1)), ("h", 2)]
        
        two_qubit_gates = 0 # 2-qubit gates look like 1-qubit gates to the logic (because the len(column) doesn't reflect the tuple of control and target bits)
        completed_two_qubit_gates = 0 # because of the above note, we need to keep track of whether or not we've calculate a 2-qubit gate into the column matrix
        for gate_information in column:
            if type(gate_information[1]) == tuple:
                two_qubit_gates+=1

        
        if len(column)+two_qubit_gates != self.num_qubits: # this loop only occurs if the circuit doesn't have gates on all qubits.
            all_qubits = set(range(0,self.num_qubits)) # using set operations
            single_gate_qubits = set([column[i][1] for i in range(0,len(column)) if type(column[i][1]) == int])
            multi_gate_qubits = set(sum([list(item[1]) for item in column if type(item[1]) == tuple],[])) # found this solution online but god damn. **********very likely to break with 3-qubit gates
            acting_qubits = single_gate_qubits.union(multi_gate_qubits)
            # print(acting_qubits)
            unaltered_qubits = all_qubits - acting_qubits
            
            for i in range(0,self.num_qubits):
                if i in unaltered_qubits:
                    column.insert(i, ("i", i)) # should be all fixed after this line
        

        # print(column)
        # print(column[0])
        i = 1 # using a while loop so that we can skip by 2 for 2-qubit gates, or by 3 for 3-qubit ones, whatever's necessary
        if type(column[0][1]) == int:
            column_matrix = self._logic(column[0][0])
        else: # starting with a 2-qubit gate
            # print("applying 2qubit gate1")
            column_matrix = self._logic(column[0][0])
            completed_two_qubit_gates+=len(column[0][1])-1
            i+=len(column[0][1])-1
        
        while i < self.num_qubits: # skip the first qubit because we need to make the column matrix at some point (before this loop); (if num_qubits == 1, the loop doesn't happen anyway.)
            gate_information = column[i-completed_two_qubit_gates]
            # print(column_matrix, i, completed_two_qubit_gates)
            # print(gate_information)
            if type(gate_information[1]) == int: # 1 qubit gates are just applied on one gate (n ∈ ℕ)
                column_matrix = np.kron(column_matrix, self._logic(gate_information[0]))
                i+=1
            else: # if it's not a 1 qubit gate then its gotta be at least a 2-qubit one
                # print("applying 2qubit gate2", "len(gate_information[1]) =", len(gate_information[1]))
                column_matrix = np.kron(column_matrix, self._logic(gate_information[0]))
                completed_two_qubit_gates+=len(gate_information[1])-1
                i+=len(gate_information[1])-1 # for a 2 qubit gate we kron with the 4x4 matrix and skip a qubit (because that's how it works yo).
                # *************** to change this to any n-qubit gate just change "i+=2" to "i+=len(gate_information[1])", i think
        # print(column_matrix, i, completed_two_qubit_gates)
        self.matrix = column_matrix @ self.matrix

    def print_matrix(self) -> None:
        init_printing()
        matrix_out = Matrix(self.matrix)
        pprint(matrix_out)
    
    def output_statevector(self) -> np.ndarray:
        statevector = np.append(np.ones(1), np.zeros((2**(self.num_qubits)-1,1)))
        output = self.matrix @ statevector
        if self.mode == "ibmq": # sort by 
            self._ibmq_sort(output)
        return output
    
    def _ibmq_sort(self, values : np.matrix) -> None:
        # we want to swap values based on the inverse of their position in the list (as a binary number)
        array = np.squeeze(np.asarray(values))
        moved_values = []
        for i in range(0, array.shape[0]): # only loop through half of the list because otherwise we'd swap back.
            new_spot = int(format(i, f"0{self.num_qubits}b")[::-1],2) # help from stack overflow
            if i not in moved_values and new_spot not in moved_values:
                array[i], array[new_spot] = array[new_spot], array[i]
                moved_values.append(i)
                moved_values.append(new_spot)
        return np.asmatrix(array)
            
    def _logic(self, strnput) -> np.ndarray: # strnpt is a portmanteau of string input
        all_outputs = { # created to make it easier to return functions, as well as keep the code clean
            "i": self._i(),
            "h": self._h(),
            "x": self._x(),
            "y": self._y(),
            "z": self._z(), # complete gate set formed by H, T, S, and CNOT (Clifford + T)
            "t": self._t(),
            "tdg": self._tdg(),
            "s": self._s(),
            "u": self._u(),
            "cx": self._cx(),
            "swap": self._swap(),
            "cz": self._cz(),
            "ccx": self._ccx()
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
    
    def _t(self) -> np.ndarray:
        T = np.matrix("1, 0; 0, 0.707+0.707j")
        return T
    
    def _tdg(self) -> np.ndarray: # t dagger
        TDG = np.matrix("1, 0; 0, 0.707-0.707j")
        return TDG
    
    def _s(self) -> np.ndarray: # phase gate
        S = np.matrix("1, 0; 0, 1j")
        return S
    
    def _u(self) -> np.ndarray: # arbitrary u gate with associated change_u() elsewhere
        U = self._u_values # starts out as I
        return U # should work for arbitrary size, because the qubits that this gate will act on are not defined here anyway.
    
    def change_u(self, new_values : np.matrix) -> None:
        self._u_values = new_values

    def _cx(self) -> np.ndarray:
        # CX = self._direct_sum(self._i(), self._x()) # works backwards for some reason :( incorrect behavior :(
        CX = np.matrix("1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 0, 1; 0, 0, 1, 0")
        # CX = np.matrix("1, 0, 0, 0; 0, 0, 0, 1; 0, 0, 1, 0; 0, 1, 0, 0")
        return CX
    
    def _swap(self) -> np.ndarray: # seems to work fine??
        # SWAP = self._cnot() @ ( np.kron(self._h(), self._h()) @ self._cnot() @ np.kron(self._h(), self._h()) ) @ self._cnot() # it's a mess but this is the breakdown of the swap matrix
        SWAP = np.matrix("1, 0, 0, 0; 0, 0, 1, 0; 0, 1, 0, 0; 0, 0, 0, 1")
        return SWAP
    
    def _cz(self) -> np.ndarray:
        CZ = np.matrix("1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, -1")
        return CZ
    
    def _ccx(self) -> np.ndarray:
        CCX = np.matrix("1, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 1, 0")
        return CCX

    def _direct_sum(self, a : np.matrix, b : np.matrix) -> np.matrix:
        c = np.block([[a, np.zeros_like(a)], [np.zeros_like(b), b]])
        return c


num_qubits = 4

qc = Quantum_Circuit(num_qubits, "ibmq")
# qc = Quantum_Circuit(num_qubits)


# qc.make_column([("h",0), ("h", 1), ("h",2)])
# qc.make_column([("z",1)])
qc.make_column([("h",0)])
qc.make_column([("cx", (0,1))])
qc.make_column([("ccx", (0, 1, 2))])

# qc.make_column([("h",0), ("h",1), ("h",2), ("h",3), ("h",4)])
# qc.make_column([("cx",(0,1)), ("z",2)])
# qc.make_column([("y",1)])

# invcx = np.matrix("0, 0, 1, 0; 0, 0, 0, 1; 0, 1, 0, 0; 1, 0, 0, 0")
# qc.change_u(invcx)
# qc.make_column([("u", (0,1))])

# qc.make_column([("h",0), ("h",1)])
# qc.make_column([("cz", (0,1))])

statevector = qc.output_statevector()

qc.print_matrix()
pprint(Matrix(statevector))
