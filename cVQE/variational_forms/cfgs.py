# Copyright 2021 Guillermo Blázquez
#
# Licensed under the Apache License, Version 2.0 (the "License")
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import List, Optional, Union

import numpy as np
from qiskit import QuantumRegister, QuantumCircuit
from qiskit.aqua.components.initial_states import InitialState
from qiskit.aqua.components.variational_forms import VariationalForm
from qiskit.aqua.utils.validation import validate_min, validate_in_set
from qiskit.circuit.library.standard_gates import RYGate

class CompressedFermionicGaussianState(VariationalForm):
    """
    An ansatz that implements a compressed Fermionic Gaussian state by
    starting with a compressed Fermionic Gaussian state and applying a
    special orthogonal rotation from the SO(n**2) group in n qubits.
    An ansatz that implements a special orthogonal rotation in n qubits, 
    i.e. a member of the SO(n**2) group with 2**n*(2**n-1)/2 parameters
    """
    def __init__(self,
                 num_qubits: int,
                 method: str = 'CRY',
                 initial_state: Optional[InitialState] = None,
                 son_registers: Optional[List[int]] = None,
                 count_extra_qubits: bool = False) -> None:
        """
        Args:
            num_qubits: number of qubits
            method: method to build the SO(n**2) operation. Must be 'CRY'.
                The 'CRY' method builds the variational form using a quantum circuit with multicontrolled RY gates.
            son_registers: register positions to which the ansatz will be applied to. Useful if the initial state has 
                more qubits than n**2
            count_extra_qubits: whether to set the variational form's num_qubits attribute to all the qubits (if True)
                or just the qubits of the special orthogonal rotation (if False, default). Counting only the n qubits
                of the SO rotation is useful if the rest of the qubits are ancillae and you want to use the variational
                form in a VQE algorithm, so that the algorithm accepts operators of just n qubits.
        """
        super().__init__()
        self._circuit_constructors = {'CRY': self._construct_cry_circuit}
        
        validate_min('num_qubits', num_qubits, 1)
        validate_in_set('method', method, self._circuit_constructors.keys())

        self._method = method
        self._initial_state = initial_state

        # TODO: make num_qubits a property
        if self._initial_state is not None and count_extra_qubits:
            self._num_qubits = self._initial_state._num_qubits
        else:
            self._num_qubits = num_qubits

        self._son_registers = son_registers
        self._num_parameters = (2**num_qubits) * (2**num_qubits - 1) // 2
        self._bounds = [(-np.pi, np.pi) for _ in range(self._num_parameters)]

    def construct_circuit(self, parameters: Union[List[float], np.ndarray],
                          q: Optional[QuantumRegister] = None) -> QuantumCircuit:
        """
        Construct the variational form, given its parameters.

        Args:
            parameters: circuit parameters
            q: Quantum Register for the circuit.

        Returns:
            QuantumCircuit: a quantum circuit with given `parameters`

        Raises:
            ValueError: if the number of parameters is incorrect.
        """
        if len(parameters) != self._num_parameters:
            raise ValueError('The number of parameters has to be {}'.format(self._num_parameters))

        if q is None:
            if self._initial_state is not None:
                q = QuantumRegister(self._initial_state._num_qubits, name='q')
            else:
                q = QuantumRegister(self._num_qubits, name='q')
        
        if self._initial_state is not None:
            circuit = self._initial_state.construct_circuit('circuit', q)
        else:
            circuit = QuantumCircuit(q)
        
        # For 1 qubit, the special orthogonal rotation is just an RY gate
        if self.num_qubits == 1:
            circuit.ry(parameters[0], self._son_registers)
            return circuit

        ansatz =  self._circuit_constructors[self._method](parameters)
        
        return circuit.compose(ansatz, self._son_registers)

    def _construct_cry_circuit(self, parameters):
        """
        Construct a quantum circuit that implements an SO(n**2) unitary operation
        on n qubits using multicontrolled RY gates to implement the Givens rotations

        Args:
            parameters: circuit parameters
        
        Returns:
            QuantumCircuit: a quantum circuit with given `parameters` implementing the
                SO(n**2) unitary using multicontrolled RY gates
        """
        qc = QuantumCircuit(self._num_qubits)
        
        def param_order():
            m = 2**self._num_qubits
            for i in range(m-1):
                for j in range(1+i, m):
                    yield i, j
        
        for param, (row, col) in enumerate(param_order()):
            qubit_pairs = {0: [], 1: [], 2: [], 3: []}
            for qubit in range(self._num_qubits):
                row_value = row & (1 << qubit) > 0
                col_value = col & (1 << qubit) > 0

                qubit_pairs[(row_value << 1) + col_value].append(qubit)

            controlled, qubit_pairs[1] = qubit_pairs[1][0], qubit_pairs[1][1:]

            if qubit_pairs[0] or qubit_pairs[1]:
                qc.x(qubit_pairs[0] + qubit_pairs[1])
            if qubit_pairs[1] or qubit_pairs[2]:
                qc.cnot(controlled, qubit_pairs[1] + qubit_pairs[2])

            qc.append(RYGate(parameters[param]).control(self._num_qubits - 1), 
                    list(filter(lambda x: x != controlled, range(self._num_qubits))) + [controlled])
            
            if qubit_pairs[1] or qubit_pairs[2]:
                qc.cnot(controlled, qubit_pairs[1] + qubit_pairs[2])
            if qubit_pairs[0] or qubit_pairs[1]:
                qc.x(qubit_pairs[0] + qubit_pairs[1])

        return qc
