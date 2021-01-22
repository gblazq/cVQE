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

from itertools import product
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
    """
    def __init__(self, num_qubits: int, method: str = 'CRY',
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

        self._son_num_qubits = num_qubits
        # TODO: make num_qubits a property
        if self._initial_state is not None and count_extra_qubits:
            self._num_qubits = self._initial_state._num_qubits
        else:
            self._num_qubits = num_qubits

        self._son_registers = son_registers or list(range(num_qubits))
        self._num_parameters = (2**num_qubits) * (2**num_qubits - 1) // 2
        self._bounds = [(-np.pi, np.pi) for _ in range(self._num_parameters)]

    def construct_circuit(self, parameters: Union[List[float], np.ndarray],
                          q: Optional[QuantumRegister] = None,
                          derivative: int = None,
                          shift: float = None) -> QuantumCircuit:
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

        derivative_ancilla = 1 if derivative is not None else 0

        if q is None:
            if self._initial_state is not None:
                q = QuantumRegister(self._initial_state._num_qubits + derivative_ancilla, name='q')
            else:
                q = QuantumRegister(self._num_qubits + derivative_ancilla, name='q')
        
        if self._initial_state is not None:
            circuit = self._initial_state.construct_circuit('circuit', q)
        else:
            circuit = QuantumCircuit(q)
        
        # For 1 qubit, the special orthogonal rotation is just an RY gate
        # TODO: derivative
        if self.num_qubits == 1:
            circuit.ry(parameters[0], self._son_registers)
            return circuit

        ansatz =  self._circuit_constructors[self._method](parameters, derivative, shift)
        
        return circuit.compose(ansatz, self._son_registers + derivative_ancilla*[len(q) - 1])

    def _construct_cry_circuit(self, parameters, derivative=None, shift=None):
        """
        Construct a quantum circuit that implements an SO(n**2) unitary operation
        on n qubits using multicontrolled RY gates to implement the Givens rotations

        Args:
            parameters: circuit parameters
            derivative: the index of the parameter whose derivative circuit we want
                to compute. If None, return the special unitary circuit
            shift: the shift in the derivative circuit
        
        Returns:
            QuantumCircuit: a quantum circuit with given `parameters` implementing the
                SO(n**2) unitary using multicontrolled RY gates
        """
        if derivative is not None:
            qc = QuantumCircuit(self._son_num_qubits + 1)
            qc.h(self._son_num_qubits)
            qc.x(self._son_num_qubits)
        else:
            qc = QuantumCircuit(self._son_num_qubits)
        
        def param_order():
            m = 2**self._son_num_qubits
            for i in range(m-1):
                for j in range(1+i, m):
                    yield i, j
        
        for param, (row, col) in enumerate(param_order()):
            qubit_pairs = {0: [], 1: [], 2: [], 3: []}
            for qubit in range(self._son_num_qubits):
                row_value = row & (1 << qubit) > 0
                col_value = col & (1 << qubit) > 0

                qubit_pairs[(row_value << 1) + col_value].append(qubit)

            controlled, qubit_pairs[1] = qubit_pairs[1][0], qubit_pairs[1][1:]

            if qubit_pairs[0] or qubit_pairs[1]:
                qc.x(qubit_pairs[0] + qubit_pairs[1])
            if qubit_pairs[1] or qubit_pairs[2]:
                qc.cnot(controlled, qubit_pairs[1] + qubit_pairs[2])
            
            if derivative is not None and param == derivative:
                qc.append(RYGate(parameters[param]).control(self._son_num_qubits),
                    list(filter(lambda x: x != controlled, range(self._son_num_qubits + 1))) + [controlled])
                qc.x(self._son_num_qubits)
                qc.append(RYGate(parameters[param] + shift).control(self._son_num_qubits),
                    list(filter(lambda x: x != controlled, range(self._son_num_qubits + 1))) + [controlled])                
                qc.h(self._son_num_qubits)
            else:
                qc.append(RYGate(parameters[param]).control(self._son_num_qubits - 1), 
                    list(filter(lambda x: x != controlled, range(self._son_num_qubits))) + [controlled])
            
            if qubit_pairs[1] or qubit_pairs[2]:
                qc.cnot(controlled, qubit_pairs[1] + qubit_pairs[2])
            if qubit_pairs[0] or qubit_pairs[1]:
                qc.x(qubit_pairs[0] + qubit_pairs[1])
                        
        return qc

    def get_gradient_callable(self, operator, quantum_instance=None):
        """
        Create a callable that evaluates the gradient of the ansatz expectation value 
        w.r.t. the ansatz parameters

        Args:
            operator: the operator whose expectation value's gradient we want to compute
            quantum_instance: the quantum instance to compute the expectation value
        
        Returns:
            callable(parameter_values): Function to compute the gradient
        """
        from qiskit.circuit import ParameterVector
        from qiskit.aqua.operators.list_ops import ListOp, SummedOp
        from qiskit.aqua.operators.operator_globals import Z
        from qiskit.aqua.operators.state_fns import StateFn
        from qiskit.aqua.operators.converters import CircuitSampler

        parameters = ParameterVector('θ', self._num_parameters)
        derivatives = []
        observable_operator = ~StateFn(Z^operator)
        for i in range(len(parameters)):
            derivative = []
            for shift in [1, -1]:
                op = observable_operator @ StateFn(self.construct_circuit(parameters, derivative=i, shift=shift*np.pi/2))
                derivative.append(shift*op/np.sqrt(2))
            derivatives.append(sum(derivative))
        gradient_operator = ListOp(derivatives)

        def gradient_fn(parameter_values):
            p_values_dict = dict(zip(parameters, parameter_values))
            if not quantum_instance:
                gradient = gradient_operator.assign_parameters(p_values_dict)
                gradient = gradient.eval()
            else:
                p_values_dict = {k: [v] for k, v in p_values_dict.items()}
                gradient = CircuitSampler(backend=quantum_instance).convert(
                    gradient_operator, p_values_dict)
                gradient = gradient.eval()[0]
            
            return np.array(list(map(np.real, gradient)))

        return gradient_fn
