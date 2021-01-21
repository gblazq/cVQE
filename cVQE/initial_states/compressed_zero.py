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

from qiskit import QuantumCircuit, QuantumRegister
from qiskit.aqua import AquaError
from qiskit.aqua.components.initial_states import InitialState
from qiskit.aqua.utils.validation import validate_min
import numpy as np

class CompressedZero(InitialState):
    r"""
    The purified, compressed |0>^n state, i.e. the state

    \frac{1}{\sqrt{2n}}\sum_{j=0}^{2^{n - 1} - 1}{\vert j \rangle \vert 0 \rangle \vert j \rangle +i\vert j\rangle \vert 1 \rangle \vert j\rangle}
    """
    def __init__(self, num_qubits):
        """
        Args:
            num_qubits: number of qubits of the state. Important: this is the final
                state number of qubits (i.e. after compression and purification), so
                it will be the compressed version of the 2^((num_qubits - 1) / 2) |0> state.
        """
        super().__init__()
        validate_min('num_qubits', num_qubits, 3)

        if num_qubits % 2 == 0:
            raise ValueError('num_qubits must be odd to have a valid compressed state')

        self._num_qubits = num_qubits

    def construct_circuit(self, mode='circuit', register=None):
        """
        Construct the statevector or circuit of desired initial state.

        Args:
            mode: `vector` or `circuit`. The `vector` mode produces the vector,
                while the `circuit` mode constructs the quantum circuit.
            register: register for circuit construction.

        Returns:
            QuantumCircuit or numpy.ndarray: quantum circuit or statevector.

        Raises:
            ValueError: when mode is not 'vector' or 'circuit'.
        """
        dim = (self._num_qubits + 1)//2

        if mode == 'vector':
            state = np.zeros((2**(2*dim - 1)), dtype=np.complex)
            for i in range(2**(dim-1)):
                j1 = i << dim
                j2 = i
                state[j1 + j2] = 1
                state[j1 + j2 + (1 << (dim-1))] = 1j
            return state/np.linalg.norm(state)
        elif mode == 'circuit':
            if register is None:
                register = QuantumRegister(2*dim - 1, name='q')
            quantum_circuit = QuantumCircuit(register)

            quantum_circuit.h(list(range(dim - 1, 2*dim - 1)))
            quantum_circuit.s(dim - 1)
            quantum_circuit.cx(list(range(dim, 2*dim - 1)), list(range(dim - 1)))

            return quantum_circuit
        else:
            raise AquaError('Mode should be either "vector" or "circuit"')
