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

from functools import reduce
from itertools import islice, product
from operator import xor

from qiskit.aqua.operators import OperatorBase
from qiskit.aqua.operators.converters import ConverterBase
from qiskit.aqua.operators.list_ops import SummedOp, ListOp
from qiskit.aqua.operators.operator_globals import I, X, Y, Z
from qiskit.aqua.operators.primitive_ops import PauliOp
import numpy as np

from .zero_coeff_nuller import ZeroCoeffNuller
from .tensoredop_distributor import TensoredOpDistributor


class QuadraticOperatorReducer(ConverterBase):
    """
    Compress a quadratic operator, given as a qubit operator on n qubits, to log(n) + 1 qubits.

    A quadratic operator in the Majorana operators can only have Z, XX, XY, YX, and YY terms
    when written as a qubit operator.

    The converter can only compress the operator if it is homogeneous (i.e. all `i` terms have the same
    coefficient, where `i` is one of the 5 valid terms).

    Note that the compression of any such operator with 2-qubit terms will result in an operator
    that has an exponential number of terms in the Pauli basis, but the operation can be done in linear space

    Args:
        density: whether to return the original operator's density or the original operator.
            Defaults to True
    """
    def __init__(self, density=True):
        super().__init__()
        self._zero_coeff_nuller = ZeroCoeffNuller(traverse=True)
        self._tensored_op_distributor = TensoredOpDistributor()
        self._density = density

    def convert(self, operator: OperatorBase):
        """
        Convert the original operator to the compressed one as tensors of Pauli operators.

        Args:
            operator: the operator to convert
        
        Returns:
            the compressed operator
        """
        coeffs = self._find_coefficients(operator)
        z, xx, xy, yx, yy = coeffs.get('Z', 0), coeffs.get('XX', 0), coeffs.get('XY', 0), coeffs.get('YX', 0), coeffs.get('YY', 0)

        h = -0.5j*z*(I^Y) + 0.25j*(yx-xy)*(Y^I) + 0.25j*(yx+xy)*(Y^Z) + 0.25j*(yy+xx)*(X^Y) + 0.25j*(yy-xx)*(Y^X)
        h = self._zero_coeff_nuller.convert(h)

        ops = {'I': I,
                'X': X,
                'Y': Y,
                'Z': Z}
        signs = self._sign_generator()
        for i in range(2, int(np.log2(self._num_qubits) + 1)):
            h = I^h

            xx_yy = sum(map(lambda s: (2*s[1][0]*yy - 2*s[1][1]*xx)/(2**(i+2))  * reduce(xor, [ops[op] for op in s[0]]),
                        list(islice(signs, 2**i))))
            xy_yx = sum(map(lambda s: (2*s[1][0]*yx - 2*s[1][1]*xy)/(2**(i+2))  * reduce(xor, [ops[op] for op in s[0]]),
                        list(islice(signs, 2**i))))
            h += xx_yy + xy_yx
        
        h = self._tensored_op_distributor.convert(self._zero_coeff_nuller.convert(h))
        h = (2j*h).reduce()

        if self._density:
            return h
        else:
            return h * operator.num_qubits

    def _find_coefficients(self, operator: OperatorBase):
        """
        Find the coefficients of the Z, XX, XY, YX and YY terms in the quadratic operators.

        Args:
            operator: the operator to find the coefficients

        Returns:
            a dict where the keys are strings identifying the 5 possible terms and the values
            are the coefficient of each term
        
        Raises:
            TypeError: if the operator is not a PauliOp or SummedOp of PauliOps
            ValueError: if the operator has different coefficients for a group of similar terms
        """
        operator = operator.reduce()

        operator_coefficients = {}
        num_operators = {}

        if isinstance(operator, PauliOp):
            operator = [operator]
        elif isinstance(operator, SummedOp) and all([isinstance(op, PauliOp) for op in operator]):
            pass
        else:
            raise TypeError('The operator must be a single PauliOp or a SummedOp of PauliOps')

        self._num_qubits = operator[0].num_qubits
        for op in operator:
            op_str = self._find_operator(op)

            if operator_coefficients.get(op_str) is None:
                operator_coefficients[op_str] = (op.coeff, str(op.primitive))
                num_operators[op_str] = 1
            else:
                num_operators[op_str] += 1

            if operator_coefficients[op_str][0] != op.coeff:
                raise ValueError(f'PauliOps {operator_coefficients[op_str][1]} and {str(op.primitive)} don\'t have the same coefficient')

        if num_operators.get('Z') != self._num_qubits and num_operators.get('Z') is not None:
            raise ValueError('The operator must have all possible Z terms or none')
        for op in ['XX', 'XY', 'YX', 'YY']:
            if num_operators.get(op) != (self._num_qubits - 1) and num_operators.get(op) is not None:
                raise ValueError(f'The operator must have all possible {op} terms or none')

        return {op: operator_coefficients[op][0] for op in operator_coefficients.keys()}
    
    @staticmethod
    def _find_operator(operator):
        """
        Detect which of the 5 possible terms (Z, XX, XY, YX or YY) a PauliOp is

        Args:
            operator: the operator to find which terms is
        
        Returns:
            a string identifying the operator
        
        Raises:
            ValueError: if the operator is any other tensor product of Paulis
        """
        x, z = operator.primitive.x.nonzero()[0], operator.primitive.z.nonzero()[0]
        if x.size == 0 and z.size == 1:
            return 'Z'
        elif z.size == 0 and x.size == 2 and x[1] == (x[0] + 1):
            return 'XX'
        elif x.size == 2 and z.size == 2 and x[1] == (x[0] + 1) and (x == z).all():
            return 'YY'
        elif x.size == 2 and x[1] == (x[0] + 1) and z.size == 1:
            if z[0] == x[0]:
                return 'XY'
            elif z[0] == x[1]:
                return 'YX'
        raise ValueError('The operator must only contain Z, XX, YY, XY or YX terms')

    @staticmethod
    def _sign_generator():
        """
        A generator for the signs in the Pauli decomposition of the compressed operator
        """
        last_signs = {'XX': [np.array([1,1,1,1], dtype=np.complex)]*2,
                      'XY': [np.array([1j,-1j,1j,-1j], dtype=np.complex)]*2,
                      'YX': [np.array([1j,1j,-1j,-1j], dtype=np.complex)]*2,
                      'YY': [np.array([-1,1,1,-1], dtype=np.complex)]*2}
        n = 3
        while True:
            new_last_signs = {}
            for p in product('XY', repeat=n):
                p = ''.join(p)
                new_signs = last_signs[p[1:]][0].copy()
                signs = np.concatenate([new_signs[2:], new_signs[:2]])
                if p[0] == 'Y':
                    new_signs[:2] = 1j*new_signs[:2]
                    new_signs[2:] = -1j*new_signs[2:]
                    signs[:2] = 1j*signs[:2]
                    signs[2:] = -1j*signs[2:]
                new_last_signs[p] = [new_signs, signs]
                
                if p.count('Y') % 2 == 1:
                    yield p, signs.copy()
            
            for p in last_signs.keys():
                if p.count('Y') % 2 == 0:
                    continue

                signs = last_signs[p][1][1:3]
                yield p + 'I', np.array([signs[0], signs[0], signs[1], signs[1]])
                yield p + 'Z', np.array([signs[0], -signs[0], signs[1], -signs[1]])
            
            last_signs = new_last_signs
            n += 1
