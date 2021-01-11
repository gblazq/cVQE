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

from qiskit.aqua.operators.converters import ConverterBase
from qiskit.aqua.operators.list_ops import TensoredOp, SummedOp
from qiskit.aqua.operators.primitive_ops import PauliOp

class TensoredOpDistributor(ConverterBase):
    """
    A converter that applies the distributive property of tensor products.

    Only works with operators that consist on TensoredOps, SummedOps or PauliOps.

    E.g. if op = (I - Z)^(X + 2*Y), the converter will return 
    (I^X) + 2(I^Z) - (Z^X) - 2(Z^Y)
    """
    #TODO: check coefficients
    def convert(self, operator):
        """
        Apply the distributive property to TensoredOps. If the operator is a SummedOp, apply
        it to each of the summands. If it's a PauliOp, return the operator.

        Args:
            operator: the Operator to convert.

        Returns:
            the converted Operator

        Raises:
            TypeError: if the operator is not a TensoredOp, SummedOp or PauliOp
        """
        if isinstance(operator, TensoredOp):
            return reduce(self._convert, operator).reduce()
        elif isinstance(operator, SummedOp):
            return SummedOp([self.convert(op) for op in operator], coeff=operator.coeff).reduce()
        elif isinstance(operator, PauliOp):
            return operator
        else:
            raise TypeError('TensoredOpDistributor can only distribute TensoredOps, SummedOps or PauliOps')

    def _convert(self, op1, op2):
        """
        Distribute the tensor product over two operands in a TensoredOp

        Args:
            op1: the first operator in the tensor product
            op2: the second operator in the tensor product
        
        Returns:
            the result of op1^op2, distributing the tensor product if possible
        
        Raises:
            TypeError: if any of the two operators is not a TensoredOp, SummedOp or PauliOp
        """
        if isinstance(op1, PauliOp) and isinstance(op2, PauliOp):
            return op1^op2
        elif (isinstance(op1, PauliOp) or isinstance(op1, TensoredOp)) and isinstance(op2, SummedOp):
            return SummedOp([self._convert(op1, o) for o in op2.oplist], coeff=op2.coeff)
        elif isinstance(op1, PauliOp) and isinstance(op2, TensoredOp):
            return self._convert(op1, self.convert(op2))
        elif isinstance(op1, SummedOp) and (isinstance(op2, PauliOp) or isinstance(op2, TensoredOp)):
            return SummedOp([self._convert(o, op2) for o in op1.oplist], coeff=op1.coeff)
        elif isinstance(op1, SummedOp) and isinstance(op2, SummedOp):
            return SummedOp([self._convert(o1, o2) for o1 in op1.oplist for o2 in op2.oplist])
        elif isinstance(op1, TensoredOp) and isinstance(op2, PauliOp):
            return self._convert(self.convert(op1), op2)
        elif isinstance(op1, TensoredOp) and isinstance(op2, TensoredOp):
            return self._convert(self.convert(op1), self.convert(op2))
        else:
            raise TypeError('TensoredOpDistributor can only distribute operators consisting on PauliOps, SummedOps or TensoredOps')
