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

from qiskit.aqua.operators.converters import ConverterBase
from qiskit.aqua.operators.list_ops import SummedOp

class ZeroCoeffNuller(ConverterBase):
    """
    A converter that nullifies operators with coeff = 0 in SummedOps
    """
    def __init__(self, traverse=True):
        """
        Args:
            traverse: if True, traverses down all elements of nested SummedOps.
                If False, only checks the first level
        """
        super().__init__()
        self._traverse = traverse
    
    def convert(self, operator):
        """
        If the operator is a SummedOp, filter out operators with coeff = 0.
        If not, return 0 if operator.coeff = 0 or return the operator otherwise.

        Args:
            operator: the operator to convert
        
        Returns:
            the converted operator
        """
        if operator.coeff == 0:
            return 0
        elif isinstance(operator, SummedOp):
            if self._traverse:
                return SummedOp(list(map(self.convert, filter(lambda o: o.coeff != 0, operator))), coeff=operator.coeff)
            else:
                return SummedOp(list(filter(lambda o: o.coeff != 0, operator)), coeff=operator.coeff)
        else:
            return operator
