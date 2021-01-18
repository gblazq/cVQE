# cVQE

cVQE (compressed Variational Quantum Eigensolver) is a collection of classes to run Variational Quantum Eigensolvers
with quadratic Hamiltonians on 2log(n) + 1 qubits.
It provides variational forms, converters and initial states to be used with Qiskit Aqua's
algorithms.

## Installation

You can install `cVQE` with `pip`:

```bash
$ pip install cVQE
```

cVQE requires Qiskit >= 0.23.1.

## Running compressed VQEs

To run a VQE algorithm in compressed space you'll need to compress the Hamiltonian and use an appropriate ansatz. cVQE provides a converter to compress quadratic homogeneous Hamiltonians, and a variational form and initial state to build an ansatz that can explore the whole set of quadratic Hamiltonians eigenstates:

```python
from qiskit.aqua.operators.operator_globals import I, X, Z
from cVQE.operators.converters import QuadraticOperatorReducer
from cVQE.initial_states import CompressedZero
from cVQE.variational_forms import CompressedFermionicGaussianState

# Set the original, fully compressed and final number of qubits 
original_num_qubits = 2
compressed_num_qubits = int(np.log2(n) + 1)
total_num_qubits = 2*compressed_num_qubits - 1

# Build the initial state
initial_state = CompressedZero(total_num_qubits)

# Create the variational form
var_form = CompressedFermionicGaussianState(compressed_num_qubits, initial_state=initial_state)

# Create a quadratic, homogeneous, 2-qubit Hamiltonian
H = (I^Z) + (Z^I) +2*(Y^Y)

# Compress the Hamiltonian
op_reducer = QuadraticOperatorReducer()
h = op_reducer.convert(H)
```

After these steps, you can use the variational form and the compressed operator to run a VQE algorithm. 

For more detailed and runnable examples, including the use of gradients, see the [examples notebook](examples/cVQE_examples.ipynb).

## Compressed space

### Matchgate circuits

A nearest-neighbour matchgate circuit on `n` qubits can either be simulated efficiently in a classical computer [1], or simulated using an exponentially smaller number of qubits, `log(n)`, in a quantum computer [2], provided that we measure an operator that is quadratic in the creation and annihilation operators.

The compression of such a circuit relies on three facts. First, a quadratic Hamiltonian on `n` qubits can be written as

![Image](images/H.png)

where the `x` are Majorana operators and `a` is a real antisymmetric matrix. `i*a` can be understood as a unitary gate acting on `log(n) + 1` qubits. Second, the set of Fermionic Gaussian states, which includes all ground state to quadratic Hamiltonians, is closed under the action of nearest-neighbour matchgate circuits. A Fermionic Gaussian state `ρ` can be completely characterized by its second moments, which are encoded in the covariance matrix

![Image](images/covariance_matrix.png)

This matrix has the same dimensionality as `a`. It can be used to define a density matrix `σ` as

![Image](images/density_matrix.png)

Finally, a matchgate circuit can be encoded in a rotation matrix `R` acting on `σ` such that

![Image](images/exp_value.png)

where 

![Image](images/compressed_operator.png)

If we are able to compress an operator `H`, and starting with any compressed Fermionic Gaussian state `σ`, we can use the previous identity to find the compressed ground state of `H` by means of a VQE algorithm where the ansantz implements a special orthogonal rotation of the qubits.

### Compressed Hamiltonian

The compression of a quadratic Hamiltonian relies on the fact that it can be written as

![Image](images/H.png)

In terms of the Pauli matrices, a quadratic Hamiltonian can only have the nearest-neighbour terms Z, XX, YY, XY and YX; and the compressed matrix `a` takes the form


where `O_{i,j}` is the coefficient of the operator `O` acting on qubits `i` and `j` in the original Hamiltonian. This is a banded matrix which is non-zero in the 3 diagonals up to the main one and the 3 down to it, and zero elsewhere.

In order to use the compressed Hamiltonian in a VQE algorithm, we need to know its decomposition in the basis of tensors of Pauli matrices. In general, this decomposition will have an exponential number of terms, and finding the coefficients by means of the basis matrices would require storing exponentially large matrices. Here we provide a classical algorithm to find the exponentially large decomposition of a homogeneous quadratic Hamiltonian.

First, notice that if the Hamiltonian is homogeneous, we can build `a` iteratively by using the fact that it's a banded matrix:


### Compressed states

TODO

### Building an ansatz

TODO

### References

TODO

## TODO

* Finish README.md.
* Add tests.

## License

[Apache License 2.0](LICENSE.txt)
