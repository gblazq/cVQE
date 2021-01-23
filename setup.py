import setuptools

long_description = """cVQE is a collection of classes to run Variational Quantum Eigensolvers
with quadratic Hamiltonians on compressed space.
It provides variational forms, converters and initial states to be used with Qiskit Aqua's
algorithms."""

requirements = [
    "qiskit>=0.23.0",
    "numpy>=1.17"
]

setuptools.setup(
    name="cVQE",
    version="0.0.2",
    author="Guillermo Blázquez",
    description="A package to run Variational Quantum Eigensolvers on compressed space in Qiskit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gblazq/cVQE",
    packages=setuptools.find_packages(),
    license='Apache-2.0',
    classifiers=(
        "Environment :: Console",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering"
    ),
    python_requires='>=3.6',
    install_requires=requirements
)
