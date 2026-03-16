# Quantum Gate Decomposition Validator

This repository contains the supplementary software artifact for the paper **"Conjuntos de Portas Nativos e Decomposições de Portas: Um Estudo sobre Processadores Quânticos Atuais"**. 

It provides a symbolic validation script written in Python that mathematically proves the algebraic equivalence between universal quantum logic gates (such as $U3$, $CNOT$, and $SWAP$) and their hardware-specific native decompositions across three major quantum computing architectures: **IBM**, **Rigetti**, and **IonQ**.

## Overview

Physical Quantum Processing Units (QPUs) do not implement universal gates natively. Compilers must translate these abstract operations into sequences of hardware-specific pulses (native gates). 

This script uses **SymPy** to symbolically compute the matrix multiplication of the decomposed circuits and compares them against the ideal matrices. It accounts for global phase differences by verifying if the product of the decomposed matrix and the adjoint of the ideal matrix results in a homogeneous diagonal matrix (proportional to the Identity).

### Supported Architectures
* **IBM (Superconducting):** $R_Z$, $SX$, $X$, $CZ$, $R_{ZZ}$, $ECR$
* **Rigetti (Superconducting):** $R_X$, $R_Z$, $iSWAP$
* **IonQ (Trapped Ions):** $GPi$, $GPi2$, $R_Z$, $MS$ (Mølmer–Sørensen)

## Prerequisites

To run this validator, you need:
* Python 3.8 or higher
* `sympy` library

## Installation

First, clone the repository and navigate to the project directory.

### Option 1: Using `uv` (Recommended)
If you use the `uv` package manager, you can quickly set up the environment and run the script:

```bash
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt
```

### Option 2: Using `pip`
If you use the `pip` package manager, you can set up the environment and install dependencies as follows:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

Once the dependencies are installed and the virtual environment is activated, simply run the Python script:

```bash
python gate_decompositions.py
```

## Understanding the Output

The script evaluates each decomposition and outputs the validation status in the terminal:

* `[PASS]`: The decomposed matrix is mathematically equivalent to the ideal gate (ignoring global phases). The topology alignment (e.g., `CNOT_01` vs `CNOT_10`) is also specified.
* `[FAIL]`: The decomposition is not mathematically equivalent to the target matrix.

### Output Example

```bash
=============================================
 IBM ARCHITECTURE
=============================================
[PASS] IBM U3   (2x SX) (Alignment: CNOT_01)
[PASS] IBM CNOT (1x CZ) (Alignment: CNOT_01)
[PASS] IBM CNOT (1x Rzz) (Alignment: CNOT_01)
[PASS] IBM CNOT (1x ECR) (Alignment: CNOT_01)
[PASS] IBM SWAP (3x ECR CNOTs) (Alignment: CNOT_01)
```

## Authors

* **Ruan Luiz Molgero Lopes** - *Universidade Federal de Santa Catarina (UFSC)*
* **Eduardo Lussi** - *Universidade Federal de Santa Catarina (UFSC)*
* **Evandro Chagas Ribeiro da Rosa** - *Universidade Federal de Santa Catarina (UFSC)*
* **Jerusa Marchi** - *Universidade Federal de Santa Catarina (UFSC)*

*Developed as part of the Quantum Computing Group (GCQ) at UFSC.*