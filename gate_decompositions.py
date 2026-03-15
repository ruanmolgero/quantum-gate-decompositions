"""
Quantum Gate Decomposition Validator.

This script uses SymPy to symbolically validate algebraic equivalences
between generic quantum operations (like U3, CNOT, and SWAP) and their
hardware-specific native decompositions across different quantum 
architectures (IBM, Rigetti, and IonQ).

The validation accounts for global phase differences by checking if 
the product of the decomposed matrix and the adjoint of the ideal matrix 
results in a homogeneous diagonal matrix (proportional to the Identity).
"""

import sympy as sp
from sympy import pi, I, exp, cos, sin, Matrix, eye, zeros
from sympy.physics.quantum import TensorProduct as kron

# ==========================================
# 1. Base Matrices
# ==========================================
I_mat = eye(2)
X = Matrix([[0, 1], [1, 0]])
Y = Matrix([[0, -I], [I, 0]])
Z = Matrix([[1, 0], [0, -1]])
H = (X + Z) / sp.sqrt(2)

def Rx(theta):
    """Returns the matrix representation of an X-axis rotation."""
    return cos(theta / 2) * I_mat - I * sin(theta / 2) * X

def Ry(theta):
    """Returns the matrix representation of a Y-axis rotation."""
    return cos(theta / 2) * I_mat - I * sin(theta / 2) * Y

def Rz(theta):
    """Returns the matrix representation of a Z-axis rotation."""
    return Matrix([[exp(-I * theta / 2), 0], [0, exp(I * theta / 2)]])

# Defined exactly as in the article's equation (Square root of X)
SX = (1/2) * Matrix([[1+I, 1-I], [1-I, 1+I]])

def Rzz(theta):
    """
    Returns the matrix representation of the parametric Rzz entanglement gate.
    Used in some IBM superconducting architectures.
    """
    return Matrix([
        [exp(-I * theta / 2), 0, 0, 0],
        [0, exp(I * theta / 2), 0, 0],
        [0, 0, exp(I * theta / 2), 0],
        [0, 0, 0, exp(-I * theta / 2)]
    ])

def GPi(phi):
    """
    Returns the matrix representation of the IonQ native GPi gate.
    Corresponds to a full Rabi oscillation (pi rotation) with laser phase phi.
    """
    return Matrix([[0, exp(-I * phi)], [exp(I * phi), 0]])

def GPi2(phi):
    """
    Returns the matrix representation of the IonQ native GPi2 gate.
    Corresponds to a partial Rabi oscillation (pi/2 rotation) with laser phase phi.
    """
    return (1/sp.sqrt(2)) * Matrix([[1, -I * exp(-I * phi)], [-I * exp(I * phi), 1]])

# Real symbolic variables for U3
theta, phi, lam = sp.symbols("theta phi lambda", real=True)

# ==========================================
# 2. Ideal Matrices (Direct and Flipped)
# ==========================================
U3_IDEAL = Matrix(
    [
        [cos(theta / 2), -exp(I * lam) * sin(theta / 2)],
        [exp(I * phi) * sin(theta / 2), exp(I * (phi + lam)) * cos(theta / 2)],
    ]
)

CNOT_01 = Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])

CNOT_10 = Matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]])

SWAP_IDEAL = Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])

# ==========================================
# 3. Symbolic Validator
# ==========================================
def check_equivalence(mat_A, ideal_name, ideal_mat_01, ideal_mat_10=None):
    """
    Validates if a decomposed matrix (mat_A) is algebraically equivalent to 
    an ideal matrix, disregarding global phase differences.

    Args:
        mat_A (sympy.Matrix): The computed matrix from the decomposition.
        ideal_name (str): Label for the validation output.
        ideal_mat_01 (sympy.Matrix): The target ideal matrix.
        ideal_mat_10 (sympy.Matrix, optional): The target ideal matrix with 
                                               reversed control/target topology.

    Returns:
        sympy.Matrix: The original mat_A, allowing for function chaining.
    """
    def verify(mat_ideal):
        # 1. Compute M = U_decomp * U_ideal^dagger
        # If they are equivalent, M should be the Identity matrix (times a global phase)
        M = mat_A * mat_ideal.adjoint()
        M = M.applyfunc(sp.simplify)

        if M.free_symbols:
            M = M.applyfunc(lambda x: sp.simplify(x.rewrite(sp.exp)))

        # 2. Check if the matrix is strictly diagonal (off-diagonal elements == 0)
        for i in range(M.rows):
            for j in range(M.cols):
                if i != j:
                    if not M[i, j].equals(0):
                        return False

        # 3. Check if the diagonal is homogeneous (ignores global phase)
        diag_0 = M[0, 0]
        if diag_0.equals(0):
            return False

        for i in range(1, M.rows):
            if not (M[i, i] - diag_0).equals(0):
                return False

        return True

    if verify(ideal_mat_01):
        print(f"[PASS] {ideal_name} (Alignment: CNOT_01)")
        return mat_A

    if ideal_mat_10 is not None and verify(ideal_mat_10):
        print(f"[PASS] {ideal_name} (Inverted Topology: CNOT_10)")
        return mat_A

    print(f"[FAIL] {ideal_name}\n       Matrices are not algebraically equivalent.")
    return mat_A

# ==========================================
# 4. IBM
# ==========================================
print("\n" + "=" * 45 + "\n IBM ARCHITECTURE\n" + "=" * 45)

u3_ibm = Rz(phi) * SX * Rz(pi - theta) * SX * Rz(lam - pi)
check_equivalence(u3_ibm, "IBM U3   (2x SX)", U3_IDEAL)

# IBM CNOT 1: Decomposition via CZ
CZ = Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
had_equiv = Rz(pi / 2) * SX * Rz(pi / 2)
cnot_ibm_cz = kron(I_mat, had_equiv) * CZ * kron(I_mat, had_equiv)
cnot_ibm_cz = check_equivalence(cnot_ibm_cz, "IBM CNOT (1x CZ)", CNOT_01, CNOT_10)

# IBM CNOT 2: Decomposition via Rzz(pi/2)
cnot_ibm_rzz = kron(Rz(-pi / 2), X * SX) * kron(I_mat, had_equiv) * Rzz(pi / 2) * kron(I_mat, had_equiv)
cnot_ibm_rzz = check_equivalence(cnot_ibm_rzz, "IBM CNOT (1x Rzz)", CNOT_01, CNOT_10)

# IBM CNOT 3: Decomposition via ECR
ECR = (eye(4) - I * kron(X, Y)) / sp.sqrt(2)
U_align = kron(Ry(pi / 2), Rz(pi / 2))
U_align_dag = kron(Ry(-pi / 2), Rz(-pi / 2))
flip = kron(X, I_mat)

core_ibm_ecr = flip * U_align_dag * ECR * U_align * flip
cnot_ibm_ecr = kron(Rz(pi / 2), Rx(pi / 2)) * core_ibm_ecr
cnot_ibm_ecr = check_equivalence(cnot_ibm_ecr, "IBM CNOT (1x ECR)", CNOT_01, CNOT_10)

# swap_ibm = cnot_ibm_ecr * (kron(H, H) * cnot_ibm_ecr * kron(H, H)) * cnot_ibm_ecr
# check_equivalence(swap_ibm, "IBM SWAP (3x ECR CNOTs)", SWAP_IDEAL)

# ==========================================
# 5. RIGETTI
# ==========================================
print("\n" + "=" * 45 + "\n RIGETTI ARCHITECTURE\n" + "=" * 45)

iSWAP = Matrix([[1, 0, 0, 0], [0, 0, I, 0], [0, I, 0, 0], [0, 0, 0, 1]])

u3_rigetti = Rz(phi) * Rx(-pi / 2) * Rz(theta) * Rx(pi / 2) * Rz(lam)
check_equivalence(u3_rigetti, "Rigetti U3   (2x Rx)", U3_IDEAL)

# Rigetti CNOT
pre_rigetti = kron(Rz(-pi / 2), Rz(pi / 2))
mid_rigetti = kron(Rx(pi / 2), I_mat)
post_rigetti = kron(I_mat, H * Rx(-pi / 2))

cnot_rigetti = post_rigetti * iSWAP * mid_rigetti * iSWAP * pre_rigetti
cnot_rigetti = check_equivalence(
    cnot_rigetti, "Rigetti CNOT (2x iSWAP)", CNOT_01, CNOT_10
)

# swap_rigetti = cnot_rigetti * (kron(H, H) * cnot_rigetti * kron(H, H)) * cnot_rigetti
# check_equivalence(swap_rigetti, "Rigetti SWAP (3x iSWAP CNOTs)", SWAP_IDEAL)

# ==========================================
# 6. IONQ
# ==========================================
print("\n" + "=" * 45 + "\n IONQ ARCHITECTURE\n" + "=" * 45)

MS = (eye(4) - I * kron(X, X)) / sp.sqrt(2)

u3_ionq = Rz(phi) * GPi2(pi) * Rz(theta) * GPi2(0) * Rz(lam)
check_equivalence(u3_ionq, "IonQ U3   (2x GPi2)", U3_IDEAL)

# IonQ CNOT
core_ionq = kron(GPi(0) * H, H) * MS * kron(H * GPi(0), H)
cz_ionq = kron(Rz(pi / 2), Rz(pi / 2)) * core_ionq
cnot_ionq = kron(I_mat, H) * cz_ionq * kron(I_mat, H)

cnot_ionq = check_equivalence(cnot_ionq, "IonQ CNOT (1x MS)", CNOT_01, CNOT_10)

# swap_ionq = cnot_ionq * (kron(H, H) * cnot_ionq * kron(H, H)) * cnot_ionq
# check_equivalence(swap_ionq, "IonQ SWAP (3x MS CNOTs)", SWAP_IDEAL)