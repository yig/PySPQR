import numpy as np
from scipy.sparse import random as sparse_random
from scipy.linalg import qr as scipy_qr
import sparseqr

def adj(A):
    """Return the complex conjugate of a matrix."""
    return A.conj().T

dim = (120,100)
A = sparse_random(m=dim[0], n=dim[1], density=0.1, dtype=np.complex128)
# call pyspqr
Q, R, P, _ = sparseqr.qr(A)

atol = 1e-14  # something close to machine precision
# assert Q is unitary
np.testing.assert_allclose((Q @ adj(Q)).toarray(), np.eye(dim[0]), atol=atol)
# assert R is right upper triangular -> values under diagonal must be zero
lower_triangle = np.tril(R.toarray(), k=-1)
np.testing.assert_allclose(lower_triangle, 0, atol=atol)
# Reconstruct A via AP = QR => A = (QR) P^T
P_matrix = sparseqr.permutation_vector_to_matrix(P)
A_reconstructed = (Q @ R) @ P_matrix.T
np.testing.assert_allclose(A_reconstructed.toarray(), A.toarray(), atol=atol)
