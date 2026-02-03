"""Tests for complex matrix support in sparseqr."""

import numpy as np
import pytest
from scipy.sparse import random as sparse_random

import sparseqr


def adj(A):
    """Return the conjugate transpose (adjoint) of a matrix."""
    return A.conj().T


class TestComplexQR:
    """Tests for QR decomposition of complex matrices."""

    ATOL = 1e-14  # close to machine precision

    def test_complex_qr_unitary_q(self):
        """Test that Q is unitary for complex matrix (Q @ Q^H = I)."""
        A = sparse_random(m=120, n=100, density=0.1, dtype=np.complex128, random_state=42)
        Q, R, P, rank = sparseqr.qr(A)

        QQH = (Q @ adj(Q)).toarray()
        np.testing.assert_allclose(QQH, np.eye(Q.shape[0]), atol=self.ATOL)

    def test_complex_qr_upper_triangular_r(self):
        """Test that R is upper triangular for complex matrix."""
        A = sparse_random(m=120, n=100, density=0.1, dtype=np.complex128, random_state=42)
        Q, R, P, rank = sparseqr.qr(A)

        lower_triangle = np.tril(R.toarray(), k=-1)
        np.testing.assert_allclose(lower_triangle, 0, atol=self.ATOL)

    def test_complex_qr_reconstruction(self):
        """Test that A can be reconstructed from QR decomposition: A = (QR) @ P^T."""
        A = sparse_random(m=120, n=100, density=0.1, dtype=np.complex128, random_state=42)
        Q, R, P, rank = sparseqr.qr(A)

        P_matrix = sparseqr.permutation_vector_to_matrix(P)
        A_reconstructed = (Q @ R) @ P_matrix.T
        np.testing.assert_allclose(A_reconstructed.toarray(), A.toarray(), atol=self.ATOL)

    @pytest.mark.parametrize("shape", [(50, 50), (100, 50), (50, 100)])
    def test_complex_qr_various_shapes(self, shape):
        """Test complex QR on square, tall, and wide matrices."""
        m, n = shape
        A = sparse_random(m=m, n=n, density=0.1, dtype=np.complex128, random_state=42)
        Q, R, P, rank = sparseqr.qr(A)

        P_matrix = sparseqr.permutation_vector_to_matrix(P)
        A_reconstructed = (Q @ R) @ P_matrix.T
        np.testing.assert_allclose(A_reconstructed.toarray(), A.toarray(), atol=self.ATOL)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
