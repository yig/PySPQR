"""Tests for complex matrix support in sparseqr."""

import numpy as np
import pytest
import scipy.sparse

import sparseqr


def adj(A):
    """Return the conjugate transpose (adjoint) of a matrix."""
    return A.conj().T


class TestComplexQR:
    """Tests for QR decomposition of complex matrices."""

    ATOL = 1e-14  # close to machine precision

    def test_complex_qr_unitary_q(self):
        """Test that Q is unitary for complex matrix (Q @ Q^H = I)."""
        A = scipy.sparse.random(m=120, n=100, density=0.1, dtype=np.complex128, random_state=42)
        Q, R, P, rank = sparseqr.qr(A)

        QQH = (Q @ adj(Q)).toarray()
        np.testing.assert_allclose(QQH, np.eye(Q.shape[0]), atol=self.ATOL)

    def test_complex_qr_upper_triangular_r(self):
        """Test that R is upper triangular for complex matrix."""
        A = scipy.sparse.random(m=120, n=100, density=0.1, dtype=np.complex128, random_state=42)
        Q, R, P, rank = sparseqr.qr(A)

        lower_triangle = np.tril(R.toarray(), k=-1)
        np.testing.assert_allclose(lower_triangle, 0, atol=self.ATOL)

    def test_complex_qr_reconstruction(self):
        """Test that A can be reconstructed from QR decomposition: A = (QR) @ P^T."""
        A = scipy.sparse.random(m=120, n=100, density=0.1, dtype=np.complex128, random_state=42)
        Q, R, P, rank = sparseqr.qr(A)

        P_matrix = sparseqr.permutation_vector_to_matrix(P)
        A_reconstructed = (Q @ R) @ P_matrix.T
        np.testing.assert_allclose(A_reconstructed.toarray(), A.toarray(), atol=self.ATOL)

    @pytest.mark.parametrize("shape", [(50, 50), (100, 50), (50, 100)])
    def test_complex_qr_various_shapes(self, shape):
        """Test complex QR on square, tall, and wide matrices."""
        m, n = shape
        A = scipy.sparse.random(m=m, n=n, density=0.1, dtype=np.complex128, random_state=42)
        Q, R, P, rank = sparseqr.qr(A)

        P_matrix = sparseqr.permutation_vector_to_matrix(P)
        A_reconstructed = (Q @ R) @ P_matrix.T
        np.testing.assert_allclose(A_reconstructed.toarray(), A.toarray(), atol=self.ATOL)


class TestComplexSolve:
    """Tests for solve function with complex matrices."""

    ATOL = 1e-10

    def test_complex_solve_exact_system(self):
        """Test solving an exact complex system where Ax=b has exact solution."""
        A = scipy.sparse.random(m=50, n=50, density=0.2, dtype=np.complex128, random_state=42)
        # Add diagonal to make it well-conditioned
        A = A + 5 * scipy.sparse.eye(50, dtype=np.complex128)
        x_true = np.random.RandomState(42).random(50) + 1j * np.random.RandomState(43).random(50)
        b = A @ x_true

        x = sparseqr.solve(A, b, tolerance=0)

        assert x is not None, "solve() returned None"
        np.testing.assert_allclose(x, x_true, atol=self.ATOL)

    def test_complex_solve_overdetermined_least_squares(self):
        """Test least-squares solution for overdetermined complex system."""
        # 100 equations, 50 unknowns
        A = scipy.sparse.random(m=100, n=50, density=0.2, dtype=np.complex128, random_state=42)
        x_true = np.random.RandomState(42).random(50) + 1j * np.random.RandomState(43).random(50)
        b = A @ x_true

        x = sparseqr.solve(A, b, tolerance=0)

        assert x is not None, "solve() returned None"
        assert x.shape == (50,), f"Solution shape wrong: {x.shape}"
        # For an exact system (b in range of A), solution should match
        np.testing.assert_allclose(x, x_true, atol=self.ATOL)

    def test_complex_solve_multiple_rhs(self):
        """Test solving AX=B with multiple complex RHS vectors."""
        A = scipy.sparse.random(m=50, n=50, density=0.2, dtype=np.complex128, random_state=42)
        A = A + 5 * scipy.sparse.eye(50, dtype=np.complex128)
        X_true = (np.random.RandomState(42).random((50, 3)) +
                  1j * np.random.RandomState(43).random((50, 3)))
        B = A @ X_true

        X = sparseqr.solve(A, B, tolerance=0)

        assert X is not None, "solve() returned None"
        assert X.shape == (50, 3), f"Solution shape wrong: {X.shape}"
        np.testing.assert_allclose(X, X_true, atol=self.ATOL)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
