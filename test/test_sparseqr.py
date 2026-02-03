"""Tests for sparseqr module.

Run with: pytest test/ -v
Or simply: python -m pytest
"""

import numpy as np
import scipy.sparse
import pytest

import sparseqr


class TestQR:
    """Tests for sparseqr.qr() function."""

    def test_qr_basic_decomposition(self):
        """Test that Q*R = M*E for a random sparse matrix."""
        M = scipy.sparse.rand(10, 10, density=0.1, random_state=42)
        Q, R, E, rank = sparseqr.qr(M)

        # Q*R should equal M*E
        E_matrix = sparseqr.permutation_vector_to_matrix(E)
        residual = abs(Q @ R - M @ E_matrix).sum()
        assert residual < 1e-10, f"QR decomposition residual too large: {residual}"

    def test_qr_rectangular_tall(self):
        """Test QR on a tall (overdetermined) matrix."""
        M = scipy.sparse.rand(20, 10, density=0.1, random_state=42)
        Q, R, E, rank = sparseqr.qr(M)

        assert Q.shape == (20, 20), f"Q shape wrong: {Q.shape}"
        assert R.shape == (20, 10), f"R shape wrong: {R.shape}"

        E_matrix = sparseqr.permutation_vector_to_matrix(E)
        residual = abs(Q @ R - M @ E_matrix).sum()
        assert residual < 1e-10

    def test_qr_rectangular_wide(self):
        """Test QR on a wide (underdetermined) matrix."""
        M = scipy.sparse.rand(10, 20, density=0.1, random_state=42)
        Q, R, E, rank = sparseqr.qr(M)

        assert Q.shape == (10, 10), f"Q shape wrong: {Q.shape}"
        assert R.shape == (10, 20), f"R shape wrong: {R.shape}"

    def test_qr_economy_mode(self):
        """Test economy=True produces smaller Q and R."""
        M = scipy.sparse.rand(20, 5, density=0.1, random_state=42)

        # Full QR
        Q_full, R_full, E, rank = sparseqr.qr(M, economy=False)
        assert Q_full.shape == (20, 20)
        assert R_full.shape == (20, 5)

        # Economy QR
        Q_econ, R_econ, E, rank = sparseqr.qr(M, economy=True)
        assert Q_econ.shape == (20, 5), f"Economy Q shape wrong: {Q_econ.shape}"
        assert R_econ.shape == (5, 5), f"Economy R shape wrong: {R_econ.shape}"

    def test_qr_identity_matrix(self):
        """Test QR decomposition of identity matrix."""
        I = scipy.sparse.eye(5, format='coo')
        Q, R, E, rank = sparseqr.qr(I)

        assert rank == 5, f"Rank should be 5, got {rank}"

        # Q and R should both be identity (up to permutation)
        E_matrix = sparseqr.permutation_vector_to_matrix(E)
        residual = abs(Q @ R - I @ E_matrix).sum()
        assert residual < 1e-10

    def test_qr_with_tolerance(self):
        """Test QR with explicit tolerance parameter."""
        M = scipy.sparse.rand(10, 10, density=0.2, random_state=42)
        Q, R, E, rank = sparseqr.qr(M, tolerance=0)

        E_matrix = sparseqr.permutation_vector_to_matrix(E)
        residual = abs(Q @ R - M @ E_matrix).sum()
        assert residual < 1e-10


class TestSolve:
    """Tests for sparseqr.solve() function."""

    def test_solve_dense_rhs_single(self):
        """Test solving Ax=b with single dense RHS vector."""
        # Create a well-conditioned matrix
        A = scipy.sparse.diags([1, 2, 3, 4, 5], format='coo', dtype=np.float64) + scipy.sparse.rand(5, 5, density=0.1, random_state=42)
        b = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

        x = sparseqr.solve(A, b, tolerance=0)

        assert x is not None, "solve() returned None"
        assert x.shape == (5,), f"Solution shape wrong: {x.shape}"

    def test_solve_dense_rhs_multiple(self):
        """Test solving AX=B with multiple dense RHS vectors."""
        A = scipy.sparse.diags([1, 2, 3, 4, 5], format='coo', dtype=np.float64) + scipy.sparse.rand(5, 5, density=0.1, random_state=42)
        B = np.random.RandomState(42).random((5, 3))

        X = sparseqr.solve(A, B, tolerance=0)

        assert X is not None, "solve() returned None"
        assert X.shape == (5, 3), f"Solution shape wrong: {X.shape}"

    def test_solve_sparse_rhs(self):
        """Test solving AX=B with sparse RHS."""
        A = scipy.sparse.diags([1, 2, 3, 4, 5], format='coo', dtype=np.float64) + scipy.sparse.rand(5, 5, density=0.1, random_state=42)
        B = scipy.sparse.rand(5, 3, density=0.3, random_state=42)

        X = sparseqr.solve(A, B, tolerance=0)

        assert X is not None, "solve() returned None"
        assert scipy.sparse.issparse(X), "Solution should be sparse for sparse RHS"
        assert X.shape == (5, 3), f"Solution shape wrong: {X.shape}"

    def test_solve_overdetermined_least_squares(self):
        """Test least-squares solution for overdetermined system."""
        # 20 equations, 10 unknowns
        A = scipy.sparse.rand(20, 10, density=0.2, random_state=42)
        b = np.random.RandomState(42).random(20)

        x = sparseqr.solve(A, b, tolerance=0)

        assert x is not None, "solve() returned None"
        assert x.shape == (10,), f"Solution shape wrong: {x.shape}"

    def test_solve_exact_system(self):
        """Test solving an exact system where Ax=b has exact solution."""
        # Create a system with known solution
        A = scipy.sparse.eye(5, format='coo')
        x_true = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        b = A @ x_true

        x = sparseqr.solve(A, b, tolerance=0)

        assert x is not None
        np.testing.assert_allclose(x, x_true, rtol=1e-10)


class TestRZ:
    """Tests for sparseqr.rz() function."""

    def test_rz_basic(self):
        """Test rz() returns expected outputs."""
        A = scipy.sparse.rand(20, 10, density=0.1, random_state=42)
        b = np.random.RandomState(42).random(20)

        Z, R, E, rank = sparseqr.rz(A, b, tolerance=0)

        assert Z is not None
        assert R is not None
        assert isinstance(rank, (int, np.integer))

    def test_rz_output_shapes(self):
        """Test rz() returns correct shapes."""
        m, n = 15, 8
        A = scipy.sparse.rand(m, n, density=0.2, random_state=42)
        b = np.random.RandomState(42).random(m)

        Z, R, E, rank = sparseqr.rz(A, b, tolerance=0)

        # Z should have n rows (one per unknown)
        assert Z.shape[0] == n, f"Z rows should be {n}, got {Z.shape[0]}"


class TestQRFactorize:
    """Tests for sparseqr.qr_factorize() and qmult() functions."""

    def test_qr_factorize_basic(self):
        """Test qr_factorize() returns a factorization object."""
        M = scipy.sparse.rand(100, 100, density=0.05, random_state=42)
        QR = sparseqr.qr_factorize(M)

        assert QR is not None

    def test_qmult_basic(self):
        """Test qmult() with Householder form QR."""
        M = scipy.sparse.rand(100, 100, density=0.05, random_state=42)
        QR = sparseqr.qr_factorize(M)

        X = np.zeros((M.shape[0], 1))
        X[-1, 0] = 1

        Y = sparseqr.qmult(QR, X)

        assert Y.shape == (100, 1), f"Y shape wrong: {Y.shape}"

    def test_qmult_methods(self):
        """Test different qmult methods (Q'X, QX, XQ', XQ)."""
        M = scipy.sparse.rand(20, 20, density=0.2, random_state=42)
        QR = sparseqr.qr_factorize(M)

        X = np.random.RandomState(42).random((20, 1))

        # method=0: Q'X
        Y0 = sparseqr.qmult(QR, X, method=0)
        assert Y0.shape == (20, 1)

        # method=1: QX (default)
        Y1 = sparseqr.qmult(QR, X, method=1)
        assert Y1.shape == (20, 1)


class TestPermutationVectorToMatrix:
    """Tests for sparseqr.permutation_vector_to_matrix() function."""

    def test_identity_permutation(self):
        """Test identity permutation [0,1,2,3,4]."""
        E = np.array([0, 1, 2, 3, 4])
        P = sparseqr.permutation_vector_to_matrix(E)

        expected = scipy.sparse.eye(5)
        diff = abs(P - expected).sum()
        assert diff < 1e-10, "Identity permutation should produce identity matrix"

    def test_swap_permutation(self):
        """Test a permutation that swaps first and last."""
        E = np.array([4, 1, 2, 3, 0])
        P = sparseqr.permutation_vector_to_matrix(E)

        # P[E[k], k] = 1, so P[4,0]=1, P[1,1]=1, P[2,2]=1, P[3,3]=1, P[0,4]=1
        assert P.shape == (5, 5)
        assert P.nnz == 5  # exactly 5 non-zeros

        # Convert to dense for easier checking
        P_dense = P.toarray()
        assert P_dense[4, 0] == 1
        assert P_dense[0, 4] == 1

    def test_permutation_orthogonality(self):
        """Test that permutation matrix is orthogonal (P @ P.T = I)."""
        E = np.array([2, 0, 4, 1, 3])
        P = sparseqr.permutation_vector_to_matrix(E)

        result = P @ P.T
        expected = scipy.sparse.eye(5)
        diff = abs(result - expected).sum()
        assert diff < 1e-10


class TestEdgeCases:
    """Tests for edge cases and special matrices."""

    def test_very_sparse_matrix(self):
        """Test with very sparse matrix."""
        M = scipy.sparse.rand(50, 50, density=0.01, random_state=42)
        Q, R, E, rank = sparseqr.qr(M)

        E_matrix = sparseqr.permutation_vector_to_matrix(E)
        residual = abs(Q @ R - M @ E_matrix).sum()
        assert residual < 1e-10

    def test_diagonal_matrix(self):
        """Test with diagonal matrix."""
        M = scipy.sparse.diags([1, 2, 3, 4, 5], format='coo', dtype=np.float64)
        Q, R, E, rank = sparseqr.qr(M)

        assert rank == 5
        E_matrix = sparseqr.permutation_vector_to_matrix(E)
        residual = abs(Q @ R - M @ E_matrix).sum()
        assert residual < 1e-10

    def test_single_element_matrix(self):
        """Test with 1x1 matrix."""
        M = scipy.sparse.coo_matrix([[3.0]])
        Q, R, E, rank = sparseqr.qr(M)

        assert rank == 1
        assert Q.shape == (1, 1)
        assert R.shape == (1, 1)

    def test_different_sparse_formats(self):
        """Test that different input formats work (coo, csr, csc)."""
        base = scipy.sparse.rand(10, 10, density=0.2, random_state=42)

        for fmt in ['coo', 'csr', 'csc', 'lil']:
            M = base.asformat(fmt)
            Q, R, E, rank = sparseqr.qr(M)
            E_matrix = sparseqr.permutation_vector_to_matrix(E)
            residual = abs(Q @ R - M @ E_matrix).sum()
            assert residual < 1e-10, f"Failed for format {fmt}"


class TestNumericalAccuracy:
    """Tests for numerical accuracy."""

    def test_orthogonality_of_q(self):
        """Test that Q is orthogonal (Q.T @ Q ≈ I for the relevant columns)."""
        M = scipy.sparse.rand(10, 10, density=0.3, random_state=42)
        Q, R, E, rank = sparseqr.qr(M)

        # For a square matrix with full rank, Q should be orthogonal
        if rank == 10:
            QtQ = (Q.T @ Q).toarray()
            I = np.eye(10)
            np.testing.assert_allclose(QtQ, I, atol=1e-10)

    def test_upper_triangular_r(self):
        """Test that R is upper triangular."""
        M = scipy.sparse.rand(10, 10, density=0.3, random_state=42)
        Q, R, E, rank = sparseqr.qr(M)

        R_dense = R.toarray()
        # Check lower triangle is zero (below diagonal)
        for i in range(R_dense.shape[0]):
            for j in range(min(i, R_dense.shape[1])):
                assert abs(R_dense[i, j]) < 1e-10, f"R[{i},{j}] = {R_dense[i,j]} should be zero"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
