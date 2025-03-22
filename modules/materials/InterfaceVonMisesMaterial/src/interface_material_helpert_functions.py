import numpy as np

def interface_geometry_system_couplings(
    I: np.ndarray,
    J: np.ndarray,
    N: np.ndarray,
    T: np.ndarray,
    L: np.ndarray,
    M: np.ndarray,
):

    newshape = np.append(L.shape, N.shape[2:])
    L = np.broadcast_to(L[:, :, :, :, np.newaxis], newshape)
    M = np.broadcast_to(M[:, :, :, :, np.newaxis], newshape)


    Q = np.einsum("aibjq,ijq->abq", L, N)

    G = np.zeros(Q.shape)

    for q in range(G.shape[-1]):
        G[:, :, q] = np.linalg.inv(Q[:, :, q])

    A = np.einsum("abq,ijq->aibjq", G, N)
    LA = np.einsum("aimnq,mnbjq->aibjq", L, A)
    LAL = np.einsum("aimnq,mnbjq->aibjq", LA, L)

    B = L - LAL

    # check that the expressions are consistent
    # assert_consistent_arrays_indices(L, A, B, M, J, I)

    return M, B, L, A, G
    
def calculate_material_matrices(
    n, I, J, N, T, C_0_aibj, C_M_aibj, C_I_aibj, S_0_aibj, S_M_aibj, S_I_aibj
):

    M_0, B_0, L_0, A_0, G_0 = interface_geometry_system_couplings(
        I, J, N, T, C_0_aibj, S_0_aibj
    )

    M_M, B_M, L_M, A_M, G_M = interface_geometry_system_couplings(
        I, J, N, T, C_M_aibj, S_M_aibj
    )

    M_I, B_I, L_I, A_I, G_I = interface_geometry_system_couplings(
        I, J, N, T, C_I_aibj, S_I_aibj
    )

    F = (
        2.0 * np.einsum("aimnq,mnbjq->aibjq", M_0, B_0)
        - np.einsum("aimnq,mnbjq->aibjq", M_M, B_M)
        - np.einsum("aimnq,mnbjq->aibjq", M_I, B_I)
    )

    Y = (
        np.einsum("aimnq,mnbjq->aibjq", L_M, A_M)
        + np.einsum("aimnq,mnbjq->aibjq", L_I, A_I)
        - 2.0 * np.einsum("aimnq,mnbjq->aibjq", L_0, A_0)
    )

    # assert_equivalent_F_Falt_Y(F, F_alt, Y)

    H = 2.0 * G_0 - G_M - G_I
    Z = B_M + B_I - 2.0 * B_0

    H_inv = np.zeros(H.shape)
    for q in range(H.shape[-1]):
        H_inv[:, :, q] = np.linalg.inv(H[:, :, q])

    nF = np.einsum("aq,aibjq->ibjq", n, F)
    Fn = np.einsum("aibjq,jq->aibq", F, n)
    Yn = np.einsum("aibjq,jq->aibq", Y, n)
    H_inv_nF = np.einsum("abq,bijq->aijq", H_inv, nF)
    Yn_H_inv_Fn = np.einsum("aimq,mnq,nbjq->aibjq", Yn, H_inv, Fn)

    return Z, H_inv, H_inv_nF, Yn_H_inv_Fn


# Utility routines compute Voigt to tensor


def voigt_to_tensor(C_voigt):
    """
    Convert 6x6 Voigt matrix to fourth-order tensor in 3D.
    The material is always treated as 3D
    """

    index_map = {0: (0, 0), 1: (1, 1), 2: (2, 2), 3: (1, 2), 4: (0, 2), 5: (0, 1)}

    C_tensor = np.zeros((3, 3, 3, 3))

    for I in range(6):
        for J in range(6):
            i, j = index_map[I]
            k, l = index_map[J]
            factor = 0.5 if I >= 3 else 1.0
            factor *= 0.5 if J >= 3 else 1.0
            C_tensor[i, j, k, l] = C_voigt[I, J] * factor
            C_tensor[j, i, k, l] = C_voigt[I, J] * factor
            C_tensor[i, j, l, k] = C_voigt[I, J] * factor
            C_tensor[j, i, l, k] = C_voigt[I, J] * factor

    return C_tensor
    
def assert_consistent_arrays_inices(
    L: np.ndarray, A: np.ndarray, B: np.ndarray, M: np.ndarray, J: np.ndarray, I: np.ndarray):

    #P_norm = np.einsum("abq,ijq->aibjq", J, N)  # aibj
    #P_tan = np.einsum("abq,ijq->aibjq", J, T)  # aibj

    JI = np.einsum("abq,ijq->aibjq", J, I)
    
    LM = np.einsum("aimn,mnbj->aibj", L, M)
    ML = np.einsum("aimn,mnbj->aibj", M, L)

    assert np.allclose(
        LM, JI, atol=10.0 ** (-8)
    ), "Arrays are not equal within the tolerance."

    LA = np.einsum("aimn,mnbj->aibj", L, A)
    BM = np.einsum("aimn,mnbj->aibj", B, M)

    assert np.allclose(
        LA + BM, JI, atol=10.0 ** (-8)
    ), "Arrays are not equal within the tolerance."

    return "Passed the assertion about proper definition of matrices A, B"


def assert_equivalent_F_Falt_Y(
    F: np.ndarray,
    F_alt: np.ndarray,
    Y: np.ndarray,
    A_0: np.ndarray,
    L_0: np.ndarray,
    A_M: np.ndarray,
    L_M: np.ndarray,
    A_I: np.ndarray,
    L_I: np.ndarray,
):

    F_alt = (
        -2.0 * np.einsum("aimn,mnbj->aibj", A_0, L_0)
        + np.einsum("aimn,mnbj->aibj", A_M, L_M)
        + np.einsum("aimn,mnbj->aibj", A_I, L_I)
    )

    assert np.allclose(
        F, F_alt, atol=10.0 ** (-8)
    ), "Arrays are not equal within the tolerance."
    assert np.allclose(
        Y.T, F, atol=10.0 ** (-8)
    ), "Arrays are not equal within the tolerance."
    return
