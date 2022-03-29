import numpy as np

"""
scaling:
w[i]  : N! / [ (N-i)!   (i)!   ] p^i     (1-p)^(N-i)
w[i-1]: N! / [ (N-i+1)! (i-1)! ] p^(i-1) (1-p)^(N-i+1)
normalisation:
r[i]  : (-1)^i     [ (1-p)/p ]^i     i!     / (-N)_i 
r[i-1]: (-1)^(i-1) [ (1-p)/p ]^(i-1) (i-1)! / (-N)_(i-1) 
"""
def norm_factors(N, p):
    p_odd = p / (1 - p)
    w = np.zeros(N)
    rho_inv = np.zeros(N)
    w[0] = (1 - p) ** (N - 1)
    rho_inv[0] = 1
    for i in range(1, N):
        w[i] = w[i - 1] * (N - i) / i * p_odd
        rho_inv[i] = 1 / (-1 / p_odd * i / (-N + i)) * rho_inv[i - 1]
    return w, rho_inv


def krawtchouk_poly(N, p):
    K = np.zeros((N, N))
    x = np.array(range(N))
    K[0, :] = 1
    K[1, :] = (1 - x / (p * (N - 1)))
    for i in range(1, N - 1):
        K[i + 1] = \
            ((N - 1) * p - 2 * i * p + i - x) * K[i, :] \
                - (i * (1 - p) * K[i - 1, :])
        K[i + 1] /= (p * (N - i - 1))
    return K


def krawtchouk_poly_weighted(N, p):
    w, rho_inv = norm_factors(N, p)
    K = krawtchouk_poly(N, p)
    K_bar = K * np.outer(np.sqrt(rho_inv), np.sqrt(w))
    return K_bar


def krawtchouk_moments(img, N, M, p):
    K_bar1 = K_bar2 = krawtchouk_poly_weighted(N, p[0])
    if (N != M) or (p[0] != p[1]):
        K_bar2 = krawtchouk_poly_weighted(M, p[1])
    Q = (K_bar2 @ (K_bar1 @ img).T).T
    return Q