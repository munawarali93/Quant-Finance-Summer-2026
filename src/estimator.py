import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yfinance as yf


def davies_harte(N, T, H):
    gamma = lambda k,H: 0.5*(np.abs(k-1)**(2*H) - 2*np.abs(k)**(2*H) + np.abs(k+1)**(2*H))
    g = [gamma(k,H) for k in range(0,N)];    r = g + [0] + g[::-1][0:N-1]

    # Step 1 (eigenvalues)
    j = np.arange(0,2*N);   k = 2*N-1
    lk = np.fft.fft(r*np.exp(2*np.pi*complex(0,1)*k*j*(1/(2*N))))[::-1]

    # Step 2 (get random variables)
    Vj = np.zeros((2*N,2), dtype=complex);
    Vj[0,0] = np.random.standard_normal();  Vj[N,0] = np.random.standard_normal()

    for i in range(1,N):
        Vj1 = np.random.standard_normal();    Vj2 = np.random.standard_normal()
        Vj[i][0] = Vj1; Vj[i][1] = Vj2; Vj[2*N-i][0] = Vj1;    Vj[2*N-i][1] = Vj2

    # Step 3 (compute Z)
    wk = np.zeros(2*N, dtype=complex)
    wk[0] = np.sqrt((lk[0]/(2*N)))*Vj[0][0];
    wk[1:N] = np.sqrt(lk[1:N]/(4*N))*((Vj[1:N].T[0]) + (complex(0,1)*Vj[1:N].T[1]))
    wk[N] = np.sqrt((lk[0]/(2*N)))*Vj[N][0]
    wk[N+1:2*N] = np.sqrt(lk[N+1:2*N]/(4*N))*(np.flip(Vj[1:N].T[0]) - (complex(0,1)*np.flip(Vj[1:N].T[1])))

    Z = np.fft.fft(wk);     fGn = Z[0:N]
    fBm = np.cumsum(fGn)*(N**(-H))
    fBm = (T**H)*(fBm)
    path = np.array([0] + list(fBm))
    return path

##Roughness Estimator


#1.   Making blocks
#2.   W_statistics
#  *   Computes the normalized p-variation statistic W(p).
#  *   Input
#        x : one-dimensional time series,
#        p : power,
#        K : number of coarse blocks,
#        T : time horizon,
#  *   Output: W(p)
#3. estimate_H: The estimator solves approximately W(p_hat) = T, and returns H_hat = 1 / p_hat.

def _make_blocks(L, K):
    K = int(K)
    K = max(2, min(K, L))
    blocks = np.linspace(0, L, K + 1, dtype=int)
    blocks = np.unique(blocks)
    if len(blocks) < 3:
        raise ValueError("Too few non-empty blocks. Try smaller K.")
    return blocks


def W_statistic(x, p, K=None, T=1.0, eps=1e-14):
    x = np.asarray(x, dtype=float).reshape(-1)
    x = x[np.isfinite(x)]

    if len(x) < 3:
        raise ValueError("Need at least 3 observations.")

    L = len(x) - 1

    if K is None:
        K = int(np.floor(np.sqrt(L)))

    blocks = _make_blocks(L, K)

    t = np.linspace(0.0, T, len(x))
    dx = np.diff(x)
    abs_dx = np.abs(dx)

    W = 0.0

    for start, end in zip(blocks[:-1], blocks[1:]):
        numerator = np.abs(x[end] - x[start]) ** p
        denominator = np.sum(abs_dx[start:end] ** p)

        if denominator <= eps or not np.isfinite(denominator):
            continue

        block_length = t[end] - t[start]
        W += (numerator / denominator) * block_length

    return float(W)


def estimate_H(
    x,
    K=None,
    T=1.0,
    p_min=1.05,
    p_max=12.0,
    n_grid=1000,
    p_grid=None,
    return_details=False
):

    x = np.asarray(x, dtype=float).reshape(-1)
    x = x[np.isfinite(x)]

    if len(x) < 3:
        raise ValueError("Need at least 3 finite observations.")

    dx = np.diff(x)
    scale = np.nanstd(dx)

    if np.isfinite(scale) and scale > 0:
        x_work = (x - np.nanmean(x)) / scale
    else:
        x_work = x - np.nanmean(x)

    # Candidate p-values
    if p_grid is None:
        p_grid = np.linspace(p_min, p_max, n_grid)
    else:
        p_grid = np.asarray(p_grid, dtype=float)

    # Compute W(p) for all p
    W_values = np.array([
        W_statistic(x_work, p, K=K, T=T)
        for p in p_grid
    ])

    valid = np.isfinite(W_values)

    if valid.sum() < 3:
        if return_details:
            details = {
                "H_hat": np.nan,
                "p_hat": np.nan,
                "K": int(K) if K is not None else int(np.floor(np.sqrt(len(x) - 1))),
                "T": T,
                "p_grid": p_grid,
                "W_values": W_values,
                "method": "failed_not_enough_valid_values",
                "crossings": [],
                "min_abs_W_minus_T": np.nan,
            }
            return np.nan, details

        return np.nan

    p_valid = p_grid[valid]
    W_valid = W_values[valid]

    y = W_valid - T

    # First: find closest grid point
    closest_idx = int(np.argmin(np.abs(y)))
    p_closest = float(p_valid[closest_idx])

    # Second: improve by linear interpolation if a crossing exists
    crossings = []

    for i in range(len(p_valid) - 1):
        y0 = y[i]
        y1 = y[i + 1]

        if y0 == 0:
            crossings.append(float(p_valid[i]))

        elif y0 * y1 < 0:
            p0, p1 = p_valid[i], p_valid[i + 1]
            w0, w1 = W_valid[i], W_valid[i + 1]

            if w1 != w0:
                p_cross = p0 + (T - w0) * (p1 - p0) / (w1 - w0)
                crossings.append(float(p_cross))

    if len(crossings) > 0:
        p_hat = min(crossings, key=lambda p: abs(p - p_closest))
        method = "linear_interpolated_crossing"
    else:
        p_hat = p_closest
        method = "closest_grid_point_no_crossing"

    H_hat = 1.0 / p_hat

    if return_details:
        details = {
            "H_hat": float(H_hat),
            "p_hat": float(p_hat),
            "K": int(K) if K is not None else int(np.floor(np.sqrt(len(x) - 1))),
            "T": T,
            "p_grid": p_grid,
            "W_values": W_values,
            "p_valid": p_valid,
            "W_valid": W_valid,
            "method": method,
            "crossings": crossings,
            "min_abs_W_minus_T": float(np.min(np.abs(y))),
        }

        return float(H_hat), details

    return float(H_hat)
