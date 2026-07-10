# Quant-Finance-Summer-2026

# Testing Roughness in Realized Volatility

This project studies whether the roughness estimated from realized volatility is a reliable indicator of the roughness of the underlying spot volatility, or whether it is affected by how realized volatility is constructed.

The main question is:

> Is the roughness observed in realized volatility a genuine market feature, or is it introduced by the estimation procedure?

## Project Overview

Volatility is central in quantitative finance because it affects option pricing, hedging, and portfolio risk. Empirical volatility paths often look highly irregular, and this irregularity is commonly summarized by a Hurst parameter $H$. Brownian-like roughness corresponds to $H \approx 0.5$, while rough-volatility models often use smaller values such as $H \approx 0.1$ to $0.3$.

However, spot volatility is not directly observed. In practice, we observe stock prices and construct realized volatility from log returns. This project tests whether the roughness of realized volatility reflects the true roughness of spot volatility.

## Main Idea

We distinguish between spot volatility and realized volatility.

### Spot Volatility

A standard log-price model is

$$
dX_t = \mu_t dt + \sigma_t dB_t,
\qquad X_t = \log S_t.
$$

The process $\sigma_t$ is the spot volatility. It is the instantaneous volatility driving the price path, but it is latent and cannot be directly observed from market prices.

### Realized Volatility

Given observed log returns

$$
\Delta X_i = X_{t_{i+1}} - X_{t_i},
$$

realized volatility over a window is


$$RV_{t,t+\Delta}
=
\left(
\sum_{t_i \in [t,t+\Delta]} (\Delta X_i)^2
\right)^{1/2}.
$$

Realized volatility is observable, but it is a constructed proxy. It is not the same object as spot volatility.

## Roughness Estimator

For observations

$$
X_0, X_1, \ldots, X_L
$$

on $[0,T]$, the path is split into $K$ coarse blocks. The normalized $p$-variation statistic is

$$
W(p)=\sum_{i=0}^{K-1}
\frac{|X_{b_{i+1}}-X_{b_i}|^p}
{\sum_{j=b_i}^{b_{i+1}-1}|X_{j+1}-X_j|^p}
(t_{b_{i+1}}-t_{b_i}).
$$

The estimated variation exponent is chosen by

$$
\widehat p
\approx
\mathrm{argmin}_{p\in[1,12]} |W(p)-T|.
$$

The estimated Hurst parameter is

$$
\widehat H = \frac{1}{\widehat p}.
$$

Interpretation:

- $H \approx 0.5$: Brownian-like roughness.
- $H < 0.5$: rougher than Brownian motion.
- Smaller $H$ means a more irregular path.

## Experiments

### Experiment 1: Fractional Brownian Motion Sanity Check

We first test the estimator on fractional Brownian motion paths, where the true Hurst parameter is known.

We simulate paths with

$$
H \in \{0.1, 0.3, 0.5, 0.8\}.
$$

The goal is to check whether the estimator can recover the known roughness.

Example result:

| True $H$ | Mean $\widehat H$ | Std. |
|---|---:|---:|
| 0.1 | 0.1060 | 0.0201 |
| 0.3 | 0.2991 | 0.0197 |
| 0.5 | 0.5036 | 0.0117 |
| 0.8 | 0.7848 | 0.0206 |

The estimator performs well on simulated fractional Brownian motion.

### Experiment 2: Fractional OU Volatility Model

We simulate a fractional Ornstein--Uhlenbeck stochastic-volatility model:

$$
dY_t = -\gamma Y_t dt + \theta dB_t^H,
\qquad
\sigma_t = \sigma_0 e^{Y_t}.
$$

Then we simulate the log-price path:

$$
dX_t = -\frac{1}{2}\sigma_t^2 dt + \sigma_t dW_t,
\qquad S_t=e^{X_t}.
$$

From the simulated price path, we compute realized volatility:

$$
RV_i =
\left(
\sum_{j=i-m+1}^{i}(\Delta X_j)^2
\right)^{1/2}.
$$

We then compare $\widehat H(\sigma)$ and $\widehat H(RV)$.

Example result using $N=90000$, $K=300$, and $m=10$:

| True $H$ | $\widehat H(\sigma)$ | $\widehat H(RV)$ |
|---|---:|---:|
| 0.10 | 0.088 | 0.224 |
| 0.20 | 0.211 | 0.218 |
| 0.30 | 0.293 | 0.201 |
| 0.40 | 0.392 | 0.145 |
| 0.50 | 0.499 | 0.166 |
| 0.60 | 0.597 | 0.164 |
| 0.70 | 0.683 | 0.163 |
| 0.80 | 0.773 | 0.158 |

The estimator recovers the roughness of the instantaneous volatility reasonably well:

$$
\widehat H(\sigma) \approx H.
$$

However, the realized-volatility roughness remains small across the tested values of $H$.

### Experiment 3: Sensitivity to Realized-Volatility Window Size

We fix the true volatility roughness at $H=0.3$ and change the realized-volatility window size $m$.

| RV window $m$ | Window length $\Delta$ | $\widehat H(\sigma)$ | $\widehat H(RV)$ |
|---:|---:|---:|---:|
| 10 | 0.00011 | 0.29056 | 0.18523 |
| 30 | 0.00033 | 0.28049 | 0.32981 |
| 100 | 0.00111 | 0.28381 | 0.49889 |
| 300 | 0.00333 | 0.30167 | 0.65755 |

The estimate from instantaneous volatility remains close to the true value $H=0.3$, while the estimate from realized volatility changes strongly with the window size.

This shows that realized-volatility roughness is sensitive to the construction of the realized-volatility path.

### Experiment 4: AAPL Market Data

We download 1-minute AAPL close prices using `yfinance`, keep only regular trading hours, and compute realized volatility using different window sizes.

Data summary:

- Ticker: AAPL
- Frequency: 1-minute prices
- Start: 2026-06-10, 09:30
- End: 2026-07-08, 15:59
- Trading days: 19
- Price observations: 7410

Estimated roughness of AAPL realized volatility:

| RV window | Length of RV path | $\widehat H(RV)$ |
|---:|---:|---:|
| 2 | 7372 | 0.1132 |
| 5 | 7315 | 0.2168 |
| 10 | 7220 | 0.3267 |
| 15 | 7125 | 0.4258 |
| 30 | 6840 | 0.5025 |
| 60 | 6270 | 0.6313 |

The estimated roughness increases as the realized-volatility window size increases. Small windows produce rougher realized-volatility paths, while larger windows smooth the path.

## Main Conclusion

Realized-volatility roughness is not a fixed property of the market path alone. It depends strongly on how realized volatility is constructed.

A low estimate of $\widehat H(RV)$ should not automatically be interpreted as proof that spot volatility is rough. The realized-volatility window size must be reported and tested carefully.

## Final Message

This project shows that realized volatility can look rough, but the estimated roughness depends strongly on the construction method. Therefore, roughness estimates from realized volatility should be interpreted with caution and reported together with robustness checks.
