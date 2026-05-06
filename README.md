# svartools

Stata tools for simulating and visualizing structural VAR models.

## Installation

```stata
net install svartools, from("https://raw.githubusercontent.com/YOUR_USERNAME/svartools/main/")
```

To update later:

```stata
ado update svartools
```

To uninstall:

```stata
ado uninstall svartools
```

## Commands

- `sim_svar` — simulate SVAR(p) time series via Mata
- `irf_dualband` — plot impulse responses with dual confidence bands

After installing, see `help sim_svar` and `help irf_dualband`.

## Example

```stata
matrix A   = (1, 0 \ 0.5, 1)
matrix Phi = (0.6, 0.1 \ 0.2, 0.7)
sim_svar, nvars(2) lags(1) coef(Phi) amat(A) nobs(500)

var y1 y2, lags(1)
svar y1 y2, lags(1) aeq(...) beq(...)
irf create svar1, set(myirf, replace) step(36)

irf_dualband, irffile("myirf.irf") impulse(y1) response(y2) ///
              irfname(svar1) normalize
```

## License

MIT
