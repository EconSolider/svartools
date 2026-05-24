*! midasreg v1.0.1  24may2026
*! Mixed Data Sampling (MIDAS) regression
*! Implements MIDAS regression following Ghysels, Santa-Clara & Valkanov (2002)
*! "The MIDAS Touch: Mixed Data Sampling Regression Models" and
*! Ghysels, Kvedaras & Zemlys (2016) JSS 72(4), "Mixed Frequency Data
*! Sampling Regression Models: The R Package midasr".
*!
*! Author: generated for Stata
*! Distribution: free for academic use

program define midasreg, eclass
    version 14.0

    /*------------------------------------------------------------------*/
    /*  Parse command syntax                                            */
    /*    midasreg depvar [lf_indepvars] ,                              */
    /*        HFvar(varname) Lags(numlist) Mratio(integer)              */
    /*        [ Weight(string) YLags(numlist) NORMalize                 */
    /*          INITial(numlist) ITERate(integer) TOLerance(real)       */
    /*          NOConstant Robust noLOG ]                               */
    /*------------------------------------------------------------------*/

    syntax varlist(min=1 numeric ts) [if] [in] ,    ///
        HFvar(varname numeric)                       ///
        Lags(numlist integer >=0 sort)               ///
        Mratio(integer)                              ///
        [ Wscheme(string)                            ///
          YLags(numlist integer >0 sort)             ///
          INITial(numlist)                           ///
          ITERate(integer 200)                       ///
          TOLerance(real 1e-7)                       ///
          NOConstant                                 ///
          Robust                                     ///
          noLOG ]

    /* depvar = first variable, lf_indepvars = remaining LF regressors */
    tokenize `varlist'
    local depvar `1'
    macro shift
    local lfvars `*'

    marksample touse
    markout `touse' `depvar' `lfvars'

    /* internal name 'weight' kept for code readability */
    local weight "`wscheme'"

    /*--- default weight scheme ---*/
    if "`weight'" == "" local weight "nealmon"
    local weight = lower("`weight'")
    if !inlist("`weight'","nealmon","beta","betann","almon","umidas") {
        di as err "wscheme() must be one of: nealmon, beta, betann, almon, umidas"
        exit 198
    }

    /*--- validate Lags() ---*/
    local nlags : word count `lags'
    if `nlags' < 2 {
        di as err "lags() must contain at least 2 integers (kmin kmax)"
        exit 198
    }
    local kmin : word 1 of `lags'
    local kmax : word `nlags' of `lags'
    local K = `kmax' - `kmin' + 1
    if `K' < 2 {
        di as err "lags() span must be >= 2 high-frequency lags"
        exit 198
    }
    if `mratio' < 1 {
        di as err "mratio() must be a positive integer (frequency ratio)"
        exit 198
    }

    /*--- AR component (low frequency lags of y) ---*/
    local nylags = 0
    if "`ylags'" != "" {
        local nylags : word count `ylags'
    }

    /*--- decide number of weight parameters & defaults ---*/
    if "`weight'" == "nealmon" {
        /* Normalized exponential Almon: theta + 2 polynomial coeffs */
        local nwpar = 2
        local wlabel "Normalized exponential Almon"
    }
    else if "`weight'" == "beta" {
        /* Normalized beta (zero-end), 2 shape params */
        local nwpar = 2
        local wlabel "Normalized beta (zero last weight)"
    }
    else if "`weight'" == "betann" {
        /* Normalized beta (non-zero end), 3 params */
        local nwpar = 3
        local wlabel "Normalized beta (non-zero last weight)"
    }
    else if "`weight'" == "almon" {
        /* Unrestricted Almon polynomial of order 2: 3 coefficients */
        local nwpar = 3
        local wlabel "Almon polynomial (order 2)"
    }
    else if "`weight'" == "umidas" {
        /* Unrestricted MIDAS: one beta per HF lag */
        local nwpar = `K'
        local wlabel "Unrestricted MIDAS (U-MIDAS)"
    }

    /*--- build the high-frequency lag matrix using Stata time-series ops.
          The HF variable must be stored in long format aligned with the
          high-frequency time index; the data must be tsset.            */

    /*--- collect the names of generated HF lag variables ---*/
    local hfvarlist ""
    forvalues j = `kmin'/`kmax' {
        tempvar hfl`j'
        qui gen double `hfl`j'' = L`j'.`hfvar'
        local hfvarlist `hfvarlist' `hfl`j''
    }

    /*--- low frequency AR lags of depvar (if any).
          ylags() are given in low-frequency units, but the data are
          tsset on the high-frequency time index, so we have to multiply
          by mratio to get the corresponding HF lag.                   */
    local ylagvars ""
    if `nylags' > 0 {
        foreach p of local ylags {
            local hfp = `p' * `mratio'
            tempvar yl`p'
            qui gen double `yl`p'' = L`hfp'.`depvar'
            local ylagvars `ylagvars' `yl`p''
        }
    }

    /*--- update touse to drop missing rows after lagging ---*/
    markout `touse' `hfvarlist' `ylagvars' `lfvars'

    /*--- count valid observations ---*/
    qui count if `touse'
    local N = r(N)
    if `N' == 0 {
        di as err "no observations"
        exit 2000
    }

    /*--- handle constant ---*/
    local hascons = cond("`noconstant'"=="", 1, 0)

    /*--- count low-frequency regressors ---*/
    local nlf : word count `lfvars'

    /*--- total number of parameters (for U-MIDAS this is K; otherwise
         we have an impact (slope) parameter + nwpar weight pars.  
         For unrestricted Almon, weight already absorbs scale, so no  
         extra slope. We follow midasr convention: a single scale     
         (impact) coefficient + weight shape parameters, except       
         U-MIDAS and Almon which directly fit one coef per lag/poly.   */

    if "`weight'"=="umidas" | "`weight'"=="almon" {
        local nslope = 0
    }
    else {
        local nslope = 1
    }
    local nmpar = `nslope' + `nwpar'        // MIDAS-block params
    local nallpar = `hascons' + `nlf' + `nylags' + `nmpar'

    /*--- initial values ---*/
    if "`initial'" != "" {
        local ninit : word count `initial'
        if `ninit' != `nallpar' {
            di as err "initial() requires `nallpar' values; got `ninit'"
            exit 198
        }
        local b0 "`initial'"
    }
    else {
        /* default starting values */
        local b0 ""
        if `hascons' local b0 "0 "
        forvalues k = 1/`nlf' {
            local b0 "`b0' 0 "
        }
        forvalues k = 1/`nylags' {
            local b0 "`b0' 0 "
        }
        /* slope */
        if `nslope' == 1 local b0 "`b0' 1 "
        /* weight params -- defaults */
        if "`weight'" == "nealmon" {
            local b0 "`b0' 1 -0.5"
        }
        else if "`weight'" == "beta" {
            local b0 "`b0' 1 5"
        }
        else if "`weight'" == "betann" {
            local b0 "`b0' 1 5 0"
        }
        else if "`weight'" == "almon" {
            local b0 "`b0' 0 0 0"
        }
        else if "`weight'" == "umidas" {
            forvalues k = 1/`K' {
                local b0 "`b0' 0 "
            }
        }
    }

    /*--- estimation: hand off to Mata NLS engine ---*/
    local dorob = ("`robust'" != "")
    local verbose = ("`log'" != "nolog")

    tempname bvec Vmat
    tempname sN srss ss2 srmse sdfr sk
    mata: _midas_fit("`depvar'", "`lfvars'", "`ylagvars'", "`hfvarlist'", ///
        "`touse'", "`weight'", `mratio', `K', `kmin', `nwpar', `nslope', ///
        `hascons', `iterate', `tolerance', `dorob', "`b0'",              ///
        "`bvec'", "`Vmat'", `verbose')

    /* capture summary scalars Mata stuffed into r() */
    scalar `sN'    = r(_N)
    scalar `srss'  = r(_rss)
    scalar `ss2'   = r(_sigma2)
    scalar `srmse' = r(_rmse)
    scalar `sdfr'  = r(_df_r)
    scalar `sk'    = r(_k)

    /*--- build coefficient labels for ereturn ---*/
    /* bvec is already a 1xk row vector from Mata's st_matrix(bname, b') */
    local cnames ""
    if `hascons' local cnames "`cnames' _cons"
    foreach v of local lfvars {
        local cnames "`cnames' `v'"
    }
    if `nylags' > 0 {
        foreach p of local ylags {
            local cnames "`cnames' L`p'.`depvar'"
        }
    }

    /* MIDAS block coefficient names */
    if `nslope' == 1 local cnames "`cnames' `hfvar'_slope"
    if "`weight'"=="nealmon" {
        local cnames "`cnames' `hfvar'_theta1 `hfvar'_theta2"
    }
    else if "`weight'"=="beta" {
        local cnames "`cnames' `hfvar'_a `hfvar'_b"
    }
    else if "`weight'"=="betann" {
        local cnames "`cnames' `hfvar'_a `hfvar'_b `hfvar'_c"
    }
    else if "`weight'"=="almon" {
        local cnames "`cnames' `hfvar'_a0 `hfvar'_a1 `hfvar'_a2"
    }
    else if "`weight'"=="umidas" {
        forvalues j = `kmin'/`kmax' {
            local cnames "`cnames' `hfvar'_b`j'"
        }
    }

    matrix colnames `bvec' = `cnames'
    matrix colnames `Vmat' = `cnames'
    matrix rownames `Vmat' = `cnames'

    /*--- post results ---*/
    ereturn post `bvec' `Vmat', depname(`depvar') obs(`N') esample(`touse')
    ereturn local cmd        "midasreg"
    ereturn local depvar     "`depvar'"
    ereturn local hfvar      "`hfvar'"
    ereturn local weight     "`weight'"
    ereturn local wlabel     "`wlabel'"
    ereturn local lfvars     "`lfvars'"
    ereturn local ylags      "`ylags'"
    ereturn local indepvars  "`lfvars'"
    ereturn local vcetype    = cond(`dorob', "Robust", "OLS")
    ereturn scalar mratio    = `mratio'
    ereturn scalar kmin      = `kmin'
    ereturn scalar kmax      = `kmax'
    ereturn scalar K         = `K'
    ereturn scalar nwpar     = `nwpar'
    ereturn scalar nslope    = `nslope'
    ereturn scalar nylags    = `nylags'
    ereturn scalar hascons   = `hascons'
    ereturn scalar rss       = `srss'
    ereturn scalar sigma2    = `ss2'
    ereturn scalar rmse      = `srmse'
    ereturn scalar df_r      = `sdfr'
    ereturn scalar k         = `sk'
    ereturn local title      "Mixed Data Sampling (MIDAS) regression"

    /*--- pretty display ---*/
    Display
end

/*=====================================================================*/
program define Display
    di
    di as txt "{hline 78}"
    di as txt "Mixed Data Sampling (MIDAS) regression"
    di as txt "  Dependent variable : " as res "`e(depvar)'"
    di as txt "  HF variable        : " as res "`e(hfvar)'" ///
        as txt "   (frequency ratio m = " as res e(mratio) as txt ")"
    di as txt "  HF lag span        : " as res e(kmin) as txt " ... " ///
        as res e(kmax) as txt "  (K = " as res e(K) as txt ")"
    di as txt "  Weighting scheme   : " as res "`e(weight)'" ///
        as txt "  (`e(wlabel)')"
    di as txt "  Observations       : " as res e(N)
    di as txt "{hline 78}"
    ereturn display
end

/*=====================================================================*/
/*                              MATA SECTION                            */
/*=====================================================================*/

version 14.0
mata:
mata set matastrict off

/*-------- Weight functions --------*/
/* All return a column vector of length K of weights ; b is the
   high-frequency lag index s = 1,...,K.  Functions are normalized
   so that the weights sum to one (except 'almon' which is the raw
   polynomial, and 'umidas' which is identity). */

real colvector _midas_w_nealmon(real vector theta, real scalar K) {
    /* w_s = exp(theta1*s + theta2*s^2) / sum_j exp(theta1*j+theta2*j^2)
       Subtract the max of the exponent before exp() for numerical
       stability (softmax-style); the ratio is invariant to that shift. */
    real colvector s, lin, num
    real scalar lmax
    s   = (1::K)
    lin = theta[1]:*s :+ theta[2]:*(s:^2)
    lmax = max(lin)
    if (lmax == . | missing(lmax)) lmax = 0
    num = exp(lin :- lmax)
    return(num :/ sum(num))
}

real colvector _midas_w_beta(real vector theta, real scalar K) {
    /* Normalized beta with last weight forced to zero (Ghysels Beta).
       theta = (a, b),  a>0, b>0.
       x_s = (s-1)/(K-1) ; small xi to avoid 0 and 1.
       Compute in log domain then exp(stable shift) for safety. */
    real colvector s, x, lognum, num
    real scalar xi, a, b, lmax
    xi = 1e-8
    s  = (1::K)
    x  = xi :+ (1-2*xi) :* (s :- 1) :/ (K - 1)
    a  = max((theta[1], 1e-3))
    b  = max((theta[2], 1e-3))
    lognum = (a-1) :* log(x) :+ (b-1) :* log(1:-x)
    lmax = max(lognum)
    if (lmax == . | missing(lmax)) lmax = 0
    num = exp(lognum :- lmax)
    return(num :/ sum(num))
}

real colvector _midas_w_betann(real vector theta, real scalar K) {
    /* Normalized beta with non-zero last weight: w_s = (psi_s + c)/sum  */
    real colvector s, x, lognum, num
    real scalar xi, a, b, lmax
    xi = 1e-8
    s  = (1::K)
    x  = xi :+ (1-2*xi) :* (s :- 1) :/ (K - 1)
    a  = max((theta[1], 1e-3))
    b  = max((theta[2], 1e-3))
    lognum = (a-1) :* log(x) :+ (b-1) :* log(1:-x)
    lmax = max(lognum)
    if (lmax == . | missing(lmax)) lmax = 0
    num = exp(lognum :- lmax) :+ theta[3]
    /* protect against negative weights when theta[3] is very negative */
    num = num :- min(num) :+ 1e-10
    return(num :/ sum(num))
}

real colvector _midas_w_almon(real vector a, real scalar K) {
    /* Raw Almon polynomial of order 2:  b_s = a0 + a1*s + a2*s^2.
       No normalization: each parameter is identified directly.       */
    real colvector s
    s = (1::K)
    return(a[1] :+ a[2]:*s :+ a[3]:*(s:^2))
}

real colvector _midas_w_umidas(real vector a, real scalar K) {
    /* return a column vector regardless of input orientation */
    if (cols(a) > rows(a)) return(a')
    return(a)
}

/*-------- dispatch on weight name --------*/
real colvector _midas_weights(string scalar wname, real vector theta,
                              real scalar K) {
    if (wname=="nealmon") return(_midas_w_nealmon(theta, K))
    if (wname=="beta")    return(_midas_w_beta(theta, K))
    if (wname=="betann")  return(_midas_w_betann(theta, K))
    if (wname=="almon")   return(_midas_w_almon(theta, K))
    if (wname=="umidas")  return(_midas_w_umidas(theta, K))
    _error("unknown weight scheme: " + wname)
}

/*-------- residuals: data layout
   Z  : N x p   matrix of LF regressors stacked as
        [ const | lfvars | yLagVars ]            (p columns)
   X  : N x K   matrix of HF lag values (already extracted)
   b  : full coefficient vector
        first p entries are LF coefficients
        next nslope=0 or 1 entry is the MIDAS scale lambda
        last nwpar entries are weight-function parameters             */

real colvector _midas_resid(real colvector b, real colvector y,
                            real matrix Z, real matrix X,
                            string scalar wname,
                            real scalar K, real scalar nwpar,
                            real scalar nslope, real scalar p) {
    real colvector beta_lf, theta, w, yhat, mid
    real scalar lambda

    if (nslope==1) {
        theta  = b[(p+2)::(p+1+nwpar)]
        w      = _midas_weights(wname, theta, K)
        lambda = b[p+1]
        mid    = lambda :* (X*w)
    }
    else {
        theta = b[(p+1)::(p+nwpar)]
        w     = _midas_weights(wname, theta, K)
        mid   = X*w
    }
    if (p > 0) {
        beta_lf = b[1::p]
        yhat = Z*beta_lf + mid
    }
    else {
        yhat = mid
    }
    return(y :- yhat)
}

real scalar _midas_sse(real colvector b, real colvector y,
                       real matrix Z, real matrix X,
                       string scalar wname,
                       real scalar K, real scalar nwpar,
                       real scalar nslope, real scalar p) {
    real colvector e
    e = _midas_resid(b, y, Z, X, wname, K, nwpar, nslope, p)
    return(e' * e)
}

/*-------- numerical Jacobian of residuals wrt b --------*/
real matrix _midas_J(real colvector b, real colvector y,
                     real matrix Z, real matrix X,
                     string scalar wname,
                     real scalar K, real scalar nwpar,
                     real scalar nslope, real scalar p) {
    real matrix Jac
    real colvector bp, bm, ep, em
    real scalar i, h, npar, sqrteps
    npar = rows(b)
    sqrteps = sqrt(epsilon(1))           /* ~ 1.49e-8 */
    Jac = J(rows(y), npar, .)
    for (i=1; i<=npar; i++) {
        /* Press et al. (Numerical Recipes) recommendation */
        h = sqrteps * max((abs(b[i]), 1))
        bp = b ; bm = b
        bp[i] = b[i] + h ; bm[i] = b[i] - h
        h = bp[i] - b[i]                 /* actual step after roundoff */
        ep = _midas_resid(bp, y, Z, X, wname, K, nwpar, nslope, p)
        em = _midas_resid(bm, y, Z, X, wname, K, nwpar, nslope, p)
        Jac[,i] = (ep :- em) :/ (2*h)
    }
    /* The Jacobian of e = y - f wrt b is -df/db, so Jac above is -df/db */
    return(Jac)
}

/*-------- Levenberg-Marquardt nonlinear least squares --------*/
real colvector _midas_lm(real colvector b0, real colvector y,
                         real matrix Z, real matrix X,
                         string scalar wname,
                         real scalar K, real scalar nwpar,
                         real scalar nslope, real scalar p,
                         real scalar maxit, real scalar tol,
                         real scalar verbose) {
    real colvector b, e, g, delta, b_new, e_new, b_best
    real matrix Jac, JJ, A
    real scalar lambda, sse, sse_new, it, npar, gnorm, ratio, accepted
    real scalar sse_best, maxstep, dnorm, scale
    b = b0
    e = _midas_resid(b, y, Z, X, wname, K, nwpar, nslope, p)
    sse = e' * e
    sse_best = sse
    b_best = b
    lambda = 1e-3
    npar = rows(b)
    /* trust-region radius: cap each L-M step's L2 norm to this initially. */
    maxstep = 2.0
    if (verbose) {
        printf("\n{txt}Levenberg-Marquardt iterations:\n")
        printf("{txt}  iter      sse          lambda      |grad|     |step|\n")
        printf("{txt} %4.0f   %12.6g   %10.4g\n", 0, sse, lambda)
    }
    for (it=1; it<=maxit; it++) {
        Jac = _midas_J(b, y, Z, X, wname, K, nwpar, nslope, p)
        if (hasmissing(Jac)) {
            if (verbose) printf("{err}  Jacobian has missing values at iteration %f; stopping\n", it)
            break
        }
        JJ = Jac' * Jac
        g  = Jac' * e
        gnorm = sqrt(g' * g)
        if (gnorm < 1e-8 * (1 + sqrt(sse))) {
            if (verbose) printf("{txt} %4.0f   %12.6g   %10.4g   %10.4g  (gradient small)\n", it, sse, lambda, gnorm)
            break
        }
        accepted = 0
        /* try up to 12 lambda increases per iteration */
        for (ratio = 0; ratio < 12; ratio++) {
            A = JJ + lambda * (diag(diagonal(JJ)) + 1e-10 * I(npar))
            delta = -lusolve(A, g)
            if (hasmissing(delta)) {
                lambda = lambda * 10
                if (lambda > 1e+15) break
                continue
            }
            /* trust-region: clip ||delta|| <= maxstep */
            dnorm = sqrt(delta' * delta)
            if (dnorm > maxstep) {
                scale = maxstep / dnorm
                delta = delta * scale
                dnorm = maxstep
            }
            b_new = b + delta
            e_new = _midas_resid(b_new, y, Z, X, wname, K, nwpar, nslope, p)
            if (hasmissing(e_new)) {
                lambda = lambda * 10
                maxstep = maxstep / 2
                if (lambda > 1e+15 | maxstep < 1e-12) break
                continue
            }
            sse_new = e_new' * e_new
            if (sse_new < sse) {
                /* accept the step */
                if (sse_new < sse_best) {
                    sse_best = sse_new
                    b_best = b_new
                }
                if (verbose & (it<=3 | mod(it,10)==0)) {
                    printf("{txt} %4.0f   %12.6g   %10.4g   %10.4g   %10.4g\n", it, sse_new, lambda, gnorm, dnorm)
                }
                if (abs(sse - sse_new) < tol*(abs(sse)+tol)) {
                    b = b_new ; e = e_new ; sse = sse_new
                    if (verbose) {
                        printf("{txt} %4.0f   %12.6g   %10.4g   %10.4g   %10.4g  (converged)\n", it, sse_new, lambda, gnorm, dnorm)
                    }
                    return(b)
                }
                b = b_new ; e = e_new ; sse = sse_new
                lambda = max((lambda/10, 1e-15))
                maxstep = min((maxstep * 2, 100))
                accepted = 1
                break
            }
            /* step did not lower sse */
            lambda = lambda * 10
            maxstep = maxstep / 2
            if (lambda > 1e+15 | maxstep < 1e-12) break
        }
        if (accepted == 0) {
            if (verbose) printf("{txt}  could not find a downhill step at iteration %f (best sse so far = %g)\n", it, sse_best)
            break
        }
    }
    if (verbose & it>=maxit) {
        printf("{err}  warning: maximum iterations (%f) reached\n", maxit)
    }
    return(b_best)
}

/*-------- main fit routine called from ado --------*/
void _midas_fit(string scalar dvname, string scalar lfnames,
                string scalar ynames,  string scalar hfnames,
                string scalar touse,   string scalar wname,
                real   scalar mratio,  real scalar K, real scalar kmin,
                real scalar nwpar,     real scalar nslope,
                real scalar hascons,   real scalar maxit,
                real scalar tol,       real scalar dorob,
                string scalar b0str,
                string scalar bname,   string scalar vname,
                real scalar verbose) {
    real colvector y, b0, b, e
    real matrix    Z, X, V, Jac, XtX_inv, JtJ, Je
    real scalar    N, p, sigma2

    /* Pull data */
    y = st_data(., dvname, touse)
    N = rows(y)

    /* Build Z = [const | lfvars | ylagvars] */
    Z = J(N, 0, .)
    if (hascons) Z = (Z, J(N,1,1))
    if (lfnames != "") {
        Z = (Z, st_data(., tokens(lfnames), touse))
    }
    if (ynames != "") {
        Z = (Z, st_data(., tokens(ynames), touse))
    }
    p = cols(Z)

    /* Build X (N x K) of high-frequency lags  */
    X = st_data(., tokens(hfnames), touse)
    if (cols(X) != K) {
        _error(sprintf("internal: HF lag matrix has %f cols, expected %f",
                       cols(X), K))
    }

    /* starting values */
    b0 = strtoreal(tokens(b0str))'
    if (rows(b0) != p + nslope + nwpar) {
        _error(sprintf("internal: starting vector length %f, expected %f",
                       rows(b0), p + nslope + nwpar))
    }

    /* run Levenberg-Marquardt */
    b = _midas_lm(b0, y, Z, X, wname, K, nwpar, nslope, p,
                  maxit, tol, verbose)

    /* residuals & variance at final b */
    e = _midas_resid(b, y, Z, X, wname, K, nwpar, nslope, p)
    if (hasmissing(e)) {
        printf("{err}  warning: estimated residuals contain missing values; results may be unreliable\n")
        e = editmissing(e, 0)
    }
    sigma2 = (e' * e) / (N - rows(b))

    Jac = _midas_J(b, y, Z, X, wname, K, nwpar, nslope, p)
    if (hasmissing(Jac)) {
        printf("{err}  warning: Jacobian at solution contains missing values; covariance matrix may be unreliable\n")
        Jac = editmissing(Jac, 0)
    }
    /* covariance:  Gauss-Newton (J'J)^{-1} * sigma2  ; or sandwich */
    JtJ = Jac' * Jac
    XtX_inv = invsym(JtJ + 1e-12 * I(rows(JtJ)))
    if (dorob) {
        /* HC0-type sandwich, computed without forming N x N diag(e^2) */
        Je = Jac :* (e:^2)
        V  = XtX_inv * (Jac' * Je) * XtX_inv
    }
    else {
        V = sigma2 :* XtX_inv
    }

    /* return to ado namespace */
    st_matrix(bname, b')
    st_matrix(vname, V)

    /* useful scalars communicated back through r() */
    st_numscalar("r(_N)",      N)
    st_numscalar("r(_rss)",    e' * e)
    st_numscalar("r(_sigma2)", sigma2)
    st_numscalar("r(_rmse)",   sqrt(sigma2))
    st_numscalar("r(_df_r)",   N - rows(b))
    st_numscalar("r(_k)",      rows(b))
}

end
