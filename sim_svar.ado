*=============================================================*
* sim_svar: 用 Mata 模拟 SVAR(p)
*
* 结构式: A y_t = c + Phi_1 y_{t-1} + ... + Phi_p y_{t-p} + B eps_t
*         eps_t ~ N(0, I_n)
*
* 必选:
*   nvars(#)     - 变量个数 n
*   lags(#)      - 滞后阶数 p
*   coef(matname)- n x (n*p) 简化式滞后系数 [Phi_1 Phi_2 ... Phi_p]
*
* 可选:
*   amat(matname)  - n x n 结构矩阵 A,默认 I
*   bmat(matname)  - n x n 结构矩阵 B,默认 I
*   const(matname) - n x 1 常数向量,默认 0
*   nobs(#)        - 保留观测数,默认 500
*   burnin(#)      - 烧入期,默认 200
*   seed(#)        - 随机种子,默认 12345
*   names(string)  - 变量名,默认 y1 y2 ...
*=============================================================*

capture program drop sim_svar
capture mata: mata drop _sim_svar()

program define sim_svar
    version 14
    syntax , Nvars(integer) Lags(integer) Coef(string)       ///
            [Amat(string) Bmat(string) Const(string)         ///
             Nobs(integer 500) Burnin(integer 200)           ///
             Seed(integer 12345) Names(string)]

    if `nvars' < 1 | `lags' < 1 {
        di as error "nvars 与 lags 必须为正整数"
        exit 198
    }

    capture confirm matrix `coef'
    if _rc {
        di as error "矩阵 `coef' 不存在"
        exit 198
    }
    local expcols = `nvars' * `lags'
    if rowsof(`coef') != `nvars' | colsof(`coef') != `expcols' {
        di as error "coef 应为 `nvars' × `expcols'"
        exit 198
    }

    if "`amat'" == "" {
        capture matrix drop _simsvar_A
        matrix _simsvar_A = I(`nvars')
        local Aname _simsvar_A
    }
    else {
        capture confirm matrix `amat'
        if _rc {
            di as error "amat `amat' 不存在"
            exit 198
        }
        local Aname `amat'
    }

    if "`bmat'" == "" {
        capture matrix drop _simsvar_B
        matrix _simsvar_B = I(`nvars')
        local Bname _simsvar_B
    }
    else {
        capture confirm matrix `bmat'
        if _rc {
            di as error "bmat `bmat' 不存在"
            exit 198
        }
        local Bname `bmat'
    }

    if "`const'" == "" {
        capture matrix drop _simsvar_C
        matrix _simsvar_C = J(`nvars', 1, 0)
        local Cname _simsvar_C
    }
    else {
        capture confirm matrix `const'
        if _rc {
            di as error "const `const' 不存在"
            exit 198
        }
        local Cname `const'
    }

    if "`names'" == "" {
        local varnames ""
        forvalues j = 1/`nvars' {
            local varnames "`varnames' y`j'"
        }
        local varnames = trim("`varnames'")
    }
    else {
        local varnames "`names'"
        local nc : word count `varnames'
        if `nc' != `nvars' {
            di as error "names 个数(`nc') ≠ nvars(`nvars')"
            exit 198
        }
    }

    set seed `seed'
    local total = `nobs' + `burnin'
    clear
    set obs `total'
    gen t = _n
    foreach v of local varnames {
        quietly gen double `v' = .
    }

    mata: _sim_svar("`Aname'", "`Bname'", "`coef'", "`Cname'", ///
                   `nvars', `lags', `total', "`varnames'")

    quietly drop if t <= `burnin'
    quietly replace t = _n
    tsset t

    di as result _n "SVAR(`lags') 模拟完成:`nvars' 变量,`nobs' 观测"
    di as text "变量名:`varnames'"
end


mata:
void _sim_svar(string scalar Aname,
               string scalar Bname,
               string scalar Phiname,
               string scalar Cname,
               real   scalar n,
               real   scalar p,
               real   scalar T,
               string scalar vnames)
{
    real matrix      A, B, Phi, Ainv, M, eps, Y, F
    real colvector   C, cred, ylag
    real rowvector   ev_abs
    real scalar      i, k, maxev, np
    string rowvector vlist

    A   = st_matrix(Aname)
    B   = st_matrix(Bname)
    Phi = st_matrix(Phiname)
    C   = st_matrix(Cname)

    printf("{txt}A: %g x %g, B: %g x %g, Phi: %g x %g, C: %g x %g\n",
           rows(A), cols(A), rows(B), cols(B),
           rows(Phi), cols(Phi), rows(C), cols(C))

    np = n * p
    F  = J(np, np, 0)
    F[(1::n), .] = Phi
    if (p > 1) {
        F[((n+1)::np), (1::(np-n))] = I(np - n)
    }
    ev_abs = abs(eigenvalues(F))
    maxev  = max(ev_abs)
    if (maxev >= 1) {
        printf("{txt}警告:伴随矩阵最大特征根模 = %9.4f ≥ 1\n", maxev)
    }

    Ainv = luinv(A)
    M    = Ainv * B
    cred = Ainv * C
    eps  = rnormal(T, n, 0, 1)

    Y = J(T, n, 0)
    for (i = p+1; i <= T; i++) {
        ylag = J(0, 1, .)
        for (k = 1; k <= p; k++) {
            ylag = ylag \ Y[i-k, .]'
        }
        Y[i, .] = (cred + Phi * ylag + M * eps[i, .]')'
    }

    vlist = tokens(vnames)
    st_store(., vlist, Y)
}
end