capture program drop irf_dualband
program define irf_dualband
    version 14
    syntax , Irffile(string) Impulse(string) Response(string) ///
            Irfname(string) ///
            [Type(string) ///
             Inner(integer 68) Outer(integer 90) ///
             NORMalize ///
             Title(string) Ytitle(string) Xtitle(string) ///
             INNERcolor(string) OUTERcolor(string) ///
             INNERalpha(integer 60) OUTERalpha(integer 40) ///
             Name(string) ///
			 Export(string)]

    if "`type'"       == "" local type       sirf
    if "`xtitle'"     == "" local xtitle     "Steps"
    if "`ytitle'"     == "" local ytitle     "Response"
    if "`title'"      == "" local title      "`impulse' → `response'"
    if "`innercolor'" == "" local innercolor "139 0 0"
    if "`outercolor'" == "" local outercolor "219 149 168"
	if "`name'" == "" local name "`impulse'_`response'"
	
    local z_inner = invnormal(1 - (1 - `inner'/100)/2)
    local z_outer = invnormal(1 - (1 - `outer'/100)/2)

    preserve
    use "`irffile'", clear

    capture confirm variable `type'
    if _rc {
        di as error "IRF 文件没有 `type' 列"
        ds *irf*
        restore
        exit 198
    }
    capture confirm variable std`type'
    if _rc {
        di as error "IRF 文件没有 std`type' 标准误列"
        restore
        exit 198
    }

    *--- 算 normalize 基准:impulse 对自己 step 0 的响应 ---*
    local scale = 1
    if "`normalize'" != "" {
        quietly summarize `type' if irfname == "`irfname'" &              ///
                                    impulse == "`impulse'" &              ///
                                    response == "`impulse'" &             ///
                                    step    == 0
        if r(N) == 0 {
            di as error "找不到 normalize 基准:`impulse' → `impulse' at step 0"
            restore
            exit 198
        }
        local scale = r(mean)
        if abs(`scale') < 1e-10 {
            di as error "normalize 基准值接近 0,无法 rescale"
            restore
            exit 198
        }
        di as text "Normalized by `impulse' own response at step 0 = `scale'"
    }

    *--- 筛选 ---*
    quietly keep if irfname == "`irfname'" & ///
                    impulse == "`impulse'" & response == "`response'"
    if _N == 0 {
        di as error "找不到匹配数据"
        restore
        exit 198
    }

    sort step
    quietly {
        gen pe       = `type' / `scale'
        gen std_adj  = std`type' / abs(`scale')
        gen lo_inner = pe - `z_inner' * std_adj
        gen up_inner = pe + `z_inner' * std_adj
        gen lo_outer = pe - `z_outer' * std_adj
        gen up_outer = pe + `z_outer' * std_adj
    }

    twoway (rarea lo_outer up_outer step,                                      ///
                color("`outercolor'%`outeralpha'") lwidth(none))               ///
           (rarea lo_inner up_inner step,                                      ///
                color("`innercolor'%`inneralpha'") lwidth(none))               ///
           (line pe step, lcolor(black) lwidth(medthick))                      ///
           (scatter pe step, mcolor(black) msymbol(Oh) msize(small)            ///
                mlcolor(black) mfcolor(white)),                                ///
           yline(0, lcolor(black) lwidth(thin))                                ///
           legend(off)                                                         ///
           title(`"`title'"', size(medium))                                    ///
           ytitle(`"`ytitle'"') xtitle(`"`xtitle'"')                           ///
           graphregion(color(white)) plotregion(color(white))                  ///
           xlabel(, grid glcolor(gs14) glpattern(dash))                        ///
           ylabel(, grid glcolor(gs14) glpattern(dash))							///
		   name(`name',replace)

    if "`export'" != "" {
        graph export "`export'", replace
    }

    restore
end