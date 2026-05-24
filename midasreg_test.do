*===============================================================*
*  midasreg_test.do
*  生成混频数据并测试 midasreg 命令
*
*  数据生成过程(DGP):
*    y 为低频(季度)被解释变量
*    x1, x2 为高频(月度)解释变量,频率比 m = 3
*
*    y_t = 2 + 0.5*y_{t-1}
*        + sum_{j=0}^{11} beta1_j * x1_{3t-j}
*        + sum_{j=0}^{11} beta2_j * x2_{3t-j}
*        + e_t
*
*  权重函数 beta_j 满足归一化指数 Almon 多项式:
*    beta_j = lambda * exp(theta1*s + theta2*s^2) / sum(...)
*===============================================================*

clear all
set more off
set seed 20260524

*---------------------------------------------------------------*
* 1. 基本参数设定
*---------------------------------------------------------------*
local n_low    = 200           // 低频观测数(季度数)
local mratio   = 3             // 频率比:每季度 3 个月
local n_high   = `n_low' * `mratio'   // 高频观测总数 = 600
local K        = 12            // MIDAS 多项式中高频滞后数(覆盖近 4 季度)

display _newline "{txt}样本设定:"
display "  低频观测数(季度) = `n_low'"
display "  高频观测数(月度) = `n_high'"
display "  频率比 m         = `mratio'"
display "  MIDAS 滞后数 K   = `K'"

*---------------------------------------------------------------*
* 2. 创建月度面板(高频时间索引)
*---------------------------------------------------------------*
set obs `n_high'
gen int month_id = _n              // 月度时间索引 1..600
tsset month_id

* 同时记录所属季度编号(每 3 个月一个季度)
gen int quarter_id = ceil(month_id / `mratio')

* 标记每个季度的最后一个月(即低频观测落点)
gen byte last_in_q = mod(month_id, `mratio') == 0

*---------------------------------------------------------------*
* 3. 生成两个高频解释变量(月度,i.i.d. 标准正态)
*---------------------------------------------------------------*
gen double x1 = rnormal()
gen double x2 = rnormal()

*---------------------------------------------------------------*
* 4. 设定真实的 MIDAS 权重(归一化指数 Almon 多项式)
*    w(s) = exp(theta1*s + theta2*s^2) / sum_j exp(...)
*    s 从 1 到 K
*---------------------------------------------------------------*
* x1 的真实参数:正向影响,衰减较快
local lambda1  = 1.5
local theta1_1 = 0.6
local theta2_1 = -0.10

* x2 的真实参数:负向影响,衰减较慢
local lambda2  = -1.0
local theta1_2 = 0.3
local theta2_2 = -0.04

* 用 Mata 计算真实权重并保存到 Stata 矩阵
mata:
    K       = `K'
    s       = (1::K)
    // x1
    num1    = exp(`theta1_1':*s :+ `theta2_1':*(s:^2))
    w1      = `lambda1' :* num1 :/ sum(num1)
    // x2
    num2    = exp(`theta1_2':*s :+ `theta2_2':*(s:^2))
    w2      = `lambda2' :* num2 :/ sum(num2)
    st_matrix("true_w1", w1)
    st_matrix("true_w2", w2)
end

matrix list true_w1, title("x1 的真实 MIDAS 系数")
matrix list true_w2, title("x2 的真实 MIDAS 系数")

*---------------------------------------------------------------*
* 5. 在低频期末构造 MIDAS 加权和:
*    MIDAS_x1_t = sum_{j=0}^{K-1} w1[j+1] * x1_{3t-j}
*---------------------------------------------------------------*
gen double midas_x1 = 0
gen double midas_x2 = 0
forvalues j = 0/`=`K'-1' {
    local s = `j' + 1                    // s 从 1 开始
    local wj1 : display %20.15f true_w1[`s', 1]
    local wj2 : display %20.15f true_w2[`s', 1]
    qui replace midas_x1 = midas_x1 + (`wj1') * L`j'.x1
    qui replace midas_x2 = midas_x2 + (`wj2') * L`j'.x2
}

* 只在每个季度的最后一个月(low_in_q==1)保留 MIDAS 加权和;
* 其余月份置为缺失。
qui replace midas_x1 = . if !last_in_q
qui replace midas_x2 = . if !last_in_q

*---------------------------------------------------------------*
* 6. 生成低频被解释变量 y(只在每季度末有非缺失值)
*    y_t = 2 + 0.5*y_{t-1} + midas_x1 + midas_x2 + e_t
*
*    注:K = 12 个高频滞后需要至少 12 个月才能算出第一个 MIDAS 加权和,
*    即从第 4 个季度(month_id = 12)开始 midas_x1/x2 才有值,
*    所以 y 也从第 4 个季度开始生成,前 3 个季度的 y 保持缺失。
*---------------------------------------------------------------*
gen double y = .
gen double e_lf = .

* 仅在每个季度的末月生成噪声项
qui replace e_lf = rnormal() if last_in_q

* 找出第一个 MIDAS 加权和非缺失的季度
qui sum quarter_id if last_in_q & !missing(midas_x1) & !missing(midas_x2)
local q_start = r(min)
display _newline "{txt}第一个可用季度 = `q_start'"

* 起始季度的 y:用截距+midas+噪声(不含 AR 项)作为初始条件
qui replace y = 2 + midas_x1 + midas_x2 + e_lf ///
    if last_in_q & quarter_id == `q_start'

* 之后按递推方式生成:L`mratio'.y 即上一季度末的 y
local q_next = `q_start' + 1
forvalues q = `q_next'/`n_low' {
    qui replace y = 2 + 0.5*L`mratio'.y ///
                    + midas_x1 + midas_x2 + e_lf       ///
                  if last_in_q & quarter_id == `q'
}

*---------------------------------------------------------------*
* 7. 数据检查
*---------------------------------------------------------------*
display _newline "{txt}数据概览:"
summarize y x1 x2 midas_x1 midas_x2

display _newline "{txt}y 的非缺失观测(应等于 `n_low'):"
count if !missing(y)

display _newline "{txt}前 10 个季度末的 y 值:"
list quarter_id month_id y midas_x1 midas_x2 if last_in_q in 1/10, ///
    abbreviate(12)

*---------------------------------------------------------------*
* 8. 保存数据集(可选)
*---------------------------------------------------------------*
label var month_id     "月度时间索引"
label var quarter_id   "季度编号"
label var last_in_q    "是否为季度末月"
label var x1           "高频解释变量 1(月度)"
label var x2           "高频解释变量 2(月度)"
label var y            "低频被解释变量(季度)"
label var midas_x1     "x1 的真实 MIDAS 加权和"
label var midas_x2     "x2 的真实 MIDAS 加权和"

save midas_test_data.dta, replace
display _newline "{txt}数据已保存为 midas_test_data.dta"

*===============================================================*
*  测试 midasreg 命令
*===============================================================*

display _newline _newline ///
    "{hline 70}" _newline ///
    "{txt}开始测试 midasreg 命令" _newline ///
    "{hline 70}"

*---------------------------------------------------------------*
* 测试 1:单变量 MIDAS,默认指数 Almon 权重
*---------------------------------------------------------------*
display _newline "{txt}>>> 测试 1: 仅用 x1,nealmon 权重,含 1 阶季度 AR"
midasreg y, hfvar(x1) lags(0 11) mratio(3) wscheme(nealmon) ylags(1)

display _newline "{txt}真实参数: lambda1 = `lambda1', theta1_1 = `theta1_1', " ///
    "theta2_1 = `theta2_1'"
display "{txt}真实截距 = 2, 真实 AR(1) = 0.5"

*---------------------------------------------------------------*
* 测试 2:加入第二个高频变量 x2 作为低频均值(对比基准)
* 注:本程序的 hfvar() 一次只接受一个高频变量。
* 真正的多高频变量 MIDAS 需要分别建模;这里仅演示 x1 的拟合,
* 同时把 midas_x2(已加权)当作低频回归量加入,作为"理想"对照。
*---------------------------------------------------------------*
display _newline "{txt}>>> 测试 2: x1 用 MIDAS,x2 加权后作为低频协变量"
midasreg y midas_x2, hfvar(x1) lags(0 11) mratio(3) ///
    wscheme(nealmon) ylags(1)

*---------------------------------------------------------------*
* 测试 3:换用 Beta 权重
*---------------------------------------------------------------*
display _newline "{txt}>>> 测试 3: x1 用 Beta 权重"
midasreg y, hfvar(x1) lags(0 11) mratio(3) wscheme(beta) ylags(1)

*---------------------------------------------------------------*
* 测试 4:无约束 MIDAS (U-MIDAS)——每个高频滞后单独估一个系数
*---------------------------------------------------------------*
display _newline "{txt}>>> 测试 4: U-MIDAS(无约束)"
midasreg y, hfvar(x1) lags(0 11) mratio(3) wscheme(umidas) ylags(1)

display _newline "{txt}对比:U-MIDAS 各系数 vs. 真实 MIDAS 权重"
matrix list true_w1, title("x1 的真实 MIDAS 系数")

*---------------------------------------------------------------*
* 测试 5:稳健标准误
*---------------------------------------------------------------*
display _newline "{txt}>>> 测试 5: nealmon 权重 + 稳健标准误"
midasreg y, hfvar(x1) lags(0 11) mratio(3) ///
    wscheme(nealmon) ylags(1) robust

*---------------------------------------------------------------*
* 测试 6:二阶 Almon 多项式
*---------------------------------------------------------------*
display _newline "{txt}>>> 测试 6: 二阶 Almon 多项式"
midasreg y, hfvar(x1) lags(0 11) mratio(3) wscheme(almon) ylags(1)

*---------------------------------------------------------------*
* 测试 7:不含 AR 项的纯 DL-MIDAS
*---------------------------------------------------------------*
display _newline "{txt}>>> 测试 7: 无 AR 项的 DL-MIDAS"
midasreg y, hfvar(x1) lags(0 11) mratio(3) wscheme(nealmon)

*---------------------------------------------------------------*
* 完成
*---------------------------------------------------------------*
display _newline "{hline 70}"
display "{txt}全部测试完成"
display "{hline 70}"
