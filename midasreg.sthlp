{smcl}
{* *! version 1.0.0  24may2026}{...}
{title:标题}

{phang}
{bf:midasreg} {hline 2} 混频数据抽样(MIDAS)回归


{title:语法}

{p 8 17 2}
{cmd:midasreg} {it:被解释变量} [{it:低频自变量}] {ifin} {cmd:,}
{cmdab:hf:var(}{it:varname}{cmd:)}
{cmdab:l:ags(}{it:numlist}{cmd:)}
{cmdab:m:ratio(}{it:#}{cmd:)}
[{it:其它选项}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :模型设定(必选)}
{synopt :{cmdab:hf:var(}{it:varname}{cmd:)}}高频解释变量(以"长格式"存储,每个高频期一条观测){p_end}
{synopt :{cmdab:l:ags(}{it:kmin kmax}{cmd:)}}MIDAS 多项式中包含的高频滞后阶范围;例如 {cmd:lags(0 11)} 表示 12 个滞后{p_end}
{synopt :{cmdab:m:ratio(}{it:#}{cmd:)}}频率比 m(月度对季度填 3,月度对年度填 12,季度对年度填 4 等){p_end}

{syntab :模型细节}
{synopt :{cmdab:ws:cheme(}{it:scheme}{cmd:)}}权重函数:{cmd:nealmon}、{cmd:beta}、{cmd:betann}、{cmd:almon} 或 {cmd:umidas};默认 {cmd:nealmon}{p_end}
{synopt :{cmdab:yl:ags(}{it:numlist}{cmd:)}}加入被解释变量的低频滞后项(构成 ADL-MIDAS){p_end}
{synopt :{cmdab:noc:onstant}}不估计常数项{p_end}

{syntab :数值优化}
{synopt :{cmdab:init:ial(}{it:numlist}{cmd:)}}用户指定的初始值(顺序与输出结果表相同){p_end}
{synopt :{cmdab:iter:ate(}{it:#}{cmd:)}}Levenberg-Marquardt 迭代最大次数;默认 {cmd:iterate(200)}{p_end}
{synopt :{cmdab:tol:erance(}{it:#}{cmd:)}}残差平方和的收敛容忍度;默认 {cmd:tolerance(1e-7)}{p_end}
{synopt :{cmdab:nolog}}不显示迭代日志{p_end}

{syntab :标准误/推断}
{synopt :{cmdab:r:obust}}使用异方差稳健(HC0 三明治型)标准误{p_end}
{synoptline}

{p 4 6 2}
{it:被解释变量}须以低频观测(每个低频期一条非缺失值)。{it:低频自变量} 是可选的额外低频回归量。{cmd:hfvar()} 指定的高频变量每个低频期含 m 个观测,数据须先用 {cmd:tsset} 在高频时间索引上声明。具体数据布局见下文{it:数据布局}一节。


{title:命令说明}

{pstd}
{cmd:midasreg} 实现 Ghysels、Santa-Clara 与 Valkanov(2002)提出的混频数据抽样(MIDAS)回归。当被解释变量 {it:y_t} 以低频观测、而一个或多个解释变量以更高频率观测时,直接把高频数据加总到低频会丢失信息。MIDAS 的思路是保留全部高频观测,通过一个简约的参数化权重函数把众多高频滞后系数约束在一起,从而避免参数爆炸。

{pstd}
程序估计的完整模型是自回归分布滞后 MIDAS(ADL-MIDAS):

{p 8 17 2}
{it:y_t = mu + sum_p alpha_p y_{t-p} + Z_t'gamma + lambda * sum_{s=kmin}^{kmax} w(s; theta) x_{t,s} + e_t}

{pstd}
其中 {it:x_{t,s}} 是低频期 {it:t} 结束前 s 个高频期的高频变量观测,{it:Z_t} 为可选的低频回归量集合,{it:w(s; theta)} 是满足 sum_s w(s; theta) = 1 的权重函数({cmd:almon} 与 {cmd:umidas} 例外,见下)。

{pstd}
估计方法为非线性最小二乘(Mata 中用 Levenberg-Marquardt 算法实现)。


{title:选项说明}

{dlgtab:模型设定}

{phang}
{cmd:hfvar(}{it:varname}{cmd:)} 指定高频解释变量。该变量按"长格式"储存——每行一条高频观测——并且要与 {it:被解释变量} 在时间上对齐:某个低频期的最后 m 行就对应该低频期的最近一次低频观测。

{phang}
{cmd:lags(}{it:kmin kmax}{cmd:)} 给出 MIDAS 多项式所用的高频滞后阶范围。{cmd:lags(0 11)} 表示包含 {it:hfvar} 的当期值及其前 11 个高频滞后;若用月度数据预测季度变量({cmd:mratio(3)}),这就涵盖了最近 4 个季度的月度信息。至少需要 2 个滞后。

{phang}
{cmd:mratio(}{it:#}{cmd:)} 是频率比 m。月度对季度填 3,季度对年度填 4,月度对年度填 12,日度对月度(工作日)填 22 等。

{dlgtab:模型细节}

{phang}
{cmd:wscheme(}{it:scheme}{cmd:)} 选择用来约束 K = kmax-kmin+1 个高频滞后系数的参数化权重函数 {it:w(s; theta)}:

{p 8 12 2}{cmd:nealmon} {hline 2} 归一化指数 Almon 多项式:{it:w_s = exp(theta1*s + theta2*s^2)/sum_j exp(...)};2 个参数;能描述多种衰减形态;默认选项。{p_end}

{p 8 12 2}{cmd:beta} {hline 2} 归一化 Beta 密度:{it:w_s ∝ x_s^(a-1) * (1-x_s)^(b-1)},其中 {it:x_s = (s-1)/(K-1)},末尾权重强制为 0;2 个参数 a>0、b>0。对单调或单峰衰减拟合极佳(Ghysels、Santa-Clara、Valkanov 2002)。{p_end}

{p 8 12 2}{cmd:betann} {hline 2} 与 {cmd:beta} 相同,但增加第三个加性参数 c,末尾权重可不为 0。{p_end}

{p 8 12 2}{cmd:almon} {hline 2} 二阶非归一化 Almon 多项式 {it:b_s = a0 + a1*s + a2*s^2};3 个参数;不归一化,不另估 slope。{p_end}

{p 8 12 2}{cmd:umidas} {hline 2} 无约束 MIDAS:每个高频滞后估一个独立系数(Foroni、Marcellino、Schumacher 2015)。等价于对高频滞后矩阵做 OLS;当 K 较小时推荐使用。{p_end}

{phang}
{cmd:ylags(}{it:numlist}{cmd:)} 把被解释变量的低频滞后加入右端,使模型由 DL-MIDAS 变为 ADL-MIDAS。例如 {cmd:ylags(1 2)} 表示加入 {it:y_{t-1}} 与 {it:y_{t-2}}。注意此处用低频单位计阶,程序内部会乘以 {cmd:mratio} 转成高频滞后。

{phang}
{cmd:noconstant} 不估计常数项 {it:mu}。

{dlgtab:数值优化}

{phang}
{cmd:initial(}{it:numlist}{cmd:)} 提供全部系数的初始值。顺序为:常数(若有)→ 低频回归量(按输入顺序)→ AR 滞后(按 {cmd:ylags()} 顺序)→ MIDAS 主系数 lambda(仅 {cmd:nealmon}/{cmd:beta}/{cmd:betann} 有)→ 权重函数参数。省略时程序使用合理默认值。

{phang}
{cmd:iterate(}{it:#}{cmd:)} 设置 Levenberg-Marquardt 算法最大迭代次数,默认 200。

{phang}
{cmd:tolerance(}{it:#}{cmd:)} 设置残差平方和的相对收敛容忍度,默认 1e-7。

{phang}
{cmd:nolog} 不显示迭代日志。

{dlgtab:标准误/推断}

{phang}
{cmd:robust} 返回基于 Gauss-Newton 雅可比矩阵计算的异方差稳健(HC0)三明治标准误。


{title:数据布局}

{pstd}
{cmd:midasreg} 要求数据以{it:高频}为单位组织——每行一个高频期(例如一个月)。低频变量({it:被解释变量}与{it:低频自变量})只在每个低频期的最后一个高频行上取非缺失值,其余高频行应为缺失。数据须先用 {cmd:tsset} 声明。

{pstd}
举例:用月度 {it:hfvar} 预测季度 {it:被解释变量} 时,每个季度只在第三个月记录 y 的取值,前两个月为缺失;同时配合 {cmd:lags(0 11)} 与 {cmd:mratio(3)},程序就会自动使用 {it:hfvar} 的 12 个月度观测(覆盖最近 4 个季度)。


{title:e() 中保存的结果}

{pstd}
{cmd:midasreg} 把以下结果存入 {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: 标量}{p_end}
{synopt:{cmd:e(N)}}观测数{p_end}
{synopt:{cmd:e(k)}}估计参数个数{p_end}
{synopt:{cmd:e(df_r)}}残差自由度{p_end}
{synopt:{cmd:e(rss)}}残差平方和{p_end}
{synopt:{cmd:e(rmse)}}均方根误差{p_end}
{synopt:{cmd:e(sigma2)}}估计的残差方差{p_end}
{synopt:{cmd:e(mratio)}}频率比 m{p_end}
{synopt:{cmd:e(K)}}高频滞后总数{p_end}
{synopt:{cmd:e(kmin)}、{cmd:e(kmax)}}高频滞后阶范围{p_end}
{synopt:{cmd:e(nwpar)}}权重函数参数个数{p_end}
{synopt:{cmd:e(nslope)}}是否单独估计 MIDAS 主系数(1=是,0=否){p_end}
{synopt:{cmd:e(nylags)}}被解释变量低频滞后个数{p_end}
{synopt:{cmd:e(hascons)}}是否含常数项(1=是,0=否){p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: 宏}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:midasreg}{p_end}
{synopt:{cmd:e(depvar)}}被解释变量名{p_end}
{synopt:{cmd:e(hfvar)}}高频变量名{p_end}
{synopt:{cmd:e(weight)}}所用权重函数{p_end}
{synopt:{cmd:e(wlabel)}}权重函数的可读标签{p_end}
{synopt:{cmd:e(lfvars)}}低频回归量列表{p_end}
{synopt:{cmd:e(ylags)}}AR 滞后阶列表{p_end}
{synopt:{cmd:e(vcetype)}}协方差类型("OLS"或"Robust"){p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: 矩阵}{p_end}
{synopt:{cmd:e(b)}}系数向量{p_end}
{synopt:{cmd:e(V)}}协方差矩阵{p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: 函数}{p_end}
{synopt:{cmd:e(sample)}}标识估计样本{p_end}


{title:示例}

{pstd}用滞后的月度就业增长预测季度 GDP 增长,使用归一化指数 Almon 权重并加入 1 阶季度 AR 项:{p_end}
{phang2}{cmd:. tsset month_id}{p_end}
{phang2}{cmd:. midasreg gdp_growth, hfvar(emp_growth) lags(0 11) mratio(3) ylags(1)}{p_end}

{pstd}同样的模型,改用 Ghysels-Santa-Clara-Valkanov Beta 权重并报告稳健标准误:{p_end}
{phang2}{cmd:. midasreg gdp_growth, hfvar(emp_growth) lags(0 11) mratio(3) ///}{p_end}
{phang2}{cmd:    wscheme(beta) ylags(1) robust}{p_end}

{pstd}用 22 个日度收益预测月度已实现波动率(类 HAR-RV 设定):{p_end}
{phang2}{cmd:. midasreg rv_m, hfvar(rv_d) lags(1 22) mratio(22) wscheme(nealmon)}{p_end}

{pstd}无约束 MIDAS——每个高频滞后单独一个系数,等价于对全部滞后做 OLS;K 较小时推荐:{p_end}
{phang2}{cmd:. midasreg gdp_growth, hfvar(emp_growth) lags(0 5) mratio(3) wscheme(umidas)}{p_end}

{pstd}用户自定义初始值:{p_end}
{phang2}{cmd:. midasreg y, hfvar(x) lags(0 11) mratio(3) wscheme(nealmon) ///}{p_end}
{phang2}{cmd:    initial(0  1  0.5 -0.05)}{p_end}


{title:说明与注意事项}

{pstd}
{it:关于识别。} 使用 {cmd:wscheme(nealmon)}、{cmd:wscheme(beta)} 或 {cmd:wscheme(betann)} 时,权重被归一化到和为 1,因此高频系数的"总和"由一个单独的主系数 {it:lambda}(在输出中记为 {it:hfvar}_slope)识别;此时权重函数参数仅决定滞后剖面的{it:形状}。而 {cmd:wscheme(almon)} 与 {cmd:wscheme(umidas)} 直接识别各滞后系数本身,不再单独估计 lambda。

{pstd}
{it:关于滞后阶数选择。} 由于权重函数把所有高频滞后系数绑在一起,K = kmax-kmin+1 可以选得比较大而不必担心自由度损失。Foroni、Marcellino、Schumacher(2015)指出,当 K 较小时无约束 MIDAS 反而常常表现更好;此时可用 {cmd:wscheme(umidas)}。

{pstd}
{it:关于初始值。} 归一化指数 Almon 或 Beta 权重的非线性最小二乘对初始值比较敏感——尤其当默认初始下 {it:lambda} 的符号方向不对时。如果不收敛或迭代日志显示 lambda 发散,请通过 {cmd:initial()} 提供更合理的起始值,或先用 {cmd:wscheme(umidas)} 跑一遍观察滞后剖面再选择合适的参数权重。


{title:已实现的权重函数一览}

{pstd}
{hline 70}{p_end}
{p 4 4 2}{bf:选项}       {bf:含义}                                              {bf:参数个数}{p_end}
{hline 70}{p_end}
{p 4 4 2}nealmon       归一化指数 Almon 多项式(默认)                            2{p_end}
{p 4 4 2}beta          归一化 Beta(末尾权重为 0)                                 2{p_end}
{p 4 4 2}betann        归一化 Beta(末尾权重非零)                                 3{p_end}
{p 4 4 2}almon         二阶非归一化 Almon 多项式                                  3{p_end}
{p 4 4 2}umidas        无约束 MIDAS(每个滞后一个系数)                            K{p_end}
{hline 70}{p_end}


{title:安装方法}

{pstd}
把 {cmd:midasreg.ado} 和 {cmd:midasreg.sthlp} 两个文件放进 Stata 的个人 ado 目录。在 Stata 命令窗口运行 {cmd:sysdir} 可查看 PERSONAL 行所指目录,常见位置如下:

{p 8 8 2}macOS/Linux: {bf:~/ado/personal/}{p_end}
{p 8 8 2}Windows: {bf:C:\ado\personal\}{p_end}

{pstd}
Stata 按命令首字母在子目录中查找 ado 文件,所以请把两个文件放进 {bf:personal} 下的 {bf:m/} 子目录(如 {bf:~/ado/personal/m/})。然后重启 Stata,或运行 {cmd:discard} 命令清空已加载程序缓存,即可使用。

{pstd}
安装好后输入 {cmd:help midasreg} 可随时查看本帮助文档。


{title:参考文献}

{phang}
Andreou, E., E. Ghysels, and A. Kourtellos. 2010. Regression models with mixed sampling frequencies. {it:Journal of Econometrics} 158: 246-261.

{phang}
Foroni, C., M. Marcellino, and C. Schumacher. 2015. Unrestricted mixed data sampling (MIDAS): MIDAS regressions with unrestricted lag polynomials. {it:Journal of the Royal Statistical Society A} 178: 57-82.

{phang}
Ghysels, E., V. Kvedaras, and V. Zemlys. 2016. Mixed frequency data sampling regression models: The R package midasr. {it:Journal of Statistical Software} 72(4): 1-35.

{phang}
Ghysels, E., P. Santa-Clara, and R. Valkanov. 2002. The MIDAS touch: Mixed data sampling regression models. Working paper, UNC and UCLA.

{phang}
Ghysels, E., P. Santa-Clara, and R. Valkanov. 2006. Predicting volatility: Getting the most out of return data sampled at different frequencies. {it:Journal of Econometrics} 131: 59-95.


{title:作者说明}

{pstd}
本 Stata 实现遵循 Ghysels、Kvedaras、Zemlys(2016)给出的 MIDAS 估计框架,欢迎反馈意见与 bug。


{title:另请参阅}

{psee}
帮助:{help regress}、{help nl}、{help arima}、{help tsset}
