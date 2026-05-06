{smcl}
{* *! version 1.0  06may2026}{...}
{vieweralsosee "var" "help var"}{...}
{vieweralsosee "svar" "help svar"}{...}
{vieweralsosee "irf" "help irf"}{...}
{title:Title}

{phang}
{bf:sim_svar} {hline 2} 用 Mata 模拟结构向量自回归 (SVAR) 时间序列


{title:Syntax}

{p 8 17 2}
{cmd:sim_svar} {cmd:,}
{cmdab:n:vars(}{it:#}{cmd:)}
{cmdab:l:ags(}{it:#}{cmd:)}
{cmdab:c:oef(}{it:matname}{cmd:)}
[{it:options}]


{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :必选项}
{synopt:{opt n:vars(#)}}变量个数 n{p_end}
{synopt:{opt l:ags(#)}}滞后阶数 p{p_end}
{synopt:{opt c:oef(matname)}}n × (n×p) 简化式滞后系数矩阵 [Phi_1 ... Phi_p]{p_end}

{syntab :可选项}
{synopt:{opt a:mat(matname)}}n × n 结构矩阵 A,默认单位阵{p_end}
{synopt:{opt b:mat(matname)}}n × n 结构矩阵 B,默认单位阵{p_end}
{synopt:{opt const(matname)}}n × 1 常数向量,默认零向量{p_end}
{synopt:{opt n:obs(#)}}保留观测数,默认 500{p_end}
{synopt:{opt b:urnin(#)}}烧入期长度,默认 200{p_end}
{synopt:{opt s:eed(#)}}随机种子,默认 12345{p_end}
{synopt:{opt na:mes(string)}}变量名列表,默认 y1 y2 ...{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:sim_svar} 根据用户给定的结构参数模拟 SVAR(p) 模型的时间序列数据。
结构方程为:

{p 8 8 2}
A y_t = c + Phi_1 y_{t-1} + ... + Phi_p y_{t-p} + B eps_t,    eps_t ~ N(0, I_n)

{pstd}
其中 A、B 是 n × n 结构矩阵,Phi_j 是 n × n 滞后系数矩阵,c 是 n × 1 常数。
程序在 Mata 中完成矩阵递推,速度远快于纯 Stata 循环,并自动检查
模型的平稳性(伴随矩阵最大特征根模)。


{title:Options}

{phang}
{opt nvars(#)} 指定 VAR 系统的变量数 n,必须为正整数。

{phang}
{opt lags(#)} 指定滞后阶数 p,必须为正整数。

{phang}
{opt coef(matname)} 简化式滞后系数矩阵,大小必须为 n × (n×p),
按 [Phi_1 | Phi_2 | ... | Phi_p] 的顺序水平拼接。例如 n=2、p=2 时,
矩阵为 2 × 4。

{phang}
{opt amat(matname)} 结构矩阵 A,n × n。如不指定,默认为单位阵
(此时退化为简化式 VAR)。常用于 Cholesky 短期识别(下三角)。

{phang}
{opt bmat(matname)} 结构矩阵 B,n × n。结构冲击通过 A^{-1} B 映射到
简化式扰动。

{phang}
{opt const(matname)} 常数项向量,n × 1 列向量或 1 × n 行向量(自动转置)。

{phang}
{opt nobs(#)} 烧入期之后保留的观测数,默认 500。

{phang}
{opt burnin(#)} 烧入期长度,默认 200。从零初值开始递推 burnin 期,
再保留之后的 nobs 期,以减少初值影响。如果系统接近不平稳,
建议加大此值。

{phang}
{opt seed(#)} 随机种子,默认 12345,用于结构冲击 eps_t 的抽样。

{phang}
{opt names(string)} 变量名列表,空格分隔。个数必须等于 nvars。
如不指定,默认 y1 y2 ... y{it:n}。


{title:Remarks}

{pstd}
{ul:平稳性检查}: 程序内部构造伴随矩阵 F,计算其特征根。如果最大模 ≥ 1,
显示警告但仍生成数据(可能不平稳)。

{pstd}
{ul:简化式协方差}: 简化式扰动 u_t = A^{-1} B eps_t 的协方差为
A^{-1} B B' (A^{-1})'。如果只想指定简化式协方差 Sigma 而不做结构识别,
可设 amat 为单位阵、bmat 为 cholesky(Sigma)。

{pstd}
{ul:数据集}: 程序会清空当前数据集,生成时间变量 t 并 tsset。
原始数据丢失前请先 preserve 或保存。


{title:Examples}

{pstd}两变量 SVAR(1),Cholesky 短期识别:{p_end}
{phang2}{cmd:. matrix A   = (1, 0 \ 0.5, 1)}{p_end}
{phang2}{cmd:. matrix B   = I(2)}{p_end}
{phang2}{cmd:. matrix Phi = (0.6, 0.1 \ 0.2, 0.7)}{p_end}
{phang2}{cmd:. matrix c   = (0.0 \ 0.0)}{p_end}
{phang2}{cmd:. sim_svar, nvars(2) lags(1) coef(Phi) amat(A) bmat(B) const(c) ///}{p_end}
{phang2}{cmd:>     nobs(1000) burnin(200) seed(42) names("output inflation")}{p_end}

{pstd}三变量 SVAR(2):{p_end}
{phang2}{cmd:. matrix Phi = ( 0.5, 0.1, 0.0,  0.1, 0.0, 0.0 \ ///}{p_end}
{phang2}{cmd:>               0.0, 0.4, 0.1,  0.0, 0.1, 0.0 \ ///}{p_end}
{phang2}{cmd:>               0.1, 0.0, 0.6,  0.0, 0.0, 0.1 )}{p_end}
{phang2}{cmd:. matrix A = (1, 0, 0 \ 0.3, 1, 0 \ 0.2, 0.4, 1)}{p_end}
{phang2}{cmd:. sim_svar, nvars(3) lags(2) coef(Phi) amat(A) nobs(800) seed(7)}{p_end}

{pstd}模拟后估计验证:{p_end}
{phang2}{cmd:. var output inflation, lags(1)}{p_end}
{phang2}{cmd:. matrix Aid = (1, 0 \ ., 1)}{p_end}
{phang2}{cmd:. matrix Bid = (., 0 \ 0, .)}{p_end}
{phang2}{cmd:. svar output inflation, lags(1) aeq(Aid) beq(Bid)}{p_end}


{title:Stored results}

{pstd}
{cmd:sim_svar} 不返回 r() 或 e() 结果。模拟数据存储在当前数据集中,
含变量 t (时间) 和 names() 指定的变量。


{title:Author}

{pstd}
陆震坤, lzkzzer@163.com


{title:Also see}

{psee}
Online: {help var}, {help svar}, {help irf}, {help dsge}
{p_end}