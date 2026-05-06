{smcl}
{* *! version 1.0  06may2026}{...}
{vieweralsosee "irf" "help irf"}{...}
{vieweralsosee "irf graph" "help irf_graph"}{...}
{vieweralsosee "var" "help var"}{...}
{vieweralsosee "svar" "help svar"}{...}
{title:Title}

{phang}
{bf:irf_dualband} {hline 2} 从 IRF 文件提取数据并绘制双置信带脉冲响应图


{title:Syntax}

{p 8 17 2}
{cmd:irf_dualband} {cmd:,}
{cmdab:irff:ile(}{it:filename}{cmd:)}
{cmdab:imp:ulse(}{it:varname}{cmd:)}
{cmdab:r:esponse(}{it:varname}{cmd:)}
{cmdab:irfn:ame(}{it:string}{cmd:)}
[{it:options}]


{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :必选项}
{synopt:{opt irff:ile(filename)}}IRF 文件路径(.irf 文件){p_end}
{synopt:{opt imp:ulse(varname)}}冲击变量名{p_end}
{synopt:{opt r:esponse(varname)}}响应变量名{p_end}
{synopt:{opt irfn:ame(string)}}IRF 集合名(由 irf create 时指定){p_end}

{syntab :IRF 类型与置信带}
{synopt:{opt t:ype(string)}}响应类型: {bf:sirf} (默认)、{bf:oirf}、{bf:irf}{p_end}
{synopt:{opt i:nner(#)}}内置信带水平 (%),默认 68{p_end}
{synopt:{opt o:uter(#)}}外置信带水平 (%),默认 90{p_end}
{synopt:{opt norm:alize}}用 impulse 对自身 step 0 响应做基准 rescale{p_end}

{syntab :图形外观}
{synopt:{opt t:itle(string)}}图标题{p_end}
{synopt:{opt y:title(string)}}纵轴标签{p_end}
{synopt:{opt x:title(string)}}横轴标签{p_end}
{synopt:{opt innerc:olor(string)}}内带颜色,默认 RGB "139 0 0" (深红){p_end}
{synopt:{opt outerc:olor(string)}}外带颜色,默认 RGB "219 149 168" (浅粉){p_end}
{synopt:{opt innera:lpha(#)}}内带透明度 (0-100),默认 60{p_end}
{synopt:{opt outera:lpha(#)}}外带透明度 (0-100),默认 40{p_end}
{synopt:{opt n:ame(string)}}图形名,默认 "{it:impulse}_{it:response}",自动 replace{p_end}
{synopt:{opt e:xport(filename)}}导出图片(支持 .png/.pdf/.eps 等){p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:irf_dualband} 从 Stata 的 IRF 文件 (.irf) 中提取指定 impulse-response 对的
脉冲响应数据,使用渐近正态置信区间(基于 std{it:type}),并在同一张图上叠加
两个不同水平的置信带,生成常见于宏观经济学论文的"双置信带 IRF 图"。

{pstd}
程序使用 IRF 文件中存储的标准误 (std{it:type}) 现场计算两个不同水平的置信
区间,因此只需要一次 {cmd:irf create} 即可生成多置信水平的图。


{title:Options}

{dlgtab:必选项}

{phang}
{opt irffile(filename)} IRF 文件路径,通常以 .irf 结尾。该文件由 {cmd:irf create}
的 {opt set()} 选项指定生成。

{phang}
{opt impulse(varname)} 结构冲击来源变量,需与 IRF 文件中的 {cmd:impulse} 列匹配。

{phang}
{opt response(varname)} 响应变量,需与 IRF 文件中的 {cmd:response} 列匹配。

{phang}
{opt irfname(string)} IRF 集合名,即 {cmd:irf create {it:name}} 时的名字。
可用 {cmd:irf describe} 查看可用名称。

{dlgtab:IRF 类型与置信带}

{phang}
{opt type(string)} 指定使用的响应类型:

{p 12 12 2}
{bf:sirf} - 结构脉冲响应(SVAR 默认),需要 IRF 文件包含此列{p_end}
{p 12 12 2}
{bf:oirf} - 正交化脉冲响应(Cholesky 分解){p_end}
{p 12 12 2}
{bf:irf}  - 简化式脉冲响应{p_end}

{phang}
{opt inner(#)} 内层(较窄)置信带的置信水平,百分数形式,默认 68。

{phang}
{opt outer(#)} 外层(较宽)置信带的置信水平,百分数形式,默认 90。

{phang}
{opt normalize} 启用归一化。开启后,所有响应除以
"impulse → impulse at step 0" 的值,使冲击大小被解释为
"让 impulse 变量在当期上升 1 个单位"。不启用时显示原始一标准差冲击的响应。

{dlgtab:图形外观}

{phang}
{opt title(string)} 图标题,默认为 "{it:impulse} → {it:response}"。

{phang}
{opt ytitle(string)} 纵轴标签,默认 "Response"。

{phang}
{opt xtitle(string)} 横轴标签,默认 "Steps"。

{phang}
{opt innercolor(string)}, {opt outercolor(string)} 置信带颜色,
可用 Stata 命名颜色(如 {bf:cranberry})、RGB(如 "139 0 0")。

{phang}
{opt inneralpha(#)}, {opt outeralpha(#)} 颜色不透明度,0 (完全透明) 到 100 (不透明)。
对 RGB 颜色生效。

{phang}
{opt name(string)} 图形窗口名,默认为 "{it:impulse}_{it:response}",
始终自动 {bf:replace}。可用于一次画多张图后再用 {cmd:graph combine} 拼图。

{phang}
{opt export(filename)} 导出图片到文件,根据扩展名自动识别格式
(.png, .pdf, .eps, .tif, .svg 等)。


{title:Remarks}

{pstd}
{ul:置信区间方法}: 程序使用渐近正态近似,即
{it:pe} ± {it:z} × {it:std} 形式的对称区间。如需 bootstrap 置信带,
请在 {cmd:irf create} 时加 {opt bs reps(#)} 选项,但此时 IRF 文件
存储的是分位数端点 (sirf_lo/sirf_up),而非标准误,本程序不直接支持
该格式(可手动改写)。

{pstd}
{ul:IRF 文件结构}: IRF 文件本质是一个特殊数据集。可用以下命令查看内容:

{p 8 8 2}
{cmd:. preserve}{break}
{cmd:. use "myirf.irf", clear}{break}
{cmd:. describe}{break}
{cmd:. restore}

{pstd}
{ul:何时用 sirf vs oirf vs irf}: 一般 SVAR 模型用 sirf,
普通 VAR 用 oirf(Cholesky 识别)。如果 {cmd:irf create} 时模型是 VAR
而非 SVAR,IRF 文件中不会有 {bf:sirf} 列。


{title:Examples}

{pstd}{ul:基本流程}: 估 SVAR、生成 IRF、画双带图。{p_end}

{phang2}{cmd:. var output inflation, lags(1)}{p_end}
{phang2}{cmd:. matrix Aid = (1, 0 \ ., 1)}{p_end}
{phang2}{cmd:. matrix Bid = (., 0 \ 0, .)}{p_end}
{phang2}{cmd:. svar output inflation, lags(1) aeq(Aid) beq(Bid)}{p_end}
{phang2}{cmd:. irf create svar1, set(myirf, replace) step(36)}{p_end}

{pstd}默认设置(68% 内带 / 90% 外带,深红+浅粉配色):{p_end}
{phang2}{cmd:. irf_dualband, irffile("myirf.irf") impulse(output) response(inflation) ///}{p_end}
{phang2}{cmd:>     irfname(svar1)}{p_end}

{pstd}归一化为 "+1 unit output shock" 的响应,导出 PNG:{p_end}
{phang2}{cmd:. irf_dualband, irffile("myirf.irf") impulse(output) response(inflation) ///}{p_end}
{phang2}{cmd:>     irfname(svar1) normalize ///}{p_end}
{phang2}{cmd:>     title("Response of inflation to +1pp output shock") ///}{p_end}
{phang2}{cmd:>     ytitle("Percentage points") ///}{p_end}
{phang2}{cmd:>     xtitle("Months after policy intervention") ///}{p_end}
{phang2}{cmd:>     export("irf_plot.png")}{p_end}

{pstd}自定义置信水平和颜色:{p_end}
{phang2}{cmd:. irf_dualband, irffile("myirf.irf") impulse(money) response(output) ///}{p_end}
{phang2}{cmd:>     irfname(svar1) inner(50) outer(95) ///}{p_end}
{phang2}{cmd:>     innercolor("navy") outercolor("ltblue") ///}{p_end}
{phang2}{cmd:>     inneralpha(70) outeralpha(40)}{p_end}

{pstd}画多张图后拼版:{p_end}
{phang2}{cmd:. irf_dualband, ..., impulse(output) response(output)    name(g1)}{p_end}
{phang2}{cmd:. irf_dualband, ..., impulse(output) response(inflation) name(g2)}{p_end}
{phang2}{cmd:. irf_dualband, ..., impulse(money)  response(output)    name(g3)}{p_end}
{phang2}{cmd:. irf_dualband, ..., impulse(money)  response(inflation) name(g4)}{p_end}
{phang2}{cmd:. graph combine g1 g2 g3 g4, cols(2)}{p_end}


{title:Stored results}

{pstd}
{cmd:irf_dualband} 不返回 r() 或 e() 结果。如开启 {opt normalize},
基准 scale 会显示在控制台。生成的图保存在内存中,名字由 {opt name()} 指定。


{title:Author}

{pstd}
陆震坤，lzkzzer@163.com


{title:Also see}

{psee}
Online: {help irf}, {help irf_create:irf create}, {help irf_graph:irf graph},
{help var}, {help svar}, {help twoway_rarea:twoway rarea}
{p_end}