# CSU-Beamer 模板审计报告

- 审计日期: 2026-09-01
- 模板仓库: https://github.com/Sakulyn/CSU-Beamer（Sapienza "sintef" 主题 fork，原作者 Federico Zenith / SINTEF）
- 克隆位置: `talk/beamer_template_test/CSU-Beamer/`
- 提交: `f97f5235a04e9f807ee6aa5b7bd1136463d090f2`（shallow, depth 1）
- 涉及文件: `csu-beamer.tex`、`beamerthemesintef.sty`（321 行）、`sintefcolor.sty`（28 行）
- 编译环境: XeTeX 3.141592653-2.6-0.999998（TeX Live 2026），字体未安装新包（保持系统现状）

---

## 一、编译与验证记录

| 阶段 | 命令 | 结果 |
|---|---|---|
| 第 1 次（原样） | `xelatex -interaction=nonstopmode csu-beamer.tex` | 退出码 1；**PDF 已产出 22 页**，但有缺字体错误 |
| 第 2 次（原样） | 同上 | 退出码 1；同上 |
| 修复后第 1 次 | 同上（改后文件） | 退出码 0 |
| 修复后第 2 次 | 同上（改后文件） | 退出码 0，日志 0 错误、0 缺字 |

### 原样编译的失败原因（缺字体，属于允许修复类）

原样编译并非真正"失败"（PDF 已产出），但 exit=1 且标题页文字缺失，根因是两处 macOS/Windows 专属字体在本 Linux 机器上不存在：

1. `csu-beamer.tex` L13 `\newfontfamily\menlo{Menlo}` — **Menlo**（macOS 自带等宽字体）不存在 → fontspec 报错。`\menlo` 实际只用于 listings 行号样式（`numberstyle=\tiny\menlo`，且 `numbers=left` 被注释，行号未启用）。
2. `csu-beamer.tex` L46 `\title{\CJKfontspec{Microsoft YaHei}报告标题}` — **Microsoft YaHei**（Windows 自带）不存在 → 标题 4 个汉字回退到 Times New Roman（无 CJK 字形）→ 日志 92 处 "Missing character"（标题页 + 每页脚注里的标题都缺字）。

### 最小修复（仅缺字体，未动任何视觉设计）

| 位置 | 原内容 | 改为 | 理由 |
|---|---|---|---|
| L13 | `\newfontfamily\menlo{Menlo}` | `\newfontfamily\menlo{Hack}` | Hack 为本机已装等宽字体，仅作行号字体占位 |
| L46 | `\CJKfontspec{Microsoft YaHei}` | `\CJKfontspec{LXGWWenKai-Regular.ttf}` | 与本模板 L15 已使用的 CJK 无衬线字体一致（LXGWWenKai） |

修复后两次编译 exit=0。最终 `csu-beamer.pdf`：**22 页**，页面尺寸 **453.54 × 255.12 pt（16:9）**。
内嵌字体（pdffonts）：Times New Roman（Regular/Bold/Italic）、LXGWWenKai-Regular、LXGWZhenKaiGB（CJK 粗体）、CMR/CMSY/CMEX（数学）、LMMono、Caladea。

### 页码映射（渲染用）

| 页 | 内容 | 类型 |
|---|---|---|
| 1 | 标题页（报告标题 / 2025年11月） | 标题页 |
| 2 / 8 / 20 | 目录（Introduction / Personalization / Summary） | 章节页（彩色 TOC） |
| 3 | Beamer vs. PowerPoint | 普通内容页 |
| 4 | Title page（命令讲解，含 block） | 内容页 |
| 5–7 | Writing a Simple Slide（1 个 frame 拆 3 页，含代码） | 内容页 |
| 10 | Blocks | block 示例 |
| 14 | Splitting in Columns | 分栏示例 |
| 16 | Side-Picture Slides | 侧图页（sidepic） |
| 22 | 报告标题 + "Thank you for listening!" | 结尾页（backmatter） |

逻辑 frame 数为 19（脚注编号 X/19），22 页 = 19 frame + 3 个 frame 拆页续页
（L320–321 `\setbeamertemplate{frametitle continuation}{}` 隐藏续页标记）。

---

## 二、十个审计问题

### 1. 标题页如何实现？
`beamerthemesintef.sty` 两处配合：
- **标题页模板**（L153–171）：`\vskip0pt plus 1filll` 推到下方 + `\hspace{-12mm}\vspace{20mm}` 微调位置，`beamercolorbox[wd=0.72\textwidth,sep=10pt,leftskip=8mm]`，内容依次为 `\inserttitle`（Huge 粗体、maincolor 色，L107）、`\insertsubtitle`、`\@courseLabel`、`\insertauthor`（+ `\@IDnumber`）、`\insertdate`。
- **\maketitle 重写**（L197–212）：若设置了 `\titlebackground`，先注入背景模板（星号形式 → `\SplitBackground` → 半图 `\TikzSplitSlide`；非星号 → 整页图 `\includegraphics[height=\paperheight]`），再 `\begin{frame}\titlepage\end{frame}`。
- 本 deck：L19 `\titlebackground*{assets/background-railway}`，L20 `\boolfalse{splittitle}` → 标题页为**整页背景图** + 左下方文字块。

### 2. 普通 frame 标题如何实现？
L136–143：`\setbeamertemplate{frametitle}` — `\vspace*{-5.5ex}` 上移，`beamercolorbox[leftskip=2cm]`（给左上 logo 让位），第一行粗体 `\insertframetitle`，第二行 `\insertframesubtitle`。
L316–318：`\pretocmd\beamer@checkframetitle` 自动给每帧加副标题 `\thesection \, \secname`（如 "2 Personalization"）。

### 3. logo 放在哪里？
- L73–75：`\pgfdeclareimage[width=0.09\paperwidth]{logo}{assets/logo_RGB}`（+ 反色版 `whitelogo` = `logo_RGB_negative`）。
- L133：headline 模板 = `\hspace{0.06\textwidth}\pgfuseimage{\@logo}` — **每页左上角**，宽为纸宽 9%。
- `\themecolor{main}` 时自动切换为反色版（L86, L93）。

### 4. 页脚 / 页码如何实现？
L44–55 footline 模板：`beamercolorbox[wd=\textwidth,ht=5mm,dp=3mm,rightskip=1cm,leftskip=1cm]`，左 `\insertframenumber/\inserttotalframenumber`，右 `\hfill` + `\@footlinepayoff`（默认 = `\insertauthor \enspace$\vert$\enspace \inserttitle`，L44–46；本 deck 未设 author → 实际显示 `| 报告标题`）。
L57–71 `\footlinecolor[1]`：参数为空 → 脚注透明（darkgray 字）+ block title 变 maincolor；参数非空 → 脚注白字、底色=参数色，block title 同步同色。
本 deck 的切换：默认空（p1–7）→ L62 `\footlinecolor{maincolor}`（p8 起）→ L67 `\footlinecolor{sintefred}`（p14 起）→ L71 `\footlinecolor{}`（p17 起，含结尾页）。

### 5. 章节页如何实现？
**chapter 环境**（L222–249）：`\themecolor{main}` 整页 maincolor 底色（或自定义色），清空脚注，可选 `\TikzSplitSlide` 右半图，标题 `\vspace*{8ex}`、宽 0.45\textwidth，文字在 0.35\textwidth minipage。
**注意**：示例 deck **未调用 chapter 环境**（sections/Personalization/6-chapter-slides.tex 只是讲解它）。实际出现的彩色整页是 **`\AtBeginSection` 目录页**（L294–302）：`\themecolor{main}` + `\begin{frame}{目录}\tableofcontents[hideallsubsections]`；当前 section 黄色 `\blacktriangleright` 高亮，其余 `white!80!gray!50`（L305–313）。共 3 页（p2/p8/p20）。

### 6. block 圆角与配色？
- L115：`\setbeamertemplate{blocks}[rounded]` — **圆角**。
- 配色由 `\themecolor` / `\footlinecolor` 联动：
  - 白色主题（默认）：block title 白字/maincolor 底（L94），block body darkgray/sintefgrey（L95）。
  - main 主题：block title maincolor 字/sintefgrey 底（L87–88）。
  - `\footlinecolor{}`：block title 白字/maincolor 底（L62）。
- `colorblock` 环境（L122–130）可自定义整块颜色（demo 的 3-using-colors 用了 `\testcolor`/`\colorbox` 而非 colorblock）。

### 7. 图文分屏页如何实现？
**sidepic 环境**（L252–272，被 7-side-picture-slides.tex 实际使用）：背景 = `\hspace*{0.6\paperwidth}\includegraphics[height=\paperheight]{#1}`（图占右 40% 全高）；frame 标题 `rightskip=0.4\textwidth`；正文 minipage 宽 0.6\textwidth。
**\TikzSplitSlide**（L174–185，标题页/章节页用）：`\rule{0.56\paperwidth}{0pt}` 占位 + tikzpicture 多边形 clip（纸宽 0.1→0.5 的斜边），`\includegraphics[height=\paperheight]` 右半图。

### 8. 16:9 在哪里定义？
`beamerthemesintef.sty` L31–33：`\RequirePackage{geometry}\geometry{paperwidth=16cm,paperheight=9cm}`（硬编码 16:9，注释 "Force 16:9 aspect ratio"）。实测 PDF 页 453.54 × 255.12 pt，符合。

### 9. 中文字体如何定义？
`csu-beamer.tex`：
- L11–12：`\usepackage{ctex}` + `\usepackage{fontspec}`（XeLaTeX 引擎）。
- L14：`\setsansfont{Times New Roman}` — 西文无衬线。
- L15：`\setCJKsansfont{LXGWWenKai-Regular.ttf}[BoldFont = LXGWZhenKaiGB-Regular.ttf]` — CJK 用**仓库自带字体文件**（无需系统安装），粗体用 ZhenKaiGB。
- L46（已修）：标题 `\CJKfontspec{...}` 指定标题 CJK 字体（原 Microsoft YaHei → 改为 LXGWWenKai-Regular.ttf）。
- 主题自带 caladea/carlito（Cambria/Calibri 克隆，L35–37）作为西文衬线/无衬线回退。

### 10. 数学字体如何定义？
`csu-beamer.tex` L17：`\usefonttheme[onlymath]{serif}` — **仅数学用衬线（Computer Modern）**，正文仍无衬线。pdffonts 确认 CMR10/CMR8/CMSY/CMEX10/yfrak 已嵌入。数学公式原生 LaTeX 排版。

---

## 三、后续替换候选（本轮只记录，不修改）

| 候选 | 内容 | 说明 |
|---|---|---|
| **A** | `assets/background-railway.jpg` 标题页背景 | 现为铁路照片整页背景（splittitle=false）；可换为我们的实验/示意图背景 |
| **B** | `logo_RGB.png` / `logo_RGB_negative.png` | 左上角每页 logo（宽 0.09 纸宽）；替换为我们的实验室/组 logo，需同时提供反色版（maincolor 页面用） |
| **C** | `sintefcolor.sty` maincolor = RGB **8,114,175** / CMYK 85,52,16,0 | 主题主蓝；需与当前 PPT 主蓝对比后决定是否换色 |
| **D** | 字体组合：Times New Roman（西文）+ LXGWWenKai（中文） | 需与当前 PPT 的微软雅黑/Microsoft YaHei 对比观感 |
| **E** | 脚注（页脚）| 现为 `页码/总页数` + `作者 | 标题`（本 deck 无作者 → `| 报告标题`）；候选简化为仅页码或"章节+页码" |
| **F** | 左上角 headline logo | 候选：保持 / 缩小 / 仅章节页显示 |
| **G** | 章节页 | 现为彩色整页"目录"（AtBeginSection TOC），主题另有 chapter 环境（未被 demo 使用）；很可能用于替换当前 PPT 的"章节高亮"导航页 —— **本轮不改** |

## 四、值得保留的机制（keep 清单）

1. **16:9 硬编码**（geometry 16cm×9cm）— 与现 PPT 比例一致
2. **section/frame 结构** + `\AtBeginSection` 自动目录页 + 自动副标题（`\thesection \secname`）
3. **无导航符号**（`\setbeamertemplate{navigation symbols}{}`）
4. **block 环境**（圆角、颜色与脚注联动）— 结构清晰、可自定义色块（colorblock）
5. **columns 分栏**（原生 beamer，demo 5-splitting-in-columns）
6. **TikZ 支持**（backgrounds/arrows/shapes/tikzmark/calc 已加载）— 画图能力强
7. **分屏图文助手**：sidepic 环境 + \TikzSplitSlide（标题页/章节页半图）
8. **章节页机制**：chapter 环境（整页色 + 可选半图），虽 demo 未用但机制完整
9. **原生 LaTeX 数学**（`\usefonttheme[onlymath]{serif}`，CM 字体嵌入）
10. **PDF-first 工作流**：xelatex 单文件编译即出最终 PDF，字体/公式/版本无兼容性问题（demo 文案原话）
11. **CJK 字体随仓库分发**（LXGWWenKai-Regular.ttf / LXGWZhenKaiGB-Regular.ttf 在 repo 内）— 不依赖系统字体

## 五、渲染产物

- 输出目录: `talk/beamer_template_test/render/`
- 全 22 页: `page-01.png` … `page-22.png`（130 dpi）
- 关键页命名副本: `title-page.png`(p1)、`normal-content.png`(p3)、`blocks.png`(p10)、`columns.png`(p14)、`chapter-slide.png`(p2 目录页)、`side-picture.png`(p16)、`ending.png`(p22)
- 联系表: `contact-sheet.png`（5×5 拼图，4125×2335）
- 原 PDF 未做任何修改之外的改动；渲染仅基于最终 `csu-beamer.pdf`。
