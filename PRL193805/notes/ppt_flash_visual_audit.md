# PPT 快速视觉审计 (2026-08-31)

审计对象（唯一真源，WPS 手工修改版）：
`ppt/exports/260414_郭晋康_组会_navigation_final.pptx`
- 路径：`/mnt/sync/My_project/03_following/PRL193805/ppt/exports/260414_郭晋康_组会_navigation_final.pptx`
- mtime：2026-08-31 18:14:25（晚于先前导出 15:44，确为 WPS 手工修改版）
- 大小：1,786,351 字节（先前导出为 2,950,268 字节）
- 页数：32

审阅依据：
- 主要：WPS Presentation GUI 实际显示（DISPLAY=:0，已用 `wpp` 打开确认可正常加载）
- 逐页：LibreOffice 渲染 PNG（几何由 OOXML 决定；字体/折行差异单独标注并对照原生 XML 判定）
- 科研图完整性：以嵌入媒体字节/哈希 + 原生 frame aspect 校验

标注规则：OK / MINOR / NEED_FIX；只记录真实视觉问题。

---

| 页 | 状态 | 问题 | 处理 |
|---|---|---|---|
| 1 | OK | 封面元素完整、居中 | — |
| 2 | OK | 四卡片等宽等距，主线栏居中，I–IV 对齐一致 | — |
| 3 | OK | 章节导航 I，高亮/灰度正确，与 8/15/21 几何一致 | — |
| 4 | OK | (a) 图完整（panel label/ηp/ηp′ 可辨），(c) 图完整 | — |
| 5 | OK | 公式与右侧 UP/LP 图结构完整 | — |
| 6 | NEED_FIX | "振子强度"为孤儿标签（独立文本框，挂在 Imχ 下方，无明确归属；按约定本页只需 P(ω)=χ(ω)E(ω)、极化、χ、Reχ、Imχ） | 删除该标签 |
| 7 | NEED_FIX | ρeg 方程与 P(ω) 之间有一根 40pt 无意义橙色竖线（原生 connector，无语义） | 删除该 connector |
| 8 | OK | 章节导航 II，同 3/15/21 | — |
| 9 | OK | 三卡片+结论栏平衡 | — |
| 10 | OK | 四框链清晰，"热平衡"措辞正确 | — |
| 11 | OK | 分母三项图例对位 | — |
| 12 | OK | T/A/R 图完整（轴、图例、标签） | — |
| 13 | OK | 两科研图完整（axes/legend 齐全） | — |
| 14 | OK | 四条件字重基线统一；序号徽章蓝/灰与卡片底色同节奏（判定为设计条纹，非缺陷）| — |
| 15 | OK | 章节导航 III | — |
| 16 | OK | 两卡+箭头布局正确 | — |
| 17 | NEED_FIX | 返回箭头为从 ρ(t)→P(t) 框左缘斜穿的橙色虚线（穿左卡间距、箭头短桩悬于 α 框左上角）；"P(t) 反馈至 α(t) 方程"标签悬在左卡 D[ρ] 旁，未形成真正闭环 | 改为右侧正交回环线（ρ,P 框右缘 → 右 → 上 → α 框右缘），标签移至回环右下角 |
| 18 | OK | ρ 矩阵足够大，三卡对齐，κ/γ/γφ 小字可读 | — |
| 19 | OK | 公式居中（analytic bridge 留白可接受） | — |
| 20 | OK | 7.51×10⁻⁹ takeaway 清晰，左右平衡 | — |
| 21 | OK | 章节导航 IV | — |
| 22 | OK | 图注与图无真实冲突：原生 bodyPr tIns=69.66pt（anchor=t），PowerPoint/WPS 中文字渲染在图下方（≈486pt 基线）；LibreOffice 忽略该 inset 产生"压图"假象（渲染基线≈435pt）；medida image6.png aspect=1.3333=frame，无拉伸 | 不修改（第三方渲染器差异） |
| 23 | NEED_FIX | 左卡公式+注释较满，右卡标题下大片空白、答案行贴卡底 → 左右失衡 | 右卡答案行上移至卡片中部（仅移动文字垂直位置） |
| 24 | OK | 嵌入媒体哈希 == 原论文 Fig.2 裁剪（sha c527cf7bec07544a），橙/黑路径为原图风格，完整无裁剪 | — |
| 25 | OK | α^(2)(1) 与 ΔT(ω) 公式、卡片基线、底部说明均正常 | — |
| 26 | OK | 示意热图轴/colorbar/示意字样齐全，左右平衡 | — |
| 27 | OK | 完整 Fig.3（a,b）+ colorbar/panel labels 齐全 | — |
| 28 | OK | 完整 Fig.3（c,d），框线强调 (c,d) | — |
| 29 | OK | 公式/三步流程/两解释框无穿字；ket/bra 小字可读 | — |
| 30 | OK | 双热图等尺寸、colorbars 完整、DQC 说明清楚 | — |
| 31 | OK | 三结论可读、两论文框对齐、schematic 不抢焦点 | — |
| 32 | OK | 结束页：原生 XML 五文本框无几何重叠（独立 sp，y 间隔 31–36pt）；"学号：252201032"20pt 单行 ~150pt < 框宽 167.8pt，PowerPoint/WPS 单行不折行；LibreOffice 折行致"数字压导师行"为渲染伪影 | 不修改（第三方渲染器差异） |

统计（修改前）：OK 28 页，MINOR 0 页，NEED_FIX 4 页（6、7、17、23）。

## 修改结果（2026-08-31 18:37 完成）

| 页 | 修改 | 复检 |
|---|---|---|
| 6 | 删除孤儿标签"振子强度"（shape-15 文本框） | ✅ 已消失，Reχ/Imχ 基线对齐不变 |
| 7 | 删除无意义橙色竖线（shape-7 connector） | ✅ 已消失，公式链完整 |
| 17 | 删除斜穿虚线+悬空箭头桩；新增右侧正交回环（ρ,P 框右缘→右→上→α 框右缘，含箭头）；标签移至回环右下角（x=1200, 白底 3 行） | ✅ 闭环清晰、无穿线、不挡方程；第一视觉结论"分子极化反过来驱动腔场"成立 |
| 23 | 右卡答案行上移（frame y 520.32→405.32，text y 570.04→455.04），居中于卡片中部 | ✅ 左右视觉平衡恢复 |

修改后统计：OK 32 / MINOR 0 / NEED_FIX 0。
第三方渲染器差异（未修改）：第 22 页图注（PowerPoint/WPS 依 tIns=69.66pt 渲染在图下方）、第 32 页结束页"学号"行（WPS/PPT 20pt 单行不折行）——均为 LibreOffice 字体/autofit 伪影，WPS 与 PowerPoint 正常，按规则未修改。

最终交付：
- 真源（WPS 手工版，未动）：`ppt/exports/260414_郭晋康_组会_navigation_final.pptx`
- 修改版：`ppt/exports/260414_郭晋康_组会_navigation_flash_review.pptx`（32 页）
- PDF / contact sheet：`ppt/exports/260414_郭晋康_组会_navigation_flash_review.pdf` / `*_contact_sheet.png`
- round-trip 工作区：`molecular_polariton_flash_review_20260831`（全新 import，未继承旧 SVG）

