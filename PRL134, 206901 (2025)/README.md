# Mechanism of Picocavity-Enhanced Raman Spectroscopy: H₂/Ag(111)

**基于第一性原理与有限元仿真的 H₂/Ag(111) 皮腔增强拉曼光谱机制研究**

![Status](https://img.shields.io/badge/Status-Completed-success)
![DFT](https://img.shields.io/badge/DFT-Quantum_ESPRESSO-blue)
![Lang](https://img.shields.io/badge/Scripts-Python_3-green)

## 📖 简介 (Introduction)

本项目旨在复现并解析发表于 **Physical Review Letters** 上的工作 [*Phys. Rev. Lett. 134, 206901 (2025)*]。

针对实验中观测到的 H₂ 与 D₂ 分子在 Ag(111) 表面截然不同的光谱响应（反常同位素效应），本项目通过 **DFT (Quantum ESPRESSO)** 计算了吸附势能面与电子结构，并结合 **FEM (COMSOL)** 理论结果，揭示了微观构型差异与宏观皮腔场增强的耦合机制。

> **注意**：本仓库主要包含 DFT 计算的输入文件 (`.in`) 和后处理脚本 (`.py`)。FEM (COMSOL) 仿真由于二进制文件过大未上传，仅包含结果数据/图片。

## 📂 项目结构 (Repository Structure)

```text
.
├── scf_scan/               # [核心] 最终版平躺构型势能面扫描 (0.25Å 步长)
│   ├── scan_flat_*.in      # 不同高度下的 SCF 输入文件
│   └── 1.py                # 数据提取与绘图脚本
│
├── scan_10_v/              # [核心] 垂直构型势能面扫描 (用于对比)
│   ├── ag_h2_*_6k.relax.in # 垂直构型输入文件
│   └── ...
│
├── charge/                 # [核心] 差分电荷密度计算 (CDD)
│   ├── charge_v/           # 垂直构型电荷密度
│   ├── H2.in, Ag.in...     # 平躺构型电荷密度计算 (Total/Slab/Mol)
│   ├── pp_*.in             # 后处理输入文件 (生成 .cube)
│   └── calc_diff.py        # 差分电荷计算脚本 (Total - Slab - Mol)
│
├── ag_pes_scan/            # 高密度扫描测试数据 (20个点)
│
├── Eslab/ & EH2/ & Etotal/ # 基准能量计算 (孤立表面、孤立分子、弛豫测试)
│
├── old1114/                # 早期测试与收敛性测试数据 (Archive)
│
├── *.png                   # [结果] 势能面扫描曲线图与结果可视化
└── README.md
```

## 🚀 主要结果 (Key Results)

### 1. 物理吸附势能面 (PES)
*   **文件**: `H2_on_Ag111_PES_Final.png`
*   **结论**: H₂ 在 Ag(111) 表面的基态为**平躺 (Flat)** 构型，吸附能约为 **-57.4 meV**。相比之下，垂直构型吸附能较浅 (~ -30 meV)，证明了分子在低温下倾向于平躺吸附，且存在显著的转动势垒。

### 2. 差分电荷密度 (CDD)
*   **路径**: `charge/`
*   **结论**: 界面处观察到明显的电子耗尽（Depletion，青色）与两侧的电子积聚（Accumulation，黄色），证实了**泡利排斥 (Pauli Repulsion)** 和电子回推效应是弱物理吸附的主要成因。

### 3. 皮腔场增强 (FEM Theory)
*   虽然未上传 `.mph` 文件，但仿真结果证实原子级突起将光场压缩至 < 1 nm 区域，产生高达 **$10^{14}$** 的电场增强 ($|E|^4$)，解释了 H₂ (低密度/自由) 与 D₂ (高密度/受限) 的光谱差异。

## 🛠️ 使用说明 (Usage)

### 环境依赖
*   Quantum ESPRESSO (v7.0+)
*   Python 3.x (`numpy`, `matplotlib`, `scipy`)

### 运行势能面扫描
进入 `scf_scan` 目录，使用脚本批量提交计算：
```bash
cd scf_scan
bash 1.sh
```

### 绘制结果图
使用根目录或子目录下的 Python 脚本读取输出数据并绘图：
```bash
python scf_scan/1.py
```

### 计算差分电荷
进入 `charge` 目录，依次运行 SCF 计算、PP 后处理，最后运行 Python 脚本做减法：
```bash
cd charge
# 1. 运行 pw.x (Total, Ag, H2)
# 2. 运行 pp.x (提取电荷)
python calc_diff.py  # 生成 rho_diff.cube
```

## 📝 引用 (References)
*   [1] A. Shiotari, et al., *Phys. Rev. Lett.* **134**, 206901 (2025).

## 🤝 作者 (Author)
*   Course Project: 2025 PhD Condensed Matter Theory
*   Author: [Your Name]