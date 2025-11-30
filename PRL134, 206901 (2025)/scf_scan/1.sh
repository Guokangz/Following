#!/bin/bash

# ==============================================================================
# 自动化势能面扫描脚本 (Potential Energy Surface Scan)
# 针对体系：H2 (平躺) on Ag(111)
# ==============================================================================

# --- 1. 并行与环境设置 ---
PW_COMMAND="pw.x"
# 根据你的机器实际情况调整核心数，这里设为 16
MPI_COMMAND="mpirun -np 16 --bind-to core --map-by core"

# 赝势文件名 (请确保这些文件在当前目录下)
AG_PSEUDO="Ag.pbe-n-rrkjus_psl.1.0.0.UPF"
H_PSEUDO="H.pbe-kjpaw.UPF"

# --- 2. 扫描参数设置 ---
# 扫描点数
N_POINTS=15
# 起始距离 (埃) - 离表面比较近，斥力区
START_DIST=2.0
# 结束距离 (埃) - 离表面比较远，接近真空
END_DIST=5.5

# --- 3. 物理常数 (基于你的晶胞设置) ---
# 晶胞 Z 轴总长度 (Angstrom)
CELL_Z=25.0
# 最上层 Ag 原子的分数坐标 Z (来自你的结构文件)
# 0.344028 * 25.0 = 8.6007 A
Z_AG_TOP=0.344028

# --- 4. 结果汇总文件 ---
SUMMARY_FILE="scan_results.dat"
echo "# Distance(A)   Total_Energy(Ry)   E_Adsorption(meV_approx)" > $SUMMARY_FILE

# --- 5. 开始循环扫描 ---
echo "======================================================"
echo "开始 H2 平躺构型 Z 轴扫描"
echo "扫描范围: $START_DIST A -> $END_DIST A"
echo "总计算点数: $N_POINTS"
echo "======================================================"

for i in $(seq 0 $((N_POINTS-1)))
do
    # 1. 计算当前步的物理距离 (使用 bc 计算浮点数)
    # dist = start + (end - start) * i / (n - 1)
    CURRENT_DIST=$(echo "scale=6; $START_DIST + ($END_DIST - $START_DIST) * $i / ($N_POINTS - 1)" | bc -l)
    
    # 2. 计算 H 分子在晶胞中的分数坐标 Z
    # Z_frac = Z_Ag_top + (Distance / Cell_Z)
    Z_H2=$(echo "scale=8; $Z_AG_TOP + $CURRENT_DIST / $CELL_Z" | bc -l)

    # 3. 定义文件名 (例如 dist_2.50.in)
    LABEL=$(printf "%.2f" $CURRENT_DIST)
    PREFIX="scan_flat_${LABEL}"
    INPUT_FILE="${PREFIX}.in"
    OUTPUT_FILE="${PREFIX}.out"

    echo "正在计算点 $((i+1))/$N_POINTS: 距离 = $CURRENT_DIST A (Z_frac = $Z_H2)..."

    # 4. 生成输入文件
    # 注意 H 原子的坐标设置：
    # - 两个 H 的 Z 坐标相同 (平躺)
    # - 两个 H 的 X 坐标拉开约 0.08 (对应键长约 0.74A)
    # - 标志位 '1 1 0' 表示: X弛豫, Y弛豫, Z固定 (关键!)
    cat > $INPUT_FILE << EOF
&CONTROL
   calculation  = 'relax',
   prefix       = '$PREFIX',
   pseudo_dir   = './',
   outdir       = './temp_dir/',
   tstress      = .true.,
   tprnfor      = .true.,
   disk_io      = 'low'
/
&SYSTEM
   ibrav=0, nat=29, ntyp=2, 
   ecutwfc=40.0, ecutrho=200.0,
   occupations='smearing', smearing='methfessel-paxton', degauss=0.02,
   vdw_corr = 'dft-d3'
/
&ELECTRONS
   conv_thr    = 1.0d-8,
   mixing_beta = 0.3
/
&IONS
   ion_dynamics = 'bfgs'
/
CELL_PARAMETERS {angstrom}
   8.850   0.000   0.000
  -4.425   7.665   0.000
   0.000   0.000  25.000
ATOMIC_SPECIES
   Ag  107.868  $AG_PSEUDO
   H    1.008   $H_PSEUDO
ATOMIC_POSITIONS {crystal}
Ag   0.2222220000   0.1111110000   0.2000000000   0 0 0
Ag   0.5555560000   0.1111110000   0.2000000000   0 0 0
Ag   0.8888890000   0.1111110000   0.2000000000   0 0 0
Ag   0.2222220000   0.4444440000   0.2000000000   0 0 0
Ag   0.5555560000   0.4444440000   0.2000000000   0 0 0
Ag   0.8888890000   0.4444440000   0.2000000000   0 0 0
Ag   0.2222220000   0.7777780000   0.2000000000   0 0 0
Ag   0.5555560000   0.7777780000   0.2000000000   0 0 0
Ag   0.8888890000   0.7777780000   0.2000000000   0 0 0
Ag   0.1111110000   0.2222220000   0.2500000000   0 0 0
Ag   0.4444440000   0.2222220000   0.2500000000   0 0 0
Ag   0.7777780000   0.2222220000   0.2500000000   0 0 0
Ag   0.1111110000   0.5555560000   0.2500000000   0 0 0
Ag   0.4444440000   0.5555560000   0.2500000000   0 0 0
Ag   0.7777780000   0.5555560000   0.2500000000   0 0 0
Ag   0.1111110000   0.8888890000   0.2500000000   0 0 0
Ag   0.4444440000   0.8888890000   0.2500000000   0 0 0
Ag   0.7777780000   0.8888890000   0.2500000000   0 0 0
Ag   0.0000029568   0.0000059136   0.3440280909   0 0 0
Ag   0.3333370134   0.0000060547   0.3440276273   0 0 0
Ag   0.6666690413   0.0000060547   0.3440276273   0 0 0
Ag   0.0000030712   0.3333395735   0.3440283118   0 0 0
Ag   0.3333365023   0.3333395735   0.3440283118   0 0 0
Ag   0.6666702255   0.3333394510   0.3440279663   0 0 0
Ag   0.0000043613   0.6666747760   0.3440281131   0 0 0
Ag   0.3333371438   0.6666752877   0.3440280446   0 0 0
Ag   0.6666704147   0.6666747760   0.3440281131   0 0 0
H      0.460000   0.500000   $Z_H2   1 1 0
H      0.540000   0.500000   $Z_H2   1 1 0
K_POINTS {automatic}
   2 2 1 0 0 0
EOF

    # 5. 执行计算
    $MPI_COMMAND $PW_COMMAND -i $INPUT_FILE > $OUTPUT_FILE

    # 6. 检查是否正常结束并提取能量
    if grep -q "JOB DONE" $OUTPUT_FILE; then
        # 提取最终能量 (Ry)
        ENERGY=$(grep "!    total energy" $OUTPUT_FILE | tail -1 | awk '{print $5}')
        echo "  -> 计算完成. Energy = $ENERGY Ry"
        # 写入汇总文件 (仅供参考，Adsorption需后续减去E_slab和E_H2)
        echo "$CURRENT_DIST   $ENERGY" >> $SUMMARY_FILE
    else
        echo "  -> 错误: 计算未正常结束! 请检查 $OUTPUT_FILE"
        exit 1
    fi

done

echo ""
echo "全部计算完成！数据已汇总至 $SUMMARY_FILE"
echo "请用 Python 读取该文件，减去 E_slab 和 E_mol 后作图。"
