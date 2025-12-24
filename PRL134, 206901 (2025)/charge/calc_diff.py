import numpy as np
import sys
import os

def read_cube(filename):
    """
    读取 Gaussian Cube 文件
    返回: 
        header: 头信息列表 (包含原子坐标等)
        data: 网格数据 (numpy array)
    """
    print(f"正在读取 {filename} ...")
    if not os.path.exists(filename):
        print(f"错误: 找不到文件 {filename}")
        sys.exit(1)

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Cube 文件的前两行是注释
    # 第三行包含原子数量 (第一列)
    try:
        natoms_line = lines[2].split()
        natoms = int(natoms_line[0])
    except:
        print(f"错误: {filename} 格式看似不正确。")
        sys.exit(1)

    # Cube 文件的头部行数 = 6 (基础信息) + 原子数
    # 注意：有时候 natoms 是负数（表示单位标志），取绝对值
    header_lines_count = 6 + abs(natoms)
    
    header = lines[:header_lines_count]
    
    # 剩下的内容是网格数据
    # 将剩余行合并成一个长字符串，然后用 numpy 快速解析
    data_str = " ".join(lines[header_lines_count:])
    data = np.fromstring(data_str, sep=' ')
    
    return header, data

def write_cube(filename, header, data):
    """
    写入 Gaussian Cube 文件
    """
    print(f"正在写入 {filename} ...")
    with open(filename, 'w') as f:
        # 1. 写入头信息
        f.writelines(header)
        
        # 2. 写入数据
        # Cube 标准格式：每行 6 个数据，科学计数法宽格式
        # 这种写法比逐个循环快得多
        for i in range(0, len(data), 6):
            chunk = data[i:i+6]
            line = "".join([f"{val:13.5E}" for val in chunk])
            f.write(line + "\n")

def main():
    # --- 配置输入文件名 ---
    file_total = "rho_total.cube"
    file_Ag    = "rho_Ag.cube"
    file_H2    = "rho_H2.cube"
    file_out   = "rho_diff.cube"

    # --- 1. 读取数据 ---
    # 我们使用 total 的 header 作为输出文件的 header
    # 这样在 VESTA 里显示时，会包含所有的原子 (Ag + H)
    header, rho_total = read_cube(file_total)
    _, rho_Ag    = read_cube(file_Ag)
    _, rho_H2    = read_cube(file_H2)

    # --- 2. 安全检查 ---
    # 确保三个文件的网格点数是一样的，否则不能相减
    if rho_total.shape != rho_Ag.shape or rho_total.shape != rho_H2.shape:
        print("错误: 输入文件的网格大小不一致！")
        print(f"Total grid: {rho_total.shape}")
        print(f"Ag grid:    {rho_Ag.shape}")
        print(f"H2 grid:    {rho_H2.shape}")
        print("请检查 scf 计算时是否使用了相同的晶胞参数和截断能。")
        sys.exit(1)

    # --- 3. 计算差分电荷密度 ---
    # 公式: Delta_rho = rho_total - rho_Ag - rho_H2
    print("正在计算差分: Total - Ag - H2 ...")
    rho_diff = rho_total - rho_Ag - rho_H2

    # --- 4. 保存结果 ---
    write_cube(file_out, header, rho_diff)
    print(f"成功! 差分文件已保存为: {file_out}")
    print("现在可以将其拖入 VESTA 进行可视化分析了。")

if __name__ == "__main__":
    main()