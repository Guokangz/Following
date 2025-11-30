from ase.io import read, write
import sys

# 检查命令行是否提供了输入文件名
if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    # 如果没有提供，就默认使用 ag.relax.out
    input_file = 'ag.relax.out'
    print(f"Usage: python convert.py <your_output_file.out>")
    print(f"No input file provided. Using default: {input_file}")

# 定义输出文件名
output_file = input_file.replace('.out', '.cif')

try:
    # 从 QE 的 .out 文件中读取最后一个（能量最低的）结构
    atoms = read(input_file, format='espresso-out')

    # 将结构写入 .cif 文件
    write(output_file, atoms)

    print(f"SUCCESS: Successfully converted '{input_file}' to '{output_file}'")

except Exception as e:
    print(f"ERROR: An error occurred during conversion.")
    print(e)
