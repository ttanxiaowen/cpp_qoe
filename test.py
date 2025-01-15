import os
import re


import os

def combine_txt_files(source_directory, output_file):
    # 检查源目录是否存在
    if not os.path.exists(source_directory):
        print("Source directory does not exist.")
        return

    # 确保输出文件的目录存在
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 打开输出文件
    with open(output_file, 'w', encoding='utf-8') as outfile:
        # 遍历源目录中的所有文件
        for filename in os.listdir(source_directory):
            # 检查文件是否为.txt文件
            if filename.endswith('.txt'):
                # 获取文件的完整路径
                file_path = os.path.join(source_directory, filename)
                # 打开并读取.txt文件内容
                with open(file_path, 'r', encoding='utf-8') as infile:
                    contents = infile.read()
                    # 将内容写入输出文件
                    outfile.write(contents + '\n')  # 添加换行符以分隔不同文件的内容

    print(f"All .txt files have been combined into '{output_file}'.")

# 使用示例
source_directory = 'build/link'  # 替换为你的目录路径
output_file = 'link.txt'  # 替换为你的输出文件路径
combine_txt_files(source_directory, output_file)
# 使用你的目录路径替换这里的 'path/to/your/directory'