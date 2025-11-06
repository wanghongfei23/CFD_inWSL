#!/bin/bash

# 获取脚本所在目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 切换到脚本所在目录
cd "$SCRIPT_DIR" || exit 1

# 查找所有以info_开头，以.txt结尾的配置文件
for config_file in info_*.txt; do
    # 检查是否存在匹配的文件
    if [[ ! -f "$config_file" ]]; then
        echo "未找到配置文件"
        continue
    fi
    
    # 获取文件名（不含路径）
    filename=$(basename "$config_file")
    
    # 提取算例名和插值方法（从文件名中）
    # 文件名格式: info_"算例名称"_"插值方法".txt
    casename=$(echo "$filename" | cut -d'_' -f2)
    intermethod=$(echo "$filename" | cut -d'_' -f3 | cut -d'.' -f1)
    
    echo "运行算例: $casename 使用插值方法: $intermethod"
    
    # 将配置文件复制为info.txt供程序使用
    cp "$filename" info.txt
    
    # 运行主程序并直接检查执行是否成功
    if ../build/congFD; then
        echo "算例 $casename ($intermethod) 执行成功"
        # 数据保留在当前目录，不需要移动
    else
        echo "算例 $casename ($intermethod) 执行失败"
    fi
    
    echo "----------------------------------------"
done

echo "所有算例执行完毕"