#!/bin/bash

# 脚本说明：此脚本用于自动运行batch_processing目录中所有配置文件的算例
# 作者：Assistant
# 日期：$(date +%Y-%m-%d)

echo "==============================================="
echo "CFD全自动多算例连续运行脚本"
echo "==============================================="

# 获取脚本所在目录作为工作目录
WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "工作目录: $WORK_DIR"

# 切换到工作目录
cd "$WORK_DIR" || exit
echo "当前工作目录: $(pwd)"

# 检查是否存在配置文件
if ! ls info_*.txt 1> /dev/null 2>&1; then
    echo "错误: 当前目录中未找到任何配置文件 (info_*.txt)"
    exit 1
fi

# 检查可执行文件是否存在（相对于项目根目录）
PROJECT_ROOT="$(dirname "$WORK_DIR")"
if [ ! -f "$PROJECT_ROOT/build/congFD" ]; then
    echo "错误: 可执行文件build/congfd不存在，请先编译项目"
    echo "编译方法:"
    echo "  cd $PROJECT_ROOT"
    echo "  mkdir -p build && cd build"
    echo "  cmake .."
    echo "  make -j$(nproc)"
    exit 1
fi

# 定义结果目录
DATA_DIR="$PROJECT_ROOT/batch_processing_data"

# 创建结果目录（如果不存在）
mkdir -p "$DATA_DIR"

# 运行前清理旧的结果文件（可选）
echo "清理旧结果文件..."
rm -f result_*.dat timeInfo_*.txt
rm -f "$DATA_DIR"/result_*.dat "$DATA_DIR"/timeInfo_*.txt
# 清理可能存在的CGNS文件（只清理结果目录）
rm -f "$DATA_DIR"/*.cgns

# 获取所有配置文件列表
config_files=(info_*.txt)

# 遍历所有配置文件并运行算例
for config_file in "${config_files[@]}"; do
    # 从配置文件名提取算例名称（去除前缀"info_"和后缀".txt"）
    case_name="${config_file#info_}"
    case_name="${case_name%.txt}"
    
    echo "==============================================="
    echo "正在运行 $case_name 算例..."
    echo "==============================================="

    # 复制配置文件
    cp "$config_file" info.txt
    echo "已加载 $case_name 配置文件"

    # 记录开始时间
    START_TIME=$(date +%s)

    # 运行程序
    echo "开始执行程序..."
    "$PROJECT_ROOT/build/congFD"

    # 记录结束时间
    END_TIME=$(date +%s)
    EXEC_TIME=$((END_TIME - START_TIME))

    # 检查并移动timeInfo.txt文件
    if [ -f "timeInfo.txt" ]; then
        mv timeInfo.txt "$DATA_DIR/timeInfo_${case_name}.txt"
        echo "时间信息已保存为 $DATA_DIR/timeInfo_${case_name}.txt"
    else
        echo "警告: 未找到 timeInfo.txt 文件"
    fi

    # 查找并移动生成的CGNS文件（使用更安全的方式处理文件名）
    for cgns_file in *.cgns; do
        if [ -f "$cgns_file" ]; then
            mv "$cgns_file" "$DATA_DIR/"
            echo "CGNS文件已保存为 $DATA_DIR/$cgns_file"
        fi
    done

    echo "$case_name 算例完成，耗时: $EXEC_TIME 秒"
    echo ""
done

# 完成提示
echo "==============================================="
echo "所有算例运行完成！"
echo "时间信息文件:"
for config_file in "${config_files[@]}"; do
    case_name="${config_file#info_}"
    case_name="${case_name%.txt}"
    if [ -f "$DATA_DIR/timeInfo_${case_name}.txt" ]; then
        echo "  - $DATA_DIR/timeInfo_${case_name}.txt"
    fi
done

echo "CGNS文件:"
for cgns_file in "$DATA_DIR"/*.cgns; do
    if [ -f "$cgns_file" ]; then
        echo "  - $cgns_file"
    fi
done
echo "==============================================="