# Python 代码分析工具

这个目录包含了用于分析 CFD 代码库的 Python 工具。

## 工具介绍

### 1. 代码分析器.py - 基础代码分析器

基础版本的代码分析工具，可以扫描 C++ 源代码并提取：
- 文件依赖关系 (#include 关系)
- 函数定义位置
- 函数调用关系
- 类定义位置

使用方法：
```bash
cd py
python3 代码分析器.py
```

默认会分析 `../src` 目录下的代码，并生成 `code_structure.dot` 文件。

### 2. 高级代码分析器.py - 高级代码分析器

增强版的代码分析工具，提供了更加准确的分析能力：
- 更准确地解析函数定义和调用
- 更好地处理命名空间和类成员函数
- 生成更清晰的关系图

使用方法：
```bash
cd py
python3 高级代码分析器.py
```

默认会分析 `../src` 目录下的代码，并生成 `advanced_code_structure.dot` 文件。

## 生成图形文件

两个工具都会生成 DOT 格式的文件，可以使用 Graphviz 工具转换为图像：

```bash
# 生成 PNG 图像
dot -Tpng code_structure.dot -o code_structure.png

# 生成 SVG 图像
dot -Tsvg code_structure.dot -o code_structure.svg

# 生成 PDF 文档
dot -Tpdf code_structure.dot -o code_structure.pdf
```

如果系统没有安装 Graphviz，可以通过以下命令安装：

Ubuntu/Debian:
```bash
sudo apt-get install graphviz
```

CentOS/RHEL:
```bash
sudo yum install graphviz
```

macOS:
```bash
brew install graphviz
```

Windows:
可以从 Graphviz 官网下载安装包: https://graphviz.org/download/

## 自定义参数

可以通过命令行参数自定义分析行为：

```bash
# 指定源代码目录和输出文件
python3 高级代码分析器.py --src-root ../src --output my_analysis.dot

# 简短形式
python3 代码分析器.py -s ../other_src -o analysis.dot
```