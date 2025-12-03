#!/usr/bin/env python3
"""
高级代码分析器，用于提取C++源代码中的详细函数调用关系和文件依赖关系，
并生成DOT格式的图形描述文件。
"""

import os
import re
import argparse
from collections import defaultdict


class AdvancedCodeAnalyzer:
    def __init__(self, src_root):
        self.src_root = src_root
        # 存储文件之间的包含关系
        self.file_dependencies = defaultdict(set)
        # 存储函数调用关系
        self.function_calls = defaultdict(lambda: defaultdict(set))
        # 存储函数定义位置
        self.function_definitions = {}
        # 存储类定义位置
        self.class_definitions = {}
        # 存储命名空间信息
        self.namespace_info = defaultdict(dict)

    def analyze(self):
        """
        分析整个源代码树
        """
        # 遍历所有头文件和源文件
        for root, dirs, files in os.walk(self.src_root):
            for file in files:
                if file.endswith(('.hpp', '.h', '.cpp', '.c')):
                    file_path = os.path.join(root, file)
                    self._analyze_file(file_path)

    def _analyze_file(self, file_path):
        """
        分析单个文件
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
        except UnicodeDecodeError:
            # 如果UTF-8解码失败，尝试其他编码
            try:
                with open(file_path, 'r', encoding='latin-1') as f:
                    content = f.read()
            except Exception:
                print(f"警告：无法读取文件 {file_path}")
                return

        # 提取#include语句（文件依赖）
        self._extract_includes(file_path, content)
        
        # 提取命名空间信息
        self._extract_namespace_info(content)
        
        # 提取函数定义
        self._extract_function_definitions(file_path, content)
        
        # 提取函数调用
        self._extract_function_calls(file_path, content)
        
        # 提取类定义
        self._extract_class_definitions(file_path, content)

    def _extract_includes(self, file_path, content):
        """
        提取#include语句
        """
        include_pattern = r'#include\s*[<"]([^>"]+)[>"]'
        matches = re.findall(include_pattern, content)
        for include in matches:
            self.file_dependencies[os.path.basename(file_path)].add(include)

    def _extract_namespace_info(self, content):
        """
        提取命名空间信息
        """
        namespace_pattern = r'namespace\s+(\w+)\s*\{'
        matches = re.findall(namespace_pattern, content)
        # 这里简化处理，实际应该跟踪命名空间范围

    def _extract_function_definitions(self, file_path, content):
        """
        提取函数定义
        """
        # 移除注释和字符串字面量以避免误匹配
        cleaned_content = self._clean_content(content)
        
        # 匹配函数定义的模式（包括在类内的函数）
        patterns = [
            # 普通函数定义
            r'(\w+(?:\s*(?:\*|&|\bconst\b))*)\s+(\w+(?:::\w+)?)\s*\([^)]*\)\s*(?:const)?\s*\{',
            # 构造函数/析构函数定义
            r'\b(\w+)::(\w+)\s*\([^)]*\)\s*(?:const)?\s*\{',
        ]
        
        for pattern in patterns:
            matches = re.finditer(pattern, cleaned_content)
            for match in matches:
                if pattern == patterns[0]:
                    return_type = match.group(1).strip()
                    full_func_name = match.group(2).strip()
                else:
                    class_name = match.group(1).strip()
                    func_name = match.group(2).strip()
                    full_func_name = f"{class_name}::{func_name}"
                    # 对于构造函数/析构函数，我们假设返回类型是void
                    return_type = "void"
                
                # 忽略控制结构关键字
                base_func_name = full_func_name.split('::')[-1]
                if base_func_name in ['if', 'for', 'while', 'switch', 'catch', 'return', 'delete', 'new']:
                    continue
                    
                self.function_definitions[full_func_name] = {
                    'file': os.path.basename(file_path),
                    'return_type': return_type
                }

    def _extract_function_calls(self, file_path, content):
        """
        提取函数调用
        """
        # 移除注释和字符串字面量以避免误匹配
        cleaned_content = self._clean_content(content)
        
        # 匹配函数调用的模式
        # 匹配一般的函数调用 func(...)
        call_pattern = r'(?<!\w)(\w+(?:::\w+)*)\s*\('
        matches = re.finditer(call_pattern, cleaned_content)
        
        source_file = os.path.basename(file_path)
        
        for match in matches:
            func_call = match.group(1)
            
            # 忽略一些常见的语言关键字
            base_name = func_call.split('::')[-1]
            keywords = {'if', 'for', 'while', 'switch', 'sizeof', 'typeof', 'return', 
                       'delete', 'new', 'catch', 'throw', 'explicit', 'inline', 'virtual',
                       'const', 'static', 'extern', 'friend', 'operator', 'template', 'typename'}
            
            if base_name in keywords:
                continue
                
            # 记录函数调用
            if func_call in self.function_definitions:
                target_file = self.function_definitions[func_call]['file']
                self.function_calls[source_file][target_file].add(func_call)
            elif '::' in func_call:
                # 处理命名空间或类成员函数调用
                base_name = func_call.split('::')[0]
                if base_name in self.function_definitions:
                    target_file = self.function_definitions[base_name]['file']
                    self.function_calls[source_file][target_file].add(func_call)
            else:
                # 查找是否有匹配的函数
                matched = False
                for defined_func in self.function_definitions:
                    if defined_func.endswith('::' + func_call):
                        target_file = self.function_definitions[defined_func]['file']
                        self.function_calls[source_file][target_file].add(func_call)
                        matched = True
                        break
                if not matched:
                    # 未能匹配到任何定义的函数，暂时跳过
                    pass

    def _extract_class_definitions(self, file_path, content):
        """
        提取类定义
        """
        class_pattern = r'(?:class|struct)\s+(\w+)'
        matches = re.findall(class_pattern, content)
        for class_name in matches:
            self.class_definitions[class_name] = os.path.basename(file_path)

    def _clean_content(self, content):
        """
        清理内容，移除注释和字符串以避免误匹配
        """
        # 移除单行注释
        content = re.sub(r'//.*?$', '', content, flags=re.MULTILINE)
        # 移除多行注释
        content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)
        # 简化处理：移除字符串字面量（但保留函数调用结构）
        content = re.sub(r'"(?:\\.|[^"\\])*"', '""', content)
        content = re.sub(r"'(?:\\.|[^'\\])*'", "''", content)
        return content

    def generate_dot_file(self, output_file):
        """
        生成DOT格式的图形描述文件
        """
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("digraph G {\n")
            f.write("  rankdir=TB;\n")  # 从上到下布局
            f.write("  node [shape=box, style=filled, fillcolor=lightgray];\n")
            f.write("  compound=true;\n\n")
            
            # 收集所有涉及的文件
            all_files = set()
            all_files.update(self.file_dependencies.keys())
            for source, targets in self.function_calls.items():
                all_files.add(source)
                all_files.update(targets.keys())
            
            # 为每个文件创建节点
            file_nodes = {}
            for i, file in enumerate(all_files):
                node_name = f"file_{i}"
                file_nodes[file] = node_name
                f.write(f'  {node_name} [label="{file}", fillcolor=white];\n')
            
            f.write("\n")
            
            # 写入文件依赖关系（边）
            f.write("  // 文件依赖关系\n")
            for file, includes in self.file_dependencies.items():
                if file in file_nodes:
                    source_node = file_nodes[file]
                    for include in includes:
                        # 只处理看起来像本地文件的依赖项
                        if not include.startswith('<') and include in file_nodes:
                            target_node = file_nodes[include]
                            f.write(f'  {source_node} -> {target_node} [color=blue, style=dashed, arrowhead=empty];\n')
            
            f.write("\n")
            
            # 写入函数调用关系
            f.write("  // 函数调用关系\n")
            for source_file, targets in self.function_calls.items():
                if source_file in file_nodes:
                    source_node = file_nodes[source_file]
                    for target_file, calls in targets.items():
                        if target_file in file_nodes:
                            target_node = file_nodes[target_file]
                            # 合并标签以减少混乱
                            call_list = ', '.join(list(calls)[:5])  # 最多显示5个调用
                            if len(calls) > 5:
                                call_list += f" ...({len(calls)} total)"
                            
                            f.write(f'  {source_node} -> {target_node} [label="{call_list}", color=red, fontsize=10];\n')
            
            f.write("\n")
            
            # 添加图例
            f.write("  subgraph cluster_legend {\n")
            f.write("    label=Legend;\n")
            f.write("    style=dotted;\n")
            f.write('    key1 [shape=box, label="文件", fillcolor=white];\n')
            f.write('    key2 [shape=box, style=dashed, label="文件包含关系", color=blue];\n')
            f.write('    key3 [shape=plaintext, label="函数调用", color=red];\n')
            f.write("  }\n")
            
            f.write("}\n")

    def print_summary(self):
        """
        打印分析摘要
        """
        total_deps = sum(len(deps) for deps in self.file_dependencies.values())
        total_calls = sum(len(calls) for targets in self.function_calls.values() for calls in targets.values())
        
        print(f"分析完成:")
        print(f"- 发现 {len(self.file_dependencies)} 个文件具有包含依赖")
        print(f"- 总共 {total_deps} 个文件包含关系")
        print(f"- 发现 {len(self.function_definitions)} 个函数定义")
        print(f"- 检测到 {total_calls} 次函数调用")
        print(f"- 发现 {len(self.class_definitions)} 个类定义")


def main():
    parser = argparse.ArgumentParser(description='高级C++代码分析器，生成函数调用和文件依赖关系图')
    parser.add_argument('--src-root', '-s', default='..', help='源代码根目录')
    parser.add_argument('--output', '-o', default='advanced_code_structure.dot', help='输出DOT文件名')
    
    args = parser.parse_args()
    
    # 检查源代码目录是否存在
    if not os.path.exists(args.src_root):
        print(f"错误：源代码目录 {args.src_root} 不存在")
        return
    
    print("开始分析代码...")
    analyzer = AdvancedCodeAnalyzer(args.src_root)
    analyzer.analyze()
    analyzer.generate_dot_file(args.output)
    analyzer.print_summary()
    
    print(f"\nDOT文件已生成: {args.output}")
    print("可以使用Graphviz工具将其转换为图像:")
    print(f"  dot -Tpng {args.output} -o advanced_code_structure.png")
    print(f"  dot -Tsvg {args.output} -o advanced_code_structure.svg")


if __name__ == '__main__':
    main()