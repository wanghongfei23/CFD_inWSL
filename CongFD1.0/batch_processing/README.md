# 批量处理配置文件说明

## 重要提示

配置文件必须严格按照16行格式编写，缺少任何一行都可能导致程序读取错误或运行异常。

## 配置文件命名规范

配置文件应按照以下格式命名：
```
info_"算例名称"_"插值方法".txt
```

例如：
- `info_Lax_WCNS5CONGZ.txt`
- `info_Lax_WHFTCNSA.txt`

## 配置文件内容格式

配置文件**必须**包含16行，每行一个参数，顺序如下：

1. 插值方法编号 (interMethod)
2. 时间步数 (endStep)
3. 输出时间间隔 (outputDt)
4. CFL数 (CFL)
5. 算例编号 (nCase)
6. 计算区域左边界 (calZone[0])
7. 计算区域右边界 (calZone[1])
8. 计算区域下边界 (calZone[2])
9. 计算区域上边界 (calZone[3])
10. 计算区域z方向起始 (calZone[4])
11. 计算区域z方向结束 (calZone[5])
12. X方向网格点数 (iMax[0])
13. Y方向网格点数 (iMax[1])
14. Z方向网格点数 (iMax[2])
15. 维度 (dim)
16. 线程数 (omp_threads)

## 参数说明

### 1. 插值方法编号 (interMethod)

| 编号 | 插值方法       | 说明              |
|------|----------------|-------------------|
| 0    | FIRSTORDER     | 一阶迎风格式      |
| 1    | MUSCL          | MUSCL格式         |
| 2    | WCNS5          | WENO-JS           |
| 3    | WCNSZ5         | WENO-Z            |
| 4    | WCNS5Char      | WCNS5特征重构     |
| 5    | WCNSZ5Char     | WCNSZ5特征重构    |
| 6    | WCNS5CONG      | WCNS5改进格式     |
| 7    | TCNSCongA      | TENO-A格式        |
| 8    | WCNS5CONGZ     | TENO-S          |
| 9    | WCNS5CONGZCT4  | TENO-S CT4版本  |
| 10   | WCNS5CONGZCT7  | TENO-S CT7版本  |
| 11   | TCNS5          | TENO            |
| 12   | TCNS5CT4       | TENO CT4版本    |
| 13   | TCNS5CT7       | TENO CT7版本    |
| 14   | LINEAR5        | 五阶线性格式      |
| 15   | NICEST5        | NICEST格式        |
| 16   | INTERMAX       | 最大值插值格式    |
| 17   | WHFTCNSA       | TENO-myA        |
| 18   | WHFTCNSASF002   | TENO-AS-myF002   |
| 19   | WHFTCNSAH002   | TENO-AS-myH002   |
| 20   | WHFTCNSASF102   | TENO-AS-myF102   |
| 21   | WHFTCNSASF103   | TENO-AS-myF103   |

### 2. 算例编号 (nCase)

#### 一维算例

| 编号 | 算例名称           | 说明             |
|------|--------------------|------------------|
| 0    | Sod                | Sod激波管        |
| 1    | ShuOsher           | Shu-Osher问题    |
| 2    | Lax                | Lax激波管        |
| 3    | sedov              | Sedov爆炸波      |
| 4    | Woodward_Colella   | Woodward-Colella问题 |
| 5    | Double_sparse_wave | 双稀疏波问题     |

#### 二维算例

| 编号 | 算例名称              | 说明              |
|------|-----------------------|-------------------|
| 0    | 2D_Riemann_1          | 二维黎曼问题1     |
| 1    | 2D_Riemann_2   | 二维黎曼问题2-涡 |
| 2    | implosion             | 内爆问题          |
| 3    | RTI                   | 瑞利-泰勒不稳定性 |
| 4    | Double_Mach           | 双马赫反射        |
| 5    | 2D_Riemann_3  | 二维黎曼问题3     |
| 6    | KHI                   | 开尔文-亥姆霍兹不稳定性 |

## 一维与二维问题参数设置差异

虽然一维和二维问题使用相同的16行配置文件格式，但在参数设置上有以下区别：

### 一维问题 (dim=1)
- calZone[2] 到 calZone[5] 应设置为 0
- iMax[1] 和 iMax[2] 通常设置为 2
- 示例：
  ```
  8        # interMethod
  100      # endStep
  0.01     # outputDt
  0.5      # CFL
  2        # nCase (Lax)
  -0.5     # calZone[0]
  0.5      # calZone[1]
  0        # calZone[2] (unused)
  0        # calZone[3] (unused)
  0        # calZone[4] (unused)
  0        # calZone[5] (unused)
  201      # iMax[0]
  2        # iMax[1] (unused)
  2        # iMax[2] (unused)
  1        # dim
  4        # omp_threads
  ```

### 二维问题 (dim=2)
- calZone[0] 到 calZone[3] 分别表示 x 和 y 方向的边界
- iMax[0] 和 iMax[1] 分别表示 x 和 y 方向的网格点数
- 示例：
  ```
  2        # interMethod
  100      # endStep
  0.01     # outputDt
  0.5      # CFL
  0        # nCase (2D_Riemann_1)
  -0.5     # calZone[0] (x_min)
  0.5      # calZone[1] (x_max)
  -0.5     # calZone[2] (y_min)
  0.5      # calZone[3] (y_max)
  0        # calZone[4] (unused)
  0        # calZone[5] (unused)
  201      # iMax[0] (x方向网格点数)
  201      # iMax[1] (y方向网格点数)
  2        # iMax[2] (unused)
  2        # dim
  4        # omp_threads
  ```

## 使用方法

1. 创建符合命名规范的配置文件
2. 将配置文件放入 `batch_processing` 目录
3. 运行批量处理脚本：
   ```bash
   cd batch_processing
   ./run_all_cases.sh
   ```
4. 结果文件将保存在 `batch_processing_data` 目录中