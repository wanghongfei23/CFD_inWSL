# 批量处理配置文件说明

## 配置文件命名规范

配置文件应按照以下格式命名：
```
info_"算例名称"_"插值方法".txt
```

例如：
- `info_Lax_WCNS5CONGZ.txt`
- `info_Lax_WHFTCNSA.txt`

## 配置文件内容格式

配置文件共有9行，每行一个参数，顺序如下：

1. 插值方法编号 (interMethod)
2. 时间步数 (endStep)
3. 输出时间间隔 (outputDt)
4. CFL数 (CFL)
5. 算例编号 (nCase)
6. 计算区域左边界 (calZone[0])
7. 计算区域右边界 (calZone[1])
8. 网格点数 (iMax[0])
9. 维度 (dim)

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
| 6    | WCNS5CONG      |                   |
| 7    | TCNSCongA      |                   |
| 8    | WCNS5CONGZ     | TENO-Z-S          |
| 9    | WCNS5CONGZCT4  |                   |
| 10   | WCNS5CONGZCT7  |                   |
| 11   | TCNS5          | TENO-Z            |
| 12   | TCNS5CT4       |                   |
| 13   | TCNS5CT7       |                   |
| 14   | LINEAR5        | 五阶线性格式      |
| 15   | NICEST5        |                   |
| 16   | WHFTCNSA       | TENO-Z-myA        |
| 17   | WHFTCNSAF002   | TENO-Z-myASF002   |

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

| 编号 | 算例名称         | 说明              |
|------|------------------|-------------------|
| 0    | 2D_Riemann_1     | 二维黎曼问题1     |
| 1    | 2D_Riemann_2_vortex | 二维黎曼问题2-涡 |
| 2    | implosion        | 内爆问题          |
| 3    | RTI              | 瑞利-泰勒不稳定性 |
| 4    | Double_Mach      | 双马赫反射        |
| 5    | 2D_Riemann_3_another | 二维黎曼问题3   |
| 6    | KHI              | 开尔文-亥姆霍兹不稳定性 |

## 使用方法

1. 创建符合命名规范的配置文件
2. 将配置文件放入 `batch_processing` 目录
3. 运行批量处理脚本：
   ```bash
   cd batch_processing
   ./run_all_cases.sh
   ```
4. 结果文件将保存在 `batch_processing_data` 目录中