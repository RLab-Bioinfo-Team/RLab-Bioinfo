# scRNA-seq_Pipeline
RLab 单核 RNA-seq (snRNA-seq) 标准分析流程。包含：分析脚本、测试数据与结果展示。

## 目录结构
- `01_pipeline/`: 核心分析脚本 (Jupyter Notebook 格式)
- `02_test_data/`: 抽样测试矩阵与同源基因映射参考
- `03_results/`: 预期输出的表格文件与 Marker 基因列表

## 环境依赖
* **Python**: Scanpy (用于 Part 1 Scrublet 双细胞核过滤)
* **R**: Seurat v5, Harmony, decontX (用于 Part 2-4 降维、去环境 RNA 与注释)

## 快速开始 (Quick Start)
1. **解压数据**：进入 `02_test_data`，将 `T1_test_feature_bc_matrix.zip` 和 `T2_test_feature_bc_matrix.zip` 解压至工作目录。
2. **运行流程**：按顺序执行 `01_pipeline` 中的 Part 1 至 Part 4 脚本。
   - **Part 1**: Scrublet 质控与过滤
   - **Part 2**: Seurat 对象创建与线粒体过滤
   - **Part 3**: Harmony 批次效应校正
   - **Part 4**: decontX 环境 RNA 去除与细胞核 (nuclei) 亚群注释
