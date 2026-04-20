# Troubleshooting & FAQ

本文件记录 RLab snRNA-seq 流程运行过程中的常见报错及解决方案。

## [ERROR 01] Seurat v5 Harmony 整合报错
**现象**：提示无法在对象中找到对应的 Layers。

**解决**：在 v5 版本中，执行 `IntegrateLayers` 前必须先确保数据已通过 `v5_obj[["RNA"]] <- split(v5_obj[["RNA"]], f = v5_obj$sample)` 进行拆分。

## [ERROR 02] 基因映射字典匹配为空
**现象**：`oni2zebra.txt` 运行后没有匹配到基因名。

**解决**：检查输入矩阵的基因 ID 是否带有版本号后缀（如 .1, .2），需先用 gsub 剔除后再进行字典匹配。
