# code: Analysis of pancreatic cancer bioinformation
                                   
## file: cfg.py

code文件夹的配置文件

## file: getexpmatrix.py:

- 合并表达文件
- 将ENS编码替换成基因名

## file: GSVA.R

- GSVA 的R语言核心文件
- 接受参数：
	- args[1]：mRNA的表达文件
	- args[2]：Gmt文件
	- args[3]：gsva生成文件的存储路径

## file: gsva.py

- gsva运行文件，生成一个gsva表达文件和聚类文件
- 接受参数：-g + gmt文件路径
- 目前的默认mRNA表达文件在CFG.EXPRESSION_FILE
- gsva结果输出在：CFG.resultpath + gsva_result.csv'
- 聚类结果文件默认保存在CFG.resultpath + ClusterResult/
	- 聚类数量：CFG.MIN_CLASS -> CFG.MAX_CLASS
- 聚类热图默认保存在CFG.resultpath + gsva_heatmap.pdf

## file: SurvivalAnalysis.py

- 生存分析的数据预处理
- 不接受参数，默认使用不同结果文件夹中ClusterResult/中的cluster_n文件，作为类别标记
- 将聚类结果添加临床数据和生存数据，结果保存在ClusterResult/中，文件名With_class_cluster_n.csv
- 生成的文件使用SurvivalAnalysis.R进行后续分析

## file: SurvivalAnalysis.R

- 生存分析的R语言运行文件
- 接受1个参数：

## file cibersort.R

- cibersortX分析免疫细胞组成
- 将结果存储在CIBERSORT_RESULT.csv中，供后续应用

## file CIBERSORT.py

- 计算cibersort所得到免疫评分的差异，需要进一步修改

## file limma.R

- 计算limma结果
- 三个参数
	- 表达文件（或者GSVA文件）
	- 分类文件（两列，第1列case_id，第2列类别号）
	- 类别数量


