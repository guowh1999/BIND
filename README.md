# BIND

BIND provides a useful tool for analyzing missing values(NAs) in a biological perspective.

It consists of 3 part, BIND_stat, BIND_classification, and BIND_ppi. The first one is BIND_stat, which works as a general statistical analysis of missing values in a proteomic dataset.

## BIND_stat

The input data are 2: expression matrix and group information.
There are 4 parameters,

---
**data_mtx_dir**: the expression data matrix file;

**grp_info_dir**: the group information file;

**grp_name**: the group type you tend to analysis, it must be one of the column names of grp_info file;

**output_dir**: analysis result output position;

---

A running example:

```
Rscript BIND_stat.R --data_mtx_dir "E:/webserver/data.csv" --grp_info_dir "E:/webserver/group_info.csv" --grp_name "Tissue" --output_dir "E:/webserver/stat"
```

Output file:

**binary_mtx.csv**: A matrix after transfering the data marix with 1/0;

**overall_stat.csv**: overall NA ratio (missing value number / total sample number);

**group_stat.csv**: NA ratio by groups.
