
## 5. Module3: NA-weighted protein co-expression analysis.

### Input

Files:
***Expression Matrix***, ***Group Information***, and ***Query Protein List***.

***Query Protein List*** is a csv file that only contains 1 column of the interested proteins for PPI analysis. The proteins should all be in the row names of ***Expression Matrix***. No column name is needed.

User-specified parameters:

***Group name***: the column name which you interested in ***Group Information*** file. Eg. *Grouptype1*

***Group***: the specific group you are interested in. Eg. *cancer*

***reward***: the reward parameter for simultaneously  appearing NA. Eg. 0.01

***penalty***: the penalty parameter for not simultaneously appearing NA. Eg. 0.01

### What does BIND do:

---
The new correlation coefficient ρBIND was calculated based on Spearman correlation by weighting the total and different types of missing value patterns, with reward and penalty. ρBIND is more sensitive than Spearman's ρ in detecting protein-protein interactions.

---

### Output

Files:

If **Group** has been defined, the output is:

**BINDppi_groupdiff.csv**: Spearman ρ(raw_raw in the file), ρBIND(weighted_rho in file) of both target group and other samples.

**BINDppi_groupdiff_filtered.csv**: After filtering by restricted conditions, it shows target group-specific ppi found by BIND(weighted_rho). Max 200 rows. This file is for constructing a PPI network.

If **Group** has not been defined, the output are:

**BINDppi_withoutgroup.csv**: Spearman ρ(raw_raw in the file), ρBIND(weighted_rho in the file) by all samples.

**BINDppi_withoutgroup_filtered.csv**: After filtering by restricted conditions, it shows specific PPI found by BIND compared with Spearman ρ. 

Figures:

* **BIND ppi network**: A protein-protein interaction network. It shows the group-specific PPI found by ρBIND(if **Group** has been defined), or the specific ppi found by BIND compared with Spearman ρ (if **Group** has not been defined).

### Output file preview:

**BINDppi_groupdiff.csv** and **BINDppi_groupdiff_filtered.csv**

<img src="./imgs/ppi_group.png">

or

**BINDppi_withoutgroup.csv** and **BINDppi_withoutgroup_filtered.csv**

<img src="./imgs/ppi_nogroup.png">