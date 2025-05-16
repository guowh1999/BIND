## ⬇️Output explanation

### <font color= #871F78>Figures explanation</font>

#### NA statistics
**NA heatmap**
The missing value heatmap visualizes the overall pattern of missing values. It came from the *binary_matrix.csv*. Red: values; Blue: NAs.

<img src="./imgs/exa1_7.png">

<br/><br/>

**NA proportion histogram**
NA proportion histogram represents the distribution of NA proportions of different proteins in the dataset. 
X axis: NA proportion
Y axis: Protein counts
<img src="./imgs/histogram.png">

**NA proportion pie chart**
Pie chart of overall NA proportions (0 / 0-0.2 / 0.2-0.8 / 0.8-1). Hover over to see the counts.
<img src="./imgs/pie.png">

**NA proportion scatter plot**
x-y scatter plot of NA proportions vs mean expression values of a protein. Linear regression results are provided.
X axis: NA proportions in the dataset.
Y axis: Mean expression value of the proteins.

<img src="./imgs/exa1_10.png">

**Expression value with NA proportion label point plot**
Given a specific sample, ranking the protein expression values, label the overall NA proportion > 0.8/ <0.2/ 0.2-0.8 proteins by color.
X axis: proteins in the sample
Y axis: protein expression value
Red point: NA proportion > 0.8
Blue point: NA proportion < 0.2

<img src="./imgs/exa1_14.png">


**Group-specific NA proportion scatter plot**
x-y scatter plot of NA proportions in a specific group vs NA proportions in other groups.
X axis: NA proportions in the query group
Y axis: NA proportions in other samples

<img src="./imgs/exa1_15.png">

**Group level NA proportion heatmap**
Heatmap of NA proportions in different groups, clustering by samples. 
<img src="./imgs/exa1_16.png">

#### NA statistics

**NA classification heatmap**
A heatmap for values, "technical"NA and "biological"NA. Blue: "biological"NA; Yellow: "technical"NA; Red: expression value.

<img src="./imgs/exa1_11.png">

**NA classification group-specific marker heatmap**
The specific markers of the user-interested group, data from *aftercalssification_groupmarkers.csv*. The same color mapping methods as the classification heatmap.

<img src="./imgs/exa2_12.png">

#### NA-assisted PPI analysis

**BIND ppi network**

A protein-protein interaction network. It shows the group-specific PPI found by ρBIND(**Group** has been defined), or the specific PPI found by BIND compared with Spearman ρ (**Group** has not been defined). Hover over to see the exact values.

<img src="./imgs/exa1_12.png">

<br/><br/>

### <font color= #871F78>Files for user download</font>

<img src="./imgs/download.png">

### <font color= #871F78>Files explaination</font>

#### NA statistics

| Output file | Description |
| ---       | ---         |
| **binary_mtx.csv** |  A matrix after transferring the data matrix with 1/0. |
| **overall_stat.csv** | Overall NA proportions and statistics (missing value number / total sample number).|
| **group_stat.csv** | NA proportions by groups.|

<br/>

**binary_mtx.csv**
<img src="./imgs/binary_mtx.png">

**overall_stat.csv**
<img src="./imgs/overall_stat.png">

**group_stat.csv**
<img src="./imgs/group_stat.png">

#### NA classification

| Output file | Description |
| ---         | ---         |
| **aftercalssification_label.csv** |  A data matrix with NA labeled as "biological" or "technical", for user download, and for next step analysis. |
| **aftercalssification_groupmarkers.csv** | Group-specific markers found by BIND, with the "biological"NAs replaced by 0, "technical"NAs replaced by 0.5 for convenient analysis. The expression values are 1.|
| **group_stat.csv** | NA proportions by groups.|

<br/>

**aftercalssification_label.csv**
<img src="./imgs/NAclass.png">

**aftercalssification_groupmarkers.csv**
<img src="./imgs/groupmarkers.png">

#### NA-assisted PPI analysis

If **Grp** and **Group_name** have been defined, the output is:

| Output file | Description |
| ---         | ---         |
| **BINDppi_groupdiff.csv** |  Spearman ρ(raw_raw in the file), ρBIND(weighted_rho in file) of both target group and other samples. |
| **BINDppi_groupdiff_filtered.csv** | After filtering by restricted conditions, it shows target group-specific ppi found by BIND(weighted_rho). Max 200 rows. This file is for constructing a PPI network.|

<br/>

If **Grp** and **Group_name** have not been defined, the output is:

| Output file | Description |
| ---         | ---         |
| **BINDppi_groupdiff.csv** |  Spearman ρ(raw_raw in the file), ρBIND(weighted_rho in the file) by all samples. |
| **BINDppi_groupdiff_filtered.csv** |  After filtering by restricted conditions, it shows specific PPI found by BIND compared with Spearman ρ.|

<br/>

**BINDppi_groupdiff.csv** and **BINDppi_groupdiff_filtered.csv**

<img src="./imgs/ppi_group.png">

<br/>
<br/>

**BINDppi_withoutgroup.csv** and **BINDppi_withoutgroup_filtered.csv**

<img src="./imgs/ppi_nogroup.png">
