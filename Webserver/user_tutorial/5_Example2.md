## 🧫Example 2

This example is an extracellular vesicle (EV) dataset in human cancer and normal tissue derived from a variety of solid tumor tissues. BIND analyzed the missing value composition of this dataset and enabled NA classification to identify cancer-specific protein markers. In addition, the interaction of cancer-specific proteins was identified, which provided a basis for related studies.

### <font color= #871F78>Step 1: Loading data</font>

Click ***Run HumanEV Demo***

<img src="./imgs/exa2_1.png">


### <font color= #871F78>Step 2: Job submission</font>

View Example data parameters / download example data files.

Click ***Run*** to run sample data.

<img src="./imgs/exa2_2.png">

### <font color= #871F78>Step 3: Job Query</font>

Successfully submitted tasks will be given a UUID, please copy or download the UUID for the task query.

<img src="./imgs/exa1_3.png">
Click ***Task Query*** on the top to query the submitted tasks.

<img src="./imgs/taskquery.png">

Paste the UUID of your task. Click ***Search***
<img src="./imgs/pasteuuid.png">

Get task information. If the task status is "running" or "pending", please wait until the task is completed. If the task status is "failed", please contact us.

<img src="./imgs/exa2_3.png">

Click ***View Task Detail*** after the task completed

<img src="./imgs/exa2_4.png">

### <font color= #871F78>Step 4: File download</font>

You can Download all output files on the Download Output Files page. See the Output explanation page for an explanation of the file.
<img src="./imgs/exa1_6.png">

### <font color= #871F78>Step 5: Visualization results</font>

**NA heatmap**
NA Heatmap is the missing value heatmap that visualizes the overall pattern of missing values. Red: values; Blue: NAs. The user can observe the overall pattern of missing values in the dataset.

<img src="./imgs/exa2_5.png">

**NA proportion histogram**
NA proportion histogram represents the distribution of NA proportions of different proteins in the dataset. Unlike Example1, the number of proteins in this dataset with NA ratios close to 1 is much larger than in the other cases.
X axis: NA proportion
Y axis: Protein counts

<img src="./imgs/exa2_6.png">

**NA proportion pie chart**
NA proportion pie chart shows the pie chart of overall NA proportions (0 / 0-0.2 / 0.2-0.8 / 0.8-1). The percentage of proteins with NA proportions of 0.8-1 was significantly higher in this dataset.

<img src="./imgs/exa2_7.png">

**NA proportion scatter plot**
NA proportion scatter plot is an x-y scatter plot of NA proportions vs mean expression values of a protein. Linear regression results are provided. In this dataset, the proportion of NA was also negatively correlated with the average protein expression.
X axis: NA proportions in the dataset.
Y axis: Mean expression value of the proteins.
<img src="./imgs/exa2_9.png">

**Expression value with NA proportion label point plot**
Given a specific sample, ranking the protein expression values, label the overall NA proportion > 0.8/ <0.2/ 0.2-0.8 proteins by color.
X axis: proteins in the sample
Y axis: protein expression value
Red point: NA proportion > 0.8
Blue point: NA proportion < 0.2
<img src="./imgs/exa2_13.png">

**Group-specific NA proportion scatter plot**
x-y scatter plot of NA proportions in a specific group vs NA proportions in other groups.
X axis: NA proportions in the query group
Y axis: NA proportions in other samples

<img src="./imgs/exa2_14.png">

**Group level NA proportion heatmap**
Heatmap of NA proportions in different groups, clustering by samples. 
<img src="./imgs/exa2_15.png">

**NA classification heatmap**
Heatmap for values, "technical"NA and "biological"NA. Blue: "biological"NA; Yellow: "technical"NA; Red: expression value.

<img src="./imgs/exa2_10.png">

**NA classification group-specific marker heatmap**
Group markers heatmap is the specific markers of the haematopoietic and lymphoid cancer group. The same color mapping methods as the classification heatmap. BIND identifies protein markers specific to EV in multiple cancer tissues.

<img src="./imgs/exa2_12.png">


**BIND ppi network**
PPI network is a protein-protein interaction network. It shows the cancerEV-specific PPI found by ρBIND.

<img src="./imgs/exa2_11.png">