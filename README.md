# LearningPanCancerGene
A novel unsupervised learning model for detecting driver genes from pan-cancer data through matrix tri-factorization framework with pairwise similarity constraints
![image](https://github.com/JianingXi/LearningPanCancerGene/blob/master/bin/splash.jpg)

Developer: Jianing Xi <xjn@mail.ustc.edu.cn> from Health Informatics Lab, School of Information Science and Technology, University of Science and Technology of China

## Instructions to LearningPanCancerGene (version 1.0.0)

# Requirement
* 4GB memory
* MATLAB R2015a or later

# Step-by-Step Tutorial

### Step 1. Input data

* Somatic mutations of patients across multiple cancer types
The file "./MutationData/mutation_matrices.mat" contains the TCGA somatic mutation data matrix and the bipartite matrix of relationships between patients and cancer types, which is downloaded from a [previous study](https://academic.oup.com/bioinformatics/article/32/11/1643/1742725/An-integrative-somatic-mutation-analysis-to) [1]. If you want to use a user-defined dataset, you can edit the structure variables `m_cdata` in file "./MutationData/mutation_matrices.mat":

  * `m_cdata.num_rows` <- the number of patients in user-defined dataset.}
  * `m_cdata.num_cols` <- the number of genes in the user-defined dataset.
  * `m_cdata.num_class` <- the number of cancer types in the user-defined dataset.
  * `m_cdata.patientID` <- the patient IDs for all the patients in the user-defined dataset.
  * `m_cdata.classID` <- the cancer type classes for all the patients in the user-defined dataset.
  * `m_cdata.geneID` <- the gene symbols for all the genes in the user-defined dataset.
  * `m_cdata.className` <- the names of the cancer types for all classes in the user-defined dataset.
  * `m_cdata.X` <- the mutation matrix (patient x gene) of the user-defined dataset.
  * `m_cdata.F` <- the binary bipartite matrix of the relationships between patients and cancer types (patient x gene) of the user-defined dataset.


* Disease similarities
The file "./DOSim/CancerSimilarity.mat" contains the similarites between the cancer types included in the data aforementioned, which is calculated based on disease ontologies [2] by a previous published tool [3]. If you want to use a user-defined disease similarity, you can edit the similarity matrix `CancerSimilarity` in file "./DOSim/CancerSimilarity.mat". The entries of the similarity matrix is constructed as follow:
  * `CancerSimilarity(i,j)` <- the similarity score of the i-th and j-th cancers obtained by DOSim tool [3], if i != j
  * `CancerSimilarity(i,j)` <- 1, if i == j
  * The indices of cancer types of matrix `CancerSimilarity` should be the same indices of cancer types of the relationships matrix `m_cdata.F`.

* Gene interaction network
The files "./network/edge_list.txt" and "./network/index_genes.txt" contain the edges and gene nodes of interaction network [iRefIndex 9](http://irefindex.org) [4]. If you want to use a user-defined network, you can replace the two default files with the files of the user-defined network:

  * Replace the file "./network/index_genes.txt" by a text file contains two columns. The first column represents the indice of the genes, and the second column represents the related gene symbols.
  * Replace the file "./network/edge_list.txt" by a text file contains three columns. The first column represents the indice of the source nodes of edges, the second column represents the target nodes of edges, and the third column represents their related affinity scores.


### Step 2. Parameter settings
* The three tuning parameter `lambda_S`, `lambda_L`, `lambda_V` can be reset in line 54-57 in the Matlab script file "./demo.m". The tuning parameter `lambda_S` is used to control the derivation of the similarity scores between cancer types. The tuning parameter `lambda_V` is used to balance the fitness of the model and the regularization term of sparsity on representation matrix V. The tuning parameter `lambda_L` is used to control the closeness between the representation vectors of the interacted genes.

  ![image](https://github.com/JianingXi/LearningPanCancerGene/blob/master/bin/formula.PNG)


  * The default values of three tuning parameter `lambda_V`, `lambda_S`, `lambda_L` are 1, 1 and 0.01 respectively.
  * The parameters can be reset by the user for their own task.

### Step 3. Run LearningPanCancerGene
* When the user-defined dataset and parameters are set, you can apply LearningPanCancerGene by running the Matlab script file "./demo.m". The results will then be automatically saved in file "./Output/Results.mat" after the program is finished. 

### Step 4. Variables in Output file 
* In file "./Output/Results.mat", the variable `Candidates_list` is driver gene candidates selected from the top 200 genes ranked by representation `V`. The matrix `V` is the representation matrix of the investigated genes. The variables `S` is the output similarity matrix of cancer types.

# References
[1] Park, Sunho and Kim, Seung-Jun and Yu, Donghyeon and Pena-Llopis, Samuel and Gao, Jianjiong and Park, Jin Suk and Chen, Beibei and Norris, Jessie and Wang, Xinlei and Chen, Min and others, An integrative somatic mutation analysis to identify pathways linked with survival outcomes across 19 cancer types, Bioinformatics 32 (11) (2015) 1643-1651.

[2] L. M. Schriml, C. Arze, S. Nadendla, Y.-W. W. Chang, M. Mazaitis, V. Felix, G. Feng, W. A. Kibbe, Disease Ontology: a backbone for disease semantic integration, Nucleic acids research 40 (D1) (2011) D940-D946.

[3] G. Yu, L.-G. Wang, G.-R. Yan, Q.-Y. He, DOSE: an R/Bioconductor package for disease ontology semantic and enrichment analysis, Bioinformatics 31 (4) (2014) 608-609.

[4] S. Razick, G. Magklaras, I. M. Donaldson, iRefIndex: a consolidated protein interaction database with provenance, BMC bioinformatics 9 (1) (2008) 1.
