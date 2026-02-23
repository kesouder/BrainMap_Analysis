# DSC Capstone Quarter 2 Project: Structural vs. Chemical: A Statistical Approach to the Multimodal Correlation Structure of the Human Brain

### Co-Authors: Kevin Huang, Kevin Souder

Our project utilizes the **[neuromaps](https://www.nature.com/articles/s41592-022-01625-w)** framework — a standardized system for comparing human brain maps across modalities, spaces, and scales.  
Our capstone project uses **neuromaps** to evaluate the relationship and interaction between the structural and chemical & functional layers of the brain. To combat spatial autocorrelation, we use the **spin test** for non-parametric spatial null modeling to ensure our results are accurate.

### Our Code

Our workflow is in 'main.ipynb' for readability. <br> This implements our entire project.
<br>
'index.html' contains the code and styling for our website, along with the images folder. The actual website is **[here](https://kevinhuang8706.github.io/dsc_capstone_website/)**. 
The rest of the files in this repo were downloaded during the installation of neuromaps. See below how to download neuromaps toolbox.

## Purpose of this Repository

This repository provides reproducible code for:

### 1. Data Standardization (Neuromaps)
* **Fetching**: Automated retrieval of publicly available maps, including 5-HT1A receptors, CMRglc metabolism, PC1 gene expression, and myelin (T1w/T2w).
* **Cross-Space Transformation**: Mapping volumetric MNI152 data and surface-based CIVET/fsaverage data onto a common **fsLR 32k mesh** using Connectome Workbench to ensure vertex-wise correspondence.

### 2. Advanced Exploratory Data Analysis (EDA)
* **Spatial Visualizations**: Generation of cortical surface projections for each modality to inspect spatial distributions and identify high/low intensity regions.
* **Distribution Profiles**: Quantitative evaluation of value distributions across the cortex to inform the use of non-parametric, rank-based statistics.

### 3. Multimodal Correlation & Spatial Statistics
* **Pairwise Spearman Correlation**: Construction of an $8 \times 8$ multimodal correlation matrix capturing the relationships across diverse biological layers.
* **Spatial Spin Tests**: Addressing spatial autocorrelation—the "inflated significance" problem in neuroimaging—using the Alexander-Bloch spin permutation method (1,000 permutations).
* **Max-T Correction**: Implementing a spin-based max-T correction procedure to strictly control the family-wise error rate (FWER) across the correlation matrix.



### 4. Dimensionality Reduction & Clustering
* **Principal Component Analysis (PCA)**: Identifying the dominant axes of variance (gradients) that characterize the brain’s multimodal organization and identifying the maps with the highest loadings.
* **Hierarchical Clustering**: Utilizing Ward’s linkage to identify natural biological groupings (subgroups) among structural, chemical, and functional maps.
* **Subgroup Connectivity Analysis**: Quantifying the internal coherence of hypothesized structural and neurochemical modules and evaluating cross-modal interaction.

## Datasets Used

The following eight brain maps were fetched and standardized for this analysis:

1.  **5-HT1A Receptor** (Savli et al., 2012) - _Neurochemical_
2.  **Functional Gradient 1** (Margulies et al., 2016) - _Functional_
3.  **Glucose Metabolism (CMRglc)** (Vaishnavi et al., 2010) - _Metabolic_
4.  **Intersubject Variability** (Mueller et al., 2013) - _Functional_
5.  **PC1 Gene Expression** (Burt et al., 2018 / Markello et al., 2021) - _Genetic_
6.  **Myelin (T1w/T2w Ratio)** (Glasser et al., 2016) - _Structural_
7.  **Cortical Thickness** (Glasser et al., 2016) - _Structural_
8.  **Allometric Scaling (NIH)** (Reardon et al., 2018) - _Structural_

## About Neuromaps

**neuromaps** (Markello & Mišić, _Nature Methods_, 2022) is an open-source Python package that standardizes the comparison of brain maps. It provides:

- **Public access** to over 40 published brain maps (e.g., myelin, metabolism, cortical thickness).
- **Standardized coordinate frameworks**, allowing easy transformation across different coordinate systems (fsLR, MNI152, CIVET).
- **Spatial statistics**, various methods for generating spatial nulls for both surface-based and volumemetric-based maps
- **Easy comparison and manipulation**, direct plotting, array-handling, and testing functions applicable for all brain maps.

**Dependencies** For python >= 3.8+, it's possible to directly use `pip` to download neuromaps using `pip install neuromaps`. But it's recommended to install it from
the source directly in order to fetch the latest version and **other required packages**.

```
git clone https://github.com/netneurolab/neuromaps
cd neuromaps
pip install .
```

For plotting, we used `matplotlib`, you can install it in Python using the command

```
pip install matplotlib
```

**IMPORTANT** In order to use tranformations and many other crucial functionalities in `neuromaps`, you must have **[Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench)** ready in your path, you can check this by the command

```
wb_command -version
```

**Data Location**
All datasets are automatically fetched and cached using the neuromaps.datasets module.
By default, data are stored at:

- **macOS/Linux**: `~/neuromaps-data/`

- **Windows**: `C:\Users\<you>\neuromaps-data\`

## References

1. **Markello, R. D., & Mišić, B. (2022).**  
   _Neuromaps: structural and functional maps of the human brain._  
   **Nature Methods, 19**, 1117–1124.  
   [https://doi.org/10.1038/s41592-022-01625-w](https://doi.org/10.1038/s41592-022-01625-w)

2. **Alexander-Bloch, A., Shou, H., Liu, S., Satterthwaite, T. D., Glahn, D. C., Shinohara, R. T., Vandekar, S. N., & Raznahan, A. (2018).**  
   _On testing for spatial correspondence between maps of human brain structure and function._  
   **NeuroImage, 178**, 540–551.  
   [https://doi.org/10.1016/j.neuroimage.2018.05.070](https://doi.org/10.1016/j.neuroimage.2018.05.070)

3. **Vaishnavi, S. N., Vlassenko, A. G., Rundle, M. M., Snyder, A. Z., Mintun, M. A., & Raichle, M. E. (2010).**  
   _Regional aerobic glycolysis in the human brain._  
   **Proceedings of the National Academy of Sciences (PNAS), 107**(41), 17757–17762.  
   [https://doi.org/10.1073/pnas.1010459107](https://doi.org/10.1073/pnas.1010459107)  
   _Covers CBV, CBF, CMRglc, and CMRO₂ metabolic maps used as target datasets._

4. **Glasser, M. F., Sotiropoulos, S. N., Wilson, J. A., Coalson, T. S., Fischl, B., Andersson, J. L., Xu, J., Jbabdi, S., Webster, M., Polimeni, J. R., Van Essen, D. C., & Jenkinson, M. (2013).**  
   _The minimal preprocessing pipelines for the Human Connectome Project._  
   **NeuroImage, 80**, 105–124.  
   [https://doi.org/10.1016/j.neuroimage.2013.04.127](https://doi.org/10.1016/j.neuroimage.2013.04.127)  
   _Reference for HCP S1200-derived myelin and cortical thickness maps._

5. **Margulies, D. S., Ghosh, S. S., Goulas, A., Falkiewicz, M., Huntenburg, J. M., Langs, G., Bezgin, G., Eickhoff, S. B., Castellanos, F. X., Petrides, M., Jefferies, E., & Smallwood, J. (2016).**  
   _Situating the default-mode network along a principal gradient of macroscale cortical organization._  
   **Proceedings of the National Academy of Sciences (PNAS), 113**(44), 12574–12579.  
   [https://doi.org/10.1073/pnas.1608282113](https://doi.org/10.1073/pnas.1608282113)  
   _Provides functional gradient maps 1–3 used in the comparison._

6. **Markello, R. D., Spreng, R. N., Luh, W. M., Anderson, K. M., Ge, T., Holmes, A. J., & Mišić, B. (2021).**  
   _Deconstructing the human connectome: Sensory, transmodal, and association networks are organized along distinct gradients of gene expression._  
   **Nature Neuroscience, 24**, 1255–1265.  
   [https://doi.org/10.1038/s41593-021-00815-8](https://doi.org/10.1038/s41593-021-00815-8)  
   _Source for the Allen Human Brain Atlas gene-expression-derived gradients._

7. **Sydnor, V. J., Larsen, B., Bassett, D. S., Alexander-Bloch, A., Fair, D. A., Liston, C., Mackey, A. P., Milham, M. P., Pines, A., Roalf, D. R., Seidlitz, J., Xu, T., Raznahan, A., Satterthwaite, T. D., & Graham, A. M. (2021).**  
   _Neurodevelopment of the association cortices: Patterns, mechanisms, and implications for psychopathology._  
   **Neuron, 109**(18), 2820–2846.e7.  
   [https://doi.org/10.1016/j.neuron.2021.06.016](https://doi.org/10.1016/j.neuron.2021.06.016)  
   _Associated with the SAaxis map of cortical surface area._

8. **Xu, T., Dong, H.-M., Zhang, S., Zuo, X.-N., & Milham, M. P. (2020).**  
   _Cross-species functional alignment reveals evolutionary hierarchy within the connectome._  
   **NeuroImage, 223**, 117346.  
   [https://doi.org/10.1016/j.neuroimage.2020.117346](https://doi.org/10.1016/j.neuroimage.2020.117346)  
   _Provides Functional Homology (FChomology) and Evolutionary Homology (evoexp) maps._

9. **Reardon, P. K., Seidlitz, J., Vandekar, S., Liu, S., Patel, R., Park, M. T. M., Alexander-Bloch, A., Clasen, L. S., Blumenthal, J. D., Lalonde, F. M., & others. (2018).**  
   _Normative brain size variation and brain shape diversity in humans._  
   **Science, 360**(6394), 1222–1227.  
   [https://doi.org/10.1126/science.aar2578](https://doi.org/10.1126/science.aar2578)  
   _Corresponds to the Allometric Scaling (PNC) and Allometric Scaling (NIH) dataset in CIVET 41k space._
10. **Mueller, S., Wang, D., Fox, M. D., Yeo, B. T. T., Sepulcre, J., Sabuncu, M. R., Shafee, R., Lu, J., & Liu, H. (2013).**  
    _Individual variability in functional connectivity architecture of the human brain._  
    **Neuron, 77**(3), 586–595.  
    [https://doi.org/10.1016/j.neuron.2012.12.028](https://doi.org/10.1016/j.neuron.2012.12.028)  
    _Provides the intersubject functional variability map._
