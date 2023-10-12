# Decrypting Orphan GPCR Drug Discovery via Multitask Learning
Wei-Cheng Huang1, Wei-Ting Lin1, Ming-Shiu Hung1, Jinq-Chyi Lee1, Chun-Wei Tung1*
1 Institute of Biotechnology and Pharmaceutical Research, National Health Research Institutes, Miaoli County, Taiwan
* Correspondence: cwtung@nhri.edu.tw 

## Abstract
The drug discovery of G protein-coupled receptors (GPCRs) superfamily using computational models is often limited by the availability of protein three-dimensional (3D) structures and chemicals with experimentally measured bioactivities. Orphan GPCRs without known ligands further complicate the process. To enable drug discovery for human orphan GPCRs, multitask models were proposed for predicting half maximal effective concentrations (EC50) of the pairs of chemicals and GPCRs. Protein multiple sequence alignment features, and physicochemical properties and fingerprints of chemicals were utilized to encode the protein and chemical information, respectively. The protein features enabled the transfer of data-rich GPCRs to orphan receptors and transferability based on similarity of protein features. The final model was trained using both agonist and antagonist data from 200 GPCRs and showed an excellent mean squared error (MSE) of 0.24 in the validation dataset. An independent test using the orphan dataset consisting of 16 receptors associated with less than 8 bioactivities showed a reasonably good MSE of 1.51 that can be further improved to 0.53 by considering transferability based on protein features. The informative features were identified and mapped to corresponding 3D structures to gain insights into the mechanism of GPCR-ligand interactions across the GPCR family. The proposed method provides a novel perspective on learning ligand bioactivity within the diverse human GPCR superfamily and can potentially accelerate the discovery of therapeutic agents for orphan GPCRs.

## Keywords
Multitask Learning; G Protein-Coupled Receptors; GPCR; Feature Selection; Ligand-Based Virtual Screening 

![image]('Figure 1.jpg')



# Installation and envirment setup
You may download whole directory with script files and compressed result files.

1. use the python_env_setup.sh file to setup all dependencies
2. use the extract_files.sh file to extract all the result files for this task.

```
git clone git@github.com:drhuangwc/GPCR.git
sh python_env_setup.sh
sh extract_files.sh
```

Script files:
1. `GPCR_MTL_dataprocess.ipynb`: scripts + original dataset
2. `GPCR_MTL_autogluontrain.ipynb`: scripts for separate protein Agonist using Autogluon + separate protein Antagonist using Autogluon + separate protein Agonist & Antagonist using Autogluon
3. `GPCR_MTL_FeatureSelection.ipynb`: scripts for using mRMR feature selection of the files + FeatureSelection_GPCRtrain + result data
4. `GPCR_MTL_TRAINevaluate.ipynb`: scripts for evaluation of validation data and orphan data + result data
5. `GPCR_MTL_similarity.ipynb`: scripts for calculate Tanimoto similarity with or without mRMR selected features in 7 parts of the GPCR
6. `GPCR_MTL_colorPDB.ipynb`: scripts for mapping the mRMR selected features into the protein structure files

