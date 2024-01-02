#!/usr/bin/sh


###### FOLDER and subfolders: colorPDB_selefeatureTOP200

cd colorPDB_selefeatureTOP200
tar -xvf GPCR_to_trainAgonistAntagonist_humanalltrainconcat_selefeatureTOP200.tar.xz

tar -xvf OriPDB.tar.xz

###### FOLDER and subfolders: GPCR_agonist_GPRC5A

cd ../GPCR_agonist_GPRC5A
cat AutogluonTRAIN0618_AgonistAntagonist_humanalltrainconcat_selefeatureTOP200_102min.tar.xz.parta* > AutogluonTRAIN0618_AgonistAntagonist_humanalltrainconcat_selefeatureTOP200_102min.tar.xz 
tar -xvf AutogluonTRAIN0618_AgonistAntagonist_humanalltrainconcat_selefeatureTOP200_102min.tar.xz


###### FOLDER and subfolders: Human_AgonistAntagonist

cd ../Human_AgonistAntagonist
tar -xvf GPCR_AgonistAntagonist20220915.tar.xz
tar -xvf GPCR_AgonistAntagonist_protein20220915.tar.xz
tar -xvf GPCR_AgonistAntagonist_protein_align_ptcode123_20220915.tar.xz

cat GPCR_AgonistAntagonist_protein_align_ptcode123_padel1d2d20220921_rdkECFP1019.tar.xz.part* > GPCR_AgonistAntagonist_protein_align_ptcode123_padel1d2d20220921_rdkECFP1019.tar.xz 
tar -xvf GPCR_AgonistAntagonist_protein_align_ptcode123_padel1d2d20220921_rdkECFP1019.tar.xz

cat GPCR_to_trainAgonistAntagonist_all.tar.xz.part* > GPCR_to_trainAgonistAntagonist_all.tar.xz 
tar -xvf GPCR_to_trainAgonistAntagonist_all.tar.xz

cat padleout_molecule_1d2ddescriptor09190015_GPCR_AgonistAntagonist_protein_align_ptcode123_20220915.tar.xz.part* > padleout_molecule_1d2ddescriptor09190015_GPCR_AgonistAntagonist_protein_align_ptcode123_20220915.tar.xz
tar -xvf padleout_molecule_1d2ddescriptor09190015_GPCR_AgonistAntagonist_protein_align_ptcode123_20220915.tar.xz

###### FOLDER and subfolders: Human_AgonistAntagonist/GPCR_similarty

cd GPCR_similarty
tar -xvf GPCR_to_trainAgonistAntagonist_humanalltrainconcat_selefeatureTOP200.tar.xz
tar -xvf similarty_proteinpairsimilarity20230624_TOP200_proteinseq20220915_align_7CUT.tar.xz


###### FOLDER and subfolders: Human_AgonistAntagonist/GPCR_to_test_Mol2Vec

cd ../GPCR_to_test_Mol2Vec
cat seperatePROTEIN_Mol2Vec.tar.xz.part* > seperatePROTEIN_Mol2Vec.tar.xz
tar -xvf seperatePROTEIN_Mol2Vec.tar.xz

###### FOLDER and subfolders: Human_AgonistAntagonist/GPCR_to_test_HUMANdeletion95valid

cd ../GPCR_to_test_HUMANdeletion95valid
cat AutogluonTRAIN0520_AgonistAntagonist_humanalltrainconcat_420min.tar.xz.part* > AutogluonTRAIN0520_AgonistAntagonist_humanalltrainconcat_420min.tar.xz
tar -xvf AutogluonTRAIN0520_AgonistAntagonist_humanalltrainconcat_420min.tar.xz

###### FOLDER and subfolders: Human_AgonistAntagonist/GPCR_to_test_HUMANmutant95valid

cd ../GPCR_to_test_HUMANmutant95valid
cat AutogluonTRAIN0520_AgonistAntagonist_humanalltrainconcat_420min.tar.xz.part* > AutogluonTRAIN0520_AgonistAntagonist_humanalltrainconcat_420min.tar.xz
tar -xvf AutogluonTRAIN0520_AgonistAntagonist_humanalltrainconcat_420min.tar.xz


###### FOLDER and subfolders: Human_AgonistAntagonist/GPCR_separatePROTEIN


cd ../GPCR_separatePROTEIN
cat Agonist.tar.xz.part* > Agonist.tar.xz
tar -xvf Agonist.tar.xz

cat AgonistAntagonist.tar.xz.part* > AgonistAntagonist.tar.xz
tar -xvf AgonistAntagonist.tar.xz

cat Antagonist.tar.xz.part* > Antagonist.tar.xz
tar -xvf Antagonist.tar.xz


###### FOLDER and subfolders: Human_AgonistAntagonist/GPCR_to_test_FeatureSelection
cd ../GPCR_to_test_FeatureSelection_tarxz_files
cat GPCR_to_test_FeatureSelection.tar.xz.part* GPCR_to_test_FeatureSelection.tar.xz
tar -xvf GPCR_to_test_FeatureSelection.tar.xz

echo "DONE all files extracted."


