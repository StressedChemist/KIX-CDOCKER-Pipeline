# KIX-CDOCKER-Pipeline

This repository implements a workflow for covalent docking of small molecules into the **CBP KIX domain**, leveraging the **Covalent CDOCKER** methodology developed in [Wu & Brooks, 2022](https://pubmed.ncbi.nlm.nih.gov/35984589/).  

## What is Covalent CDOCKER?  
Covalent CDOCKER extends the CHARMM-based **Rigid CDOCKER** docking algorithm by introducing a **covalent bond grid potential** that models the free energy change of covalent bond formation between ligand and receptor. This allows the docking procedure to capture both the **initial reversible association** and the **subsequent covalent bond step**—a key requirement for accurately modeling targeted covalent inhibitors (TCIs).  

## Pipeline Overview  
1. **Protein preparation** — Clean and prepare the CBP KIX domain structure; generate interaction grids.  
2. **Ligand preparation** — Process ligands in their *pre-reaction forms* with defined covalent warheads (aldehyde, ketone, nitrile, Michael acceptor, disulfide).  
3. **Docking** — Run covalent and non-covalent docking with CDOCKER; place ligands via MD-based simulated annealing.  
4. **Scoring** — Evaluate poses using van der Waals, electrostatics, ligand internal energy, and the **covalent energy term** from the grid potential.  
5. **Pose filtering** — Identify native-like poses (RMSD ≤ 2 Å) and remove false positives.  
6. **Analysis** — Cluster poses, extract pharmacophore features, and evaluate binding interactions.  
7. **Benchmarking/Screening** — Scale to larger ligand sets for retrospective screening or lead identification across TCI chemistries.  

## References  
📖 Methodology details: [Covalent Docking in CDOCKER](https://pubmed.ncbi.nlm.nih.gov/35984589/)  
⚙️ Installation instructions: [Installation guide](https://github.com/StressedChemist/KIX-CDOCKER-Pipeline/blob/main/Installation)  

