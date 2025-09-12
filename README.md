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

## 0. Installing needed tools for CDOCKER

In order to use **CDOCKER** you will need to:

- Create a **conda environment** capable of building **CHARMM** (and optionally **pyCHARMM**).
- (Recommended for the paper's workflow) Install the **MMTSB ToolSet** for pose clustering; follow the MMTSB installation instructions.
- Obtain **CHARMM** (free to academics and government labs) from Academic CHARMM. Follow the directions below to build a conda environment capable of installing **CHARMM/pyCHARMM**.
- Install **CHARMM** and (optionally) **pyCHARMM**.

Reference for what the paper actually used: Wu & Brooks (2022), *Covalent Docking in CDOCKER*, JCAMD.

### Create a conda environment

You will need a base Anaconda/Miniconda installation.

Make a conda environment (you can also use `conda env create -f <env>.yml` in §2):

```bash
conda create -y -n <name_of_environment> python=3.9  # Python only needed if you'll use pyCHARMM
```

Activate this environment:

```bash
conda activate <name_of_environment>
```

**Install mamba** as a faster conda:

```bash
conda install -y -c conda-forge mamba
```

**(Optional, GPU)** Install CUDA from NVIDIA for OpenMM GPU acceleration

```bash
# Common choices—pick ONE that matches your system driver (see table below)
mamba install -y -c nvidia cuda                                   # typically installs CUDA 12.1.1
mamba install -y -c "nvidia/label/cuda-12.0.0" cuda                # installs CUDA 12.0
```

You can see available CUDA Toolkit packages on the NVIDIA channels.

(Driver/Toolkit compatibility table below mirrors the example format.)

### Install needed packages to build CHARMM (and match the paper's CDOCKER workflow)

**Minimal + paper-faithful** (no MPI, no FFTDOCK/BLaDE flags):

```bash
mamba install -y -c conda-forge \
  gcc gxx gfortran make cmake binutils \
  openmm rdkit
```

- **gcc/gfortran/make/binutils/cmake** — toolchain to compile CHARMM.
- **OpenMM** — enables CHARMM/OpenMM parallel simulated annealing used in the paper.
- **RDKit** — ligand conformer generation (paper used ETKDG).

You'll obtain **CGenFF/ParamChem** externally (not on conda) and install **MMTSB ToolSet** separately if you want the same clustering step.

### Note on CUDA Toolkit/Driver and Compiler Compatibilities

Check the installed **CUDA driver** on your GPU node with:

```bash
nvidia-smi
```

Example (top of the output):

```
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 525.85.05  Driver Version: 525.85.05  CUDA Version: 12.0        |
|-------------------------------+----------------------+----------------------+
```

A driver **525.85.05** is compatible with **CUDA 12.0**, **GCC ≥ 12.1**, or **Intel Compilers 2021.6**.

| Toolkit Version | Minimum Required Driver | Recommended GCC | Recommended Intel Compilers |
| --- | --- | --- | --- |
| CUDA 12.1.x | ≥ 530.30.02 | 12.2 | 2021.7 |
| CUDA 12.0.0 | ≥ 525.85.05 | 12.1 | 2021.6 |
| CUDA 11.8.x | ≥ 520.61.05 | 11 | 2021 |
| CUDA 11.7.x | ≥ 515.48.07 | 11 | 2021 |
| CUDA 11.6.x | ≥ 510.47.03 | 11 | 2021 |
| CUDA 11.5.x | ≥ 495.29.05 | 11 | 2021 |
| CUDA 11.4.x | ≥ 470.82.01 | 9.x | 19.1 |
| CUDA 11.3.x | ≥ 465.19.01 | 9.x | 19.1 |
| CUDA 11.2.x | ≥ 460.32.03 | 9.x | 19.1 |
| CUDA 11.1 | ≥ 455.32 | 9.x | 19.1 |
| CUDA 11.0 | ≥ 450.51.06 | 9.x | 19.1 |
| CUDA 10 | ≥ 440.33 | 10.2 | 18.0 |
| CUDA 9 | ≥ 396.37 | 4.8.5 | 17.0 |
| CUDA 8 | ≥ 375.26 | 4.8.2 | 15, 16 |

## 2. Building the CHARMM/pyCHARMM compatible environment with a YAML file

**cdocker_wcuda12.yml**

```yaml
name: cdocker_wcuda12      # Name of your conda environment
channels:
  - conda-forge
  - nvidia/label/cuda-12.0.0 # pin CUDA 12.0 exactly; change as needed
  # nvidia/label/cuda-12.1.1 # alt CUDA line (comment/uncomment as appropriate)
dependencies:
  - python=3.9 # only needed if you'll use pyCHARMM
  - mamba
  - cuda # brings in the CUDA toolkit for GPU OpenMM
  - ca-certificates
  - certifi
  - openssl
  # Toolchain
  - gcc
  - gxx
  - gfortran
  - make
  - cmake
  - binutils
  # Paper-faithful runtime
  - openmm # CHARMM/OpenMM parallel SA during docking
  - rdkit # conformer generation
  # (Optional analytics/vis—commented out to keep this lean)
  # - pandas
  # - mdtraj
  # - biopython
  # - py3dmol
  # - pymol-open-source
  # - jupyterlab
  # - jupytext
  # prefix: /path/to/.conda/envs/cdocker_wcuda12 # optional explicit path
```

**Install from the YAML:**

```bash
conda env create -f cdocker_wcuda12.yml
```

(You can edit this YAML to change CUDA versions by switching the NVIDIA channel line, as shown above.)

## 3. CHARMM and pyCHARMM installation

*once conda environment is installed and active*

### Build CHARMM (CDOCKER-ready)

**Minimal CPU build** (no MPI, no FFTDock/BLaDE):

```bash
conda activate cdocker_wcuda12
cd <charmm_root>
mkdir -p build_charmm && cd build_charmm
# If OpenMM is installed in this conda env, detection is usually automatic.
# If needed, you can help the build find it:
export OPENMM_HOME="$CONDA_PREFIX"
# Minimal CPU build for CDOCKER (no FFTDOCK/BLaDE/MPI)
../configure -u -p <charmm_install_path>
make -j <n> install
```

**OpenMM-enabled build** (recommended to match Wu & Brooks, 2022):

```bash
conda activate cdocker_wcuda12
cd <charmm_root>/build_charmm
# Help configure find OpenMM from conda if necessary
export OPENMM_HOME="$CONDA_PREFIX"
# Build with OpenMM support, still without FFTDock/BLaDE/MPI
../configure --with-openmm -u -p <charmm_install_path>
make -j <n> install
```

- `<charmm_root>` is the path to the CHARMM source tree.
- `<charmm_install_path>` is where you want the CHARMM installation to live.
- `<n>` is the number of build threads.

### Why no `--with-fftdock` or `--with-blade` here?

Those are separate CHARMM features not required for **Rigid CDOCKER**; the paper did *not* rely on them. The recommended path is to enable **OpenMM** so you can use CHARMM/OpenMM **parallel simulated annealing** during docking.

### Build pyCHARMM (optional)

```bash
conda activate cdocker_wcuda12
cd <charmm_root>/build_charmm
rm -rf *   # Clean the build directory
# As a shared library (for pyCHARMM); keep it lean and consistent with above
export OPENMM_HOME="$CONDA_PREFIX"
../configure --as-library --without-mpi --with-openmm -u -p <pycharmm_install_path>
make -j <n> install
# Install the Python bindings
cd <charmm_root>
pip install "$(pwd)/tool/pycharmm"
# Point pyCHARMM at your CHARMM lib
export CHARMM_LIB_DIR="<pycharmm_install_path>/lib"
```

`<pycharmm_install_path>` is where you want the pyCHARMM-ready CHARMM library to reside.

### Final notes

- **MMTSB ToolSet** (for pose clustering) and **CGenFF/ParamChem** (for ligand topology/parameters) are outside this conda env; install/obtain them separately to replicate the paper's workflow.
- If you only run CDOCKER from CHARMM scripts (no Python wrapper), you can omit Python/pyCHARMM entirely.
- If you are CPU-only, skip the CUDA bits; OpenMM still works on CPU.
