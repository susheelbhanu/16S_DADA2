# 16S_DADA2
# About

- This repository can be used for the analyses of 16S sequencing data via `DADA2` using a `snakemake` workflow

# Set up

- You might want to adjust some settings in the files `config/sbatch.sh` and `config/config.yaml`, e.g.
the name of the `snakemake` `conda` environment, paths and number of cores/threads.

## Conda

[Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions

# create venv
conda env create -f requirements.yml -n "YourEnvName"
```

# Run the analysis

```bash
# in an interactive session
./config/sbatch.sh
# as a slurm job
sbatch ./config/sbatch.sh
```
