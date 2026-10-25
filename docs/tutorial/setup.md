# 1. Set up

You need three things: Java, Nextflow and one way of providing the software.

## Java and Nextflow

Nextflow needs Java 17 or later (Java 17 or 21 is the safe choice).

```bash
java -version
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/          # or anywhere on your PATH
nextflow -version           # AMRScan needs 25.04.6 or later
```

## Choose a software engine

AMRScan uses one software stack (AMRFinderPlus 4.2.7 with Python 3.12) and
can get it three ways. Pick the one that fits your machine and pass it as a
profile.

| Where you are | Profile | What you need installed |
| --- | --- | --- |
| Laptop or desktop | `docker` | Docker Desktop or Docker Engine |
| HPC cluster | `apptainer` or `singularity` | Apptainer or Singularity, usually as a module |
| No containers allowed | `conda` | Conda or Mamba |
| NCI Gadi | `gadi` | nothing extra; see part 5 |

Docker on Apple silicon runs the image under emulation. It is slower but gives
the same results.

## Get the code

```bash
git clone https://github.com/taowis/amrscan.git
cd amrscan
```

Check that Nextflow can read the pipeline:

```bash
nextflow run main.nf --help
```
