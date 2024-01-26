# GENEFLOW: Genomic Examination and Nucleotide Evaluation For Laboratory Operations Workflow

The GENEFLOW application serves as the latest update to the Illumina sequencing processing workflow at NYU's Center for Genomics and Systems Biology. This pipeline, developed with Nextflow, encompasses a comprehensive set of procedures:

1. **Archive the Run Directory**
2. **Basecalling**
3. **Demultiplexing (Optional)**
4. **Demultiplexing Reports**
5. **Data Merging**
6. **FastQC Reports**
7. **MultiQC Report**
8. **Data Delivery**

The pipeline is designed to interface seamlessly with TuboWeb, a web-based platform for NGS data analysis and visualization, for the retrieval of metadata and customization of run parameters.

## Configuration Instructions

To successfully deploy and run the GENEFLOW pipeline, follow these setup steps:

### 1. Update `launch.sh`
Modify the launch script (`launch.sh`) to include your email address:
```bash
#SBATCH --mail-user=your_netID@nyu.edu
```

### 2. Update `nextflow.config`
Global Configuration in Nextflow
Configure the global variables in `nextflow.config` as follows:
- `alpha`: The primary work directory for the pipeline.
- `tmp_dir`: Temporary working directory for Picard tools.
- `fastqc_path`: Destination for rsyncing FastQC files (e.g., web server).
- `archive_path`: Destination for archived run directories.
- `admin_email`: Email for pipeline administration notifications.
Set up the module paths
Specify the `workDir`
Configure email settings

### 3. Update config.py
TuboWeb API Configuration
- API path
- User credentials
- API key

File Delivery and Storage Paths
- `delivery_folder_root`: Destination for FastQ files.
- `raw_run_dir_delivery_root`: Destination for raw run directories.
- `raw_run_root`: Storage location for raw run directories.
- `alpha`: The primary work directory for the pipeline (as in `nextflow.config`).

Gmail Credentials
Set the Gmail user and password for email notifications.

## Launching the Pipeline

For ease of use, a `launch.sh` script is provided to initiate the pipeline. This script requires two essential parameters and one optional parameter:

1. **Run Directory Path**
2. **Flowcell ID**
3. **Optional**: Entry point for the pipeline (used to resume the pipeline from a specific step)

### Usage Examples

- **Basic launch:**
  ```bash
  launch.sh /scratch/gencore/sequencers/{machine_name}/{run_dir_name} {fcid}

- **Specific Example Launch:**
  ```bash
  launch.sh /scratch/gencore/sequencers/NB502067/240124_NB502067_0578_AHKFT5BGXV HKFT5BGXV

- **Launch with Entry Point:**
  For resuming at a specific step like demultiplexing (e.g., 'demux'):
   ```bash
  launch.sh /scratch/gencore/sequencers/{machine_name}/{run_dir_name} {fcid} {entry}

- **Example with Entry Point:**
  ```bash
  launch.sh /scratch/gencore/sequencers/NB502067/240124_NB502067_0578_AHKFT5BGXV HKFT5BGXV demux

### Production Deployment

In a production environment, `launch.sh` is typically submitted as an SBATCH job in Slurm. Make sure the directories for error and output files are created beforehand (required by SLURM).
```
  mkdir -p /scratch/gencore/GENEFLOW/alpha/logs/HKFT5BGXV/pipeline
  sbatch --output=/scratch/gencore/GENEFLOW/alpha/logs/HKFT5BGXV/pipeline/slurm-%j.out \
         --error=/scratch/gencore/GENEFLOW/alpha/logs/HKFT5BGXV/pipeline/slurm-%j.err \
         --job-name=GENEFLOW_MANAGER_(HKFT5BGXV) \
         launch.sh /scratch/gencore/sequencers/NB502067/240124_NB502067_0578_AHKFT5BGXV HKFT5BGXV
```

