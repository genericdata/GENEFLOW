# GENEFLOW: Genomic Examination and Nucleotide Evaluation For Laboratory Operations Workflow

The GENEFLOW application serves as the latest update to the Illumina sequencing processing workflow at NYU's Center for Genomics and Systems Biology. This pipeline, developed with Nextflow, encompasses a comprehensive set of procedures:

1. **Archive the Run Directory**
   - Utilize tools like tar, bzip, and rsync to archive the run directory.

2. **Basecalling**
   - Basecall the data using Picard Tools `IlluminaBasecallsToFastq`.

3. **Demultiplexing (Optional)**
   - Optionally demultiplex using PHENIQS.

4. **Demultiplexing Reports**
   - Generate detailed reports post-demultiplexing.

5. **Data Merging**
   - Merge data as necessary (e.g., for NextSeq runs with a single pool across multiple lanes).

6. **FastQC Reports**
   - Create FastQC reports for quality assessment.

7. **MultiQC Report**
   - Compile a comprehensive MultiQC report.

8. **Data Delivery**
   - Efficiently deliver data to end-users.


The pipeline is designed to interface seamlessly with TuboWeb for the retrieval of metadata and customization of run parameters.

## Configuration Instructions

To successfully deploy and run the GENEFLOW pipeline, follow these setup steps:

### 1. Email Configuration in Launch Script
Modify the launch script (`launch.sh`) to include your email address:
```bash
#SBATCH --mail-user=your_netID@nyu.edu
```

### 2. Global Configuration in Nextflow
Configure the global variables in `nextflow.config` as follows:
- `alpha`: The primary work directory for the pipeline.
- `tmp_dir`: Temporary working directory for Picard tools.
- `fastqc_path`: Destination for rsyncing FastQC files (e.g., web server).
- `archive_path`: Destination for archived run directories.
- `admin_email`: Email for pipeline administration notifications.

### 3. Module Paths
Set up the module paths in `nextflow.config`.

### 4. Nextflow Work Directory
Specify the `workDir` in `nextflow.config`.

### 5. Email Settings in Nextflow
Configure email settings within `nextflow.config`.

### 6. TuboWeb API Configuration
In `config.py`, configure the TuboWeb API settings:
- API path
- User credentials
- API key

### 7. File Delivery and Storage Paths
Define the following paths in `config.py`:
- `delivery_folder_root`: Destination for FastQ files.
- `raw_run_dir_delivery_root`: Destination for raw run directories.
- `raw_run_root`: Storage location for raw run directories.
- `alpha`: The primary work directory for the pipeline (as in `nextflow.config`).

### 8. Gmail Credentials
Set the Gmail user and password in `config.py` for email notifications.

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

