# SRD5A2 InSilicoAssay

# SRD5A2 InSilicoAssay

**Description:**

SRD5A2 InSilicoAssay is a computational tool designed for the identification of inhibitors of the SRD5A2 enzyme through virtual screening. This README provides an overview of how to use the tool, including setup, input requirements, execution, and additional background on the target's physiology.

**Target Physiology:**

The SRD5A2 gene encodes the 3-oxo-5α-steroid 4-dehydrogenase 2 enzyme, also known as 5α-reductase type 2 (5αR2). This enzyme is one of three isozymes of 5α-reductase. Understanding the physiology of SRD5A2 is crucial for the context of this virtual screening tool:

- **Catalytic Role:** 5αR2 catalyzes the conversion of the male sex hormone testosterone into dihydrotestosterone (DHT). DHT is a more potent androgen with a critical role in various physiological processes, especially those related to male sexual development and maintenance.

- **Tissue Expression:** 5αR2 is primarily expressed in androgen-sensitive tissues, with high levels found in the prostate. It plays a crucial role in regulating the androgenic effects in these tissues.

- **Deficiencies:** Deficiencies in 5αR2 activity can lead to a medical condition known as 5α-reductase 2 deficiency. This deficiency results in a condition known as 46,XY Disorder of Sex Development (DSD), characterized by atypical male genitalia. Therefore, understanding and modulating 5αR2 activity is not only of pharmacological interest but also has clinical relevance.

**Usage:**

## Setup

1. First, activate the conda environment dedicated to SRD5A2 InSilicoAssay:

    ```bash
    conda activate SRD5A2_InSilicoAssay
    ```

2. Set the SCHRODINGER environment variable (assuming you have Schrödinger suite installed):

    ```bash
    export SCHRODINGER=/opt/schrodinger2020-2
    ```

3. Ensure that the LM_LICENSE_FILE environment variable is set. This variable should be configured with your Schrödinger license information.

## Input

Molecules for screening should be supplied in the SMILES format and placed in the 'input' folder. The tool expects a file called 'Cosmetics.smi' to be present in the input folder. This file should contain a list of the compounds to be screened, provided in the SMILES format.

## Execution

To start the virtual screening process, execute the following command:

```bash
nohup snakemake --cores 1 --config ligands=Cosmetics.smi&
```

- `nohup` is used to run the process in the background and ensure it continues running even if you log out.
- `snakemake` is the workflow management system used to execute the virtual screening.
- `--cores 1` specifies the number of CPU cores to use. You can adjust this based on your system's capacity.
- `--config ligands=Cosmetics.smi` configures the ligands to be screened using the 'Cosmetics.smi' file in the input folder.

**Output:**

The virtual screening process will generate results that can be found in the output directory. The results typically include information about potential SRD5A2 inhibitors and their associated properties.

**Note:**

- Make sure to configure your Schrödinger license information properly.
- Ensure that the 'Cosmetics.smi' file contains the necessary compound information in SMILES format.
- The tool can be customized further as needed to meet specific screening requirements.

For more information and troubleshooting, please refer to the documentation or contact the developers of SRD5A2 InSilicoAssay.

**References:**

Manuscript in preparation
