# TractoEval pipeline
===================================

Runs the TractoEval pipeline. This Nextflow process runs RecobundlesX in a manner very similar to what has been done in Tractoinferno. This allows you to compare your own tractogram with Tractoinferno's, for any subject. TractoEval will segment your tractogram into bundles. You will obtain scores such as the Dice score for each bundle, the overlap and overreach ratio, and so on.

To access the RBX-flow example atlases: https://zenodo.org/record/4630660#.YJvmwXVKhdU.

If you use this pipeline, please cite:

```
Poulin, P., Theaud, G., Rheault, F., St-Onge, E., Bore, A., Renauld, E., ... 
& Descoteaux, M. (2022). TractoInferno-A large-scale, open-source, multi-site 
database for machine learning dMRI tractography. Scientific Data, 9(1), 725.
```

Tractoeval uses many published processes:

```
Rheault, Francois. Analyse et reconstruction de faisceaux de la matière blanche.
page 137-170, (2020), https://savoirs.usherbrooke.ca/handle/11143/17255

Garyfallidis, Eleftherios, et al. Recognition of white matter bundles using
local and global streamline-based registration and clustering.
NeuroImage 170 (2018) https://doi.org/10.1016/j.neuroimage.2017.07.015

Kurtzer GM, Sochat V, Bauer MW Singularity: Scientific containers for
mobility of compute. PLoS ONE 12(5): e0177459 (2017)
https://doi.org/10.1371/journal.pone.0177459

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows.
Nature Biotechnology 35, 316–319 (2017) https://doi.org/10.1038/nbt.3820
```

> [!CAUTION]
> Beware! Even when using the option `--fully_reproducible`, this code is a very difficult code to use for reproducible results! 
>
> RBX depends strongly on many values, such as the random seed (rng) or the order of streamlines in the tractogram. For instance, running again this code on wholebrain tractograms of subjects in the tractoinferno database with shuffled streamlines can lead to scores very far from 100% dice scores, and as low as 60% for some bundles.


> [!CAUTION]
> Be careful with your interpretation! The final bundles are *not* the group of streamlines included in your tractogram that look similar to the atlas's bundles. The segmentation depends on many parameters. In the case of this Population Atlas, parameters are quite strict and do not result in very large bundles.
> 
> At the time of creating Tractoinferno, RBX was state-of-the-art. More recent versions of RBX include BundleSeg (2023: see St-Onge et. al.), which is still not fully reproducible but contains less parameters and runs faster.


Requirements
------------

The main requirement is nextflow, version 21.*.

- [Nextflow](https://www.nextflow.io)

For instance, you may find version 21.10.0 from the Github releases: https://github.com/nextflow-io/nextflow/releases/tag/v21.10.0. Nextflow also requires Java. For version 21.10.0, version up to 17 were supported. To verify your Java version: `java -version`.

```
curl -sL https://github.com/nextflow-io/nextflow/releases/download/v21.10.0/nextflow-21.10.0-all -o nextflow
chmod +x nextflow
```

Then, you can either run nextflow from this local artifact (e.g. `./nextflow -version`), or add this path in your `$PATH` in your .bashrc, or move it to a place already on `$PATH` such as `/usr/local/bin/` or `$VIRTUAL_ENV/bin/`. 

Nextflow will run code that also requires other libraries listed below, but to make sure that you use the correct versions (i.e. the same one as was used when creating tractoinferno), we recommand running the code within a singularity (see below for explanation). 

- [scilpy](https://github.com/scilus/scilpy)
- [ANTs](https://github.com/ANTsX/ANTs)


Singularity/Docker
-----------

We recommend using the command in an environment already prepared by us, where all requirements (ex, scilpy, ants) are already installed. This can be done either with Singularity or Docker. In both cases, the environment name is scilus:1.2.0.


**Singularity**

If you are on Linux, we recommand using Singularity (which is now renamed Apptainer).

- [Singularity](https://apptainer.org/docs/admin/main/installation.html)

Then you may use singularity to build the environment into an image (.sif):

`apptainer build scilus_1.2.0.sif docker://scilus/scilus:1.2.0`

Then you can run our code through Nextlow with your favorite options (see the file [usage](https://github.com/scil-vital/TractoInferno/blob/main/USAGE)) and add this in the command line:

`-with-singularity scilus_1.2.0.sif`

**Docker**

If you are on MacOS or Windows, we recommend using the Docker container to run rbx_flow pipeline.

- [Docker](https://www.docker.com/)

Then, launch your Nextflow command with your favorite options (see the file [usage](https://github.com/scil-vital/TractoInferno/blob/main/USAGE)) and add this in the command line. It will find and pull the environment automatically.

`-with-docker scilus/scilus:1.2.0`

Usage
-----

See [usage](https://github.com/scil-vital/TractoInferno/blob/main/USAGE) or run `nextflow run main.nf --help`.

Hints:

- It is common usage to run the command line directly from the output folder. The work folder created by Nextflow will be created in the directory from where you lauch the command line.

- The final output will contain symlinks to files inside the work folder. To create a folder with direct files, you may use: `rsync -rvL $output_dir_links $output_dir_files`

Output
------

The output will contain folders Aggregate_Measures and one folder per subject.  The Aggregate_Measures present the statistics for all subjects.

Reminder: there are a total of 39 bundles in the atlas. Bundles available as --reference for each subject are the recovered bundles that passed visual QC.  Missing bundles were either not reconstructed in Tractoinferno processing, or not of enough quality to be kept. If you did recover them in your own tractogram, they will *not* participate in the scoring.

To view all bundles recovered by RBX for your tractogram, even those that were not included in the scoring, view outputs SubjX/Clean_Bundles/*.trk.