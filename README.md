# TractoEval pipeline
===================================

Run the TractoEval pipeline.

To access the RBX-flow example atlases:
https://zenodo.org/record/4630660#.YJvmwXVKhdU

If you use this pipeline, please cite:

```
Poulin, Philippe, TractoInferno, ???

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

Requirements
------------

- [Nextflow](https://www.nextflow.io)
- [scilpy](https://github.com/scilus/scilpy)
- [ants](https://github.com/ANTsX/ANTs)

Singularity/Docker
-----------
If you are on Linux, we recommend using the Singularity to run the tractoeval_flow pipeline.
If you have Singularity == 3.*, launch your Nextflow command with:
`-with-singularity scilus/scilus:1.2.0_tractoinferno-1.0.0`

If you are on MacOS or Windows, we recommend using the Docker container to run rbx_flow pipeline.
Launch your Nextflow command with:
`-with-docker scilus/scilus:1.2.0_tractoinferno-1.0.0`

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`

