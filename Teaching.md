______________________________________________________________________
## General course idea
______________________________________________________________________

### Bioinformatics fundamentals for pipeline development

We follow a live example of the pipeline development with example of a case study. The take a raw data from the paper and develop the pipeline with these tools.

* Version control
    * Git
    * Environment Management
        * Conda
	* Workflow Management
	    * Snakemake and Nextflow
	    * Report and Documentation
	        * R markdown and Jupyter
		* Containerazation
		    * Docker and Singularity
		    ______________________________________________________________________
		    Git
		    * Tool for version control and collaborating on code
		    * Track changes and edits
		    * Commit changes and
		    * pushing the commit to the remote repository
		    * Tracking all the edits and handle potential conflicts
		    * Using cloud based hosting such as Github or bitbucket
		    * Distribute the code

Conda
* Package and environment management tool
* conda install ...
* Recreate the system that is used to generate results
* conda create environment
* hosted and downloaded from channel
* widely used channels are conda-forge and Bioconda
* Conda environment allows to add PATH

Snakemake
* Workflow management systems (WMS) performs and monitors a defined sequence of computational tasks.
* Snakemake is a bioinformatics community based on python
* But python is not needed
* Can be scaled up to desktop/Laptop to server, cluster, grid or cloud
* Small subset of data and then real one in cluster
* Works on files (rather than for streams, reading/ writing from database)
* Specially with next-generation sequencing data which is computational expensive operations on large files
* WMS are tools for reproducbible 