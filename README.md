# Read this first!!

To access Juicer 1.6 (last stable release), please see [the Github Release](https://github.com/aidenlab/juicer/releases/tag/1.6). If you clone the Juicer repo directly from Github, it will clone Juicer 2, which is under active development. If you encounter any bugs, please let us know.

ENCODE's Hi-C uniform processing pipeline based on Juicer can be found [here](https://github.com/ENCODE-DCC/hic-pipeline).

# About Juicer

Juicer is a platform for analyzing kilobase resolution Hi-C data. In this
distribution, we include the pipeline for generating Hi-C maps from fastq raw
data files and command line tools for feature annotation on the Hi-C maps.

The beta release for Juicer version 1.6 can be accessed via [the Github Release](https://github.com/aidenlab/juicer/releases/tag/1.6). The main repository on Github is now focused on the Juicer 2.0 release and is under active development.
For general questions, please use
[the Google Group](https://groups.google.com/forum/#!forum/3d-genomics).

If you are interested in running Juicer in the cloud, you may want to check out the dockerized version 
of Juicer hosted by [ENCODE](https://github.com/ENCODE-DCC/hic-pipeline).

If you have any difficulties using Juicer, please do not
hesitate to contact us (aidenlab@bcm.edu)

**If you use Juicer in your research, please cite:
Neva C. Durand, Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. "Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments." Cell Systems 3(1), 2016.**

# Documentation
Please see [the wiki](https://github.com/aidenlab/juicer/wiki) for extensive documentation.

# Questions?
For FAQs, or for asking new questions, please see our forum: [aidenlab.org/forum.html](http://aidenlab.org/forum.html).

------------
Distribution
------------

In this repository, we include the scripts for running Juicer on AWS, LSF,
Univa Grid Engine, SLURM, and a single CPU.

The SLURM and CPU scripts are the most up to date. For cloud computing, we recommend the [ENCODE uniform processing pipeline based on Juicer](https://github.com/ENCODE-DCC/hic-pipeline)

/SLURM - scripts for running pipeline and postprocessing on SLURM

/CPU - scripts for running pipeline and postprocessing on a single CPU

/AWS - scripts for running pipeline and postprocessing on AWS **Deprecated**

/UGER - scripts for running pipeline and postprocessing on UGER (Univa) **Deprecated**

/LSF - scripts for running pipeline and postprocessing on LSF **Deprecated**

/misc - miscellaneous helpful scripts

----------------------------------
Hardware and Software Requirements
----------------------------------
Juicer is a pipeline optimized for parallel computation on a cluster. Juicer
consists of two parts: the pipeline that creates Hi-C files from raw data,
and the post-processing command line tools.

### Cluster requirements:

Juicer requires the use of a cluster or the cloud, with ideally >= 4 cores (min 1 core)
and >= 64 GB RAM (min 16 GB RAM)

Juicer currently works with the following resource management software:
- [OpenLava](http://www.openlava.org/)
- [LSF](https://www.ibm.com/systems/spectrum-computing/products/lsf)
- [SLURM](https://slurm.schedmd.com/download.html)
- GridEngine (Univa, etc. any flavor)

We recommend [ENCODE's Hi-C processing pipeline, based on Juicer](https://github.com/ENCODE-DCC/hic-pipeline) to run in the cloud; the AWS scripts are out of date.

### Juicer tools requirements

The minimum software requirement to run Juicer is a working Java installation
(version >= 1.8) on Windows, Linux, and Mac OSX.  We recommend using the
latest Java version available, but please do not use the Java Beta Version.
Minimum system requirements for running Java can be found at
https://java.com/en/download/help/sysreq.xml

To download and install the latest Java Runtime Environment (JRE), please go
to https://www.java.com/download


### GNU CoreUtils

The latest version of GNU coreutils can be downloaded from
https://www.gnu.org/software/coreutils/manual/

### Burrows-Wheeler Aligner (BWA)

The latest version of BWA should be installed from
http://bio-bwa.sourceforge.net/

### CUDA (for HiCCUPS peak calling)

You must have an NVIDIA GPU to install CUDA.

Instructions for installing the latest version of CUDA can be found on the
[NVIDIA Developer site](https://developer.nvidia.com/cuda-downloads).

The native libraries included with Juicer are compiled for CUDA 7 or CUDA 7.5.
See the [download page for Juicer
Tools](https://github.com/theaidenlab/juicer/wiki/Download).

Other versions of CUDA can be used, but you will need to download the
respective native libraries from
[JCuda](http://www.jcuda.org/downloads/downloads.html).

For best performance, use a dedicated GPU. You may also be able to obtain
access to GPU clusters through Amazon Web Services, Google cloud, or a local research
institution.

If you cannot access a GPU, you can run the [CPU version of HiCCUPS](https://github.com/aidenlab/juicer/wiki/CPU-HiCCUPS) directly using the `.hic` file and Juicer Tools.

### Building new jars

See the Juicebox documentation at <https://github.com/theaidenlab/Juicebox> for
details on building new jars of the juicer_tools.

------------
Quick Start
------------
Run the Juicer pipeline on your cluster of choice with "juicer.sh [options]"

```
Usage: juicer.sh [-g genomeID] [-d topDir] [-q queue] [-l long queue] [-s site]
                 [-a about] [-R end] [-S stage] [-p chrom.sizes path]
                 [-y restriction site file] [-z reference genome file]
                 [-C chunk size] [-D Juicer scripts directory]
                 [-Q queue time limit] [-L long queue time limit] [-e] [-h] [-x]
* [genomeID] must be defined in the script, e.g. "hg19" or "mm10" (default
  "hg19"); alternatively, it can be defined using the -z command
* [topDir] is the top level directory (default
  "/Users/nchernia/Downloads/neva-muck/UGER")
     [topDir]/fastq must contain the fastq files
     [topDir]/splits will be created to contain the temporary split files
     [topDir]/aligned will be created for the final alignment
* [queue] is the queue for running alignments (default "short")
* [long queue] is the queue for running longer jobs such as the hic file
  creation (default "long")
* [site] must be defined in the script, e.g.  "HindIII" or "MboI"
  (default "none")
* [about]: enter description of experiment, enclosed in single quotes
* [stage]: must be one of "chimeric", "merge", "dedup", "final", "postproc", or "early".
    -Use "chimeric" when alignments are done but chimeric handling has not finished
    -Use "merge" when alignment has finished but the merged_sort file has not
     yet been created.
    -Use "dedup" when the files have been merged into merged_sort but
     merged_nodups has not yet been created.
    -Use "final" when the reads have been deduped into merged_nodups but the
     final stats and hic files have not yet been created.
    -Use "postproc" when the hic files have been created and only
     postprocessing feature annotation remains to be completed.
    -Use "early" for an early exit, before the final creation of the stats and
     hic files
* [chrom.sizes path]: enter path for chrom.sizes file
* [restriction site file]: enter path for restriction site file (locations of
  restriction sites in genome; can be generated with the script
  (misc/generate_site_positions.py) )
* [reference genome file]: enter path for reference sequence file, BWA index
  files must be in same directory
* [chunk size]: number of lines in split files, must be multiple of 4
  (default 90000000, which equals 22.5 million reads)
* [Juicer scripts directory]: set the Juicer directory,
  which should have scripts/ references/ and restriction_sites/ underneath it
  (default /broad/aidenlab)
* [queue time limit]: time limit for queue, i.e. -W 12:00 is 12 hours
  (default 1200)
* [long queue time limit]: time limit for long queue, i.e. -W 168:00 is one week
  (default 3600)
* -f: include fragment-delimited maps from hic file creation
* -e: early exit
* -h: print this help and exit
```

------------
Juicer Usage
------------
- **Running Juicer with no arguments** will run it with genomeID hg19 and site MboI
- **Providing a genome ID**: if not defined in the script, you can either directly modify the script or provide the script with the files needed. You would provide the script with the files needed via "-z reference_sequence_path" (needs to have the BWA index files in same directory), "-p chrom_sizes_path" (these are the chromosomes you want included in .hic file), and "-s site_file" (this is the listing of all the restriction site locations, one line per chromosome). Note that ligation junction won't be defined in this case.  The script (misc/generate_site_positions.py) can help you generate the file
- **Providing a restriction enzyme**: if not defined in the script, you can either directly modify the script or provide the files needed via the "-s site_file" flag, as above.  Alternatively, if you don't want to do any fragment-level analysis (as with a DNAse experiment), you should assign the site "none", as in `juicer.sh -s none`
- **Directory structure**: Juicer expects the fastq files to be stored in a directory underneath the top-level directory. E.g. HIC001/fastq.  By default, the top-level directory is the directory where you are when you launch Juicer; you can change this via the -d flag. Fastqs can be zipped. [topDir]/splits will be created to contain the temporary split files and should be deleted once your run is completed.  [topDir]/aligned will be created for the final files, including the hic files, the statistics, the valid pairs (merged_nodups), the collisions, and the feature annotations.
- **Queues** are complicated and it's likely that you'll have to modify the script for your system, though we did our best to avoid this.  By default there's a short queue and a long queue.  We also allow you to pass in wait times for those queues; this is currently ignored by the UGER and SLURM versions.  The short queue should be able to complete alignment of one split file.  The long queue is for jobs that we expect to take a while, like writing out the merged_sort file
- **Chunk size** is intimitely associated with your queues; a smaller chunk size means more alignment jobs that complete in a faster time.  If you have a hard limit on the number of jobs, you don't want too small of a chunk size.  If your short queue has a very limited runtime ceiling, you don't want too big of a chunk size.  Run time for alignment will also depend on the particulars of your cluster.  We launch ~5 jobs per chunk.  Chunk size must be a multiple of 4.
-  **Relaunch** via the same script. Type `juicer.sh [options] -S stage` where "stage" is one of merge, dedup, final, postproc, or early. "merge" is for when alignment has finished but merged_sort hasn't been created; "dedup" is for when merged_sort is there but not merged_nodups (this will relaunch all dedup jobs); "final" is for when merged_nodups is there and you want the stats and hic files; "postproc" is for when you have the hic files and just want feature annotations; and "early" is for early exit, before hic file creation. If your jobs failed at the alignment stage, run `relaunch_prep.sh` and then run juicer.sh.
- **Miscelleaneous options** include -a 'experiment description', which will add the experiment description to the statistics file and the meta data in the hic file; -r, which allows you to use bwa aln instead of bwa mem, useful for shorter reads; -R [end], in case you have one read end that's short and one that's long and you want to align the short end with bwa aln and the long end with bwa mem; and -D [Juicer scripts directory], to set an alternative Juicer directory; must have scripts/, references/, and restriction_sites/ underneath it

------------------------
Command Line Tools Usage
------------------------
Detailed documentation about the command line tools can be found on the  wiki:

* [Annotating features with Arrowhead, HiCCUPS, MotifFinder, APA, Eigenvector, and Pearsons](https://github.com/aidenlab/juicer/wiki/Feature-Annotation)
* [Creating .hic with Pre](https://github.com/aidenlab/juicer/wiki/Pre)
* [Extracting data from .hic files with straw](https://github.com/aidenlab/straw)

To launch the command line tools, use the shell script “juicer_tools” on Unix/MacOS
or type
```
java -jar juicer_tools.jar (command...) [flags...] <parameters...>`
```

In the command line tools, there are several analysis functions:

1.  `apa` for conducting aggregate peak analysis
1.  `hiccups` for annotating loops
1.  `motifs` for finding CTCF motifs
1.  `arrowhead` for annotating contact domains
1.  `eigenvector` for calculating the eigenvector (first PC) of the Pearson's
1.  `pearsons` for calculating the Pearson's

The `juicer_tools` (Unix/MacOS) script can be used in place of the unwieldy
		`java -Djava.library.path=path/to/natives/ -jar juicer_tools.jar`
