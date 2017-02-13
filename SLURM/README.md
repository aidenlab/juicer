#Quick Start

We've provided a step-by-step guide to showcase some of the features of
Juicer. If you run into problems, see below for more detailed documentation.
This example runs on SLURM clusters, but you can also install the pipeline
on any LSF or Univa Grid Engine cluster.

1. Log in to the SLURM cluster
2. Create a custom directory (e.g. "mkdir -p /custom/filepath/MyHIC")
3. Go to the directory (e.g. "cd /custom/filepath/MyHIC")
4. Create a folder for fastq files ("mkdir fastq")
5. Use globus, other ftp software, soft linking, etc to move your fastq files to /custom/filepath/MyHIC/fastq
6. Run the Juicer pipeline from the MyHIC folder with "/local/path/scripts/juicer.sh [options]"
   where /local/path refers to the folder containing the scripts folder bundling the necessary files included
   with this distribution.
7. Do not exit the screen or kill the script until you see a message saying that "all jobs" have been submitted.
8. You can use "squeue" to check the status of jobs
9. Eventually there will be no more jobs in the queue, and the ./debug folder will have a "Pipeline successfully completed" message.
10. Results are available in the aligned directory. The Hi-C maps are ininter.hic (for MAPQ > 0) and inter_30.hic (for MAPQ >= 30). The Hi-C maps can be loaded in Juicebox and explored. They can also be used for automatic feature annotation and to extract matrices at specific resolutions.

   These results also include automatic feature annotation. The output files include a genome-wide annotation of loops and, whenever possible, the CTCF motifs that anchor them (identified using the HiCCUPS algorithm). The files also include a genome-wide annotation of contact domains (identified using the Arrowhead algorithm). The formats of these files are described in the Juicebox tutorial online; both files can be loaded into Juicebox as a 2D annotation.
11.To download a file (e.g. inter.hic) from AWS to load into Juicebox, type:
```
      sftp -i Juicer_AWS.pem ubuntu@<given IP address>
      cd /opt/juicer/work/MBR19/aligned
      get inter.hic
      get inter_30.hic
      get ... (each of hiccups, apa, motifs, and arrowhead output files)
```
