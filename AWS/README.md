#Quick Start

Extensive instructions, including how to create and launch an instance, are 
available at <http://aidenlab.org/juicer/aws.html>

Below is the set of instructions we gave to reviewers. You can obtain the 
MBR19 test data set at  
<http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/MBR19/fastq/chr19_R1.fastq.gz>  
<http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/MBR19/fastq/chr19_R1.fastq.gz>  
The HIC003 test data set is available at
<http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R1_001.fastq.gz>  
<http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R2_001.fastq.gz>  

You will need to launch your own instance, use your own Juicer_AWS.pem, and
use the public IP associated with your instance.

------------------------------------------------------------------------------
We've provided a step-by-step guide to showcase some of the features of
Juicer. If you run into problems, see below for more detailed documentation.
This example runs on Amazon Web Services, but you can install the pipeline
on any LSF, Univa Grid Engine, or SLURM cluster.

1. Make sure you're in the top-level directory, with the Juicer_AWS.pem file.
2. You were given an anonymous IP address. At a command line prompt, type:
     ` ssh -i Juicer_AWS.pem ubuntu@<given IP address>`
3. This will log you into an AWS instance that contains all the software
   needed to run the pipeline. Type
      `cd /opt/juicer/work/`
4. We will run the pipeline on a test dataset of a single chromosome of the primary+
   replicate map from (Rao+Huntley et al., 2014). Type:
      `cd MBR19`
5. Run the Juicer pipeline on the raw data, which is stored in the fastq
   directory:
      `/opt/juicer/scripts/juicer.sh -g hg19 -s MboI`
6. You will see a series of messages sending jobs to the cluster. Do not
   kill the script or close the server connection until you see:
      `“(-: Finished adding all jobs... please wait while processing.”`
7. At this point you can close the connection and come back later. 
   To see the progress of the pipeline as it works, type:
      `bjobs -w`
8. Eventually the bjobs command will report “No unfinished job found”. Type:
      `tail lsf.out`
   You should see “(-: Pipeline successfully completed (-:”
9. Results are available in the aligned directory. The Hi-C maps are in
   inter.hic (for MAPQ > 0) and inter_30.hic (for MAPQ >= 30). The Hi-C maps
   can be loaded in Juicebox and explored. They can also be used for
   automatic feature annotation and to extract matrices at specific
   resolutions.
   These results also include automatic feature annotation. The output files include 
   a genome-wide annotation of loops and, whenever possible, the CTCF motifs that anchor 
   them (identified using the HiCCUPS algorithm). The files also include a genome-wide 
   annotation of contact domains (identified using the Arrowhead algorithm). The formats 
   of these files are described in the Juicebox tutorial online; both files can be loaded 
   into Juicebox as a 2D annotation.
10. To download a file (e.g. inter.hic) from AWS to load into Juicebox, type:

   ```
      sftp -i Juicer_AWS.pem ubuntu@<given IP address>  
      cd /opt/juicer/work/MBR19/aligned  
      get inter.hic  
      get inter_30.hic  
      get ... (each of hiccups, apa, motifs, and arrowhead output files)  
   ```

11. You can also run the pipeline on genome-wide dataset that is lower resolution. Type
      `cd /opt/juicer/work/HIC003`
   Then
      `/opt/juicer/scripts/juicer.sh -g hg19 -s MboI`
   Again the pipeline will run. The results will be available in the aligned directory.
   Because this is not a deeply sequenced map, loop lists and domain lists will not be 
   produced.
