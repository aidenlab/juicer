#Quick Start
Run the Juicer pipeline on Univa with "[Univa directory]/scripts/juicer.sh [options]"
```
 Usage: juicer.sh -g genomeID [-d topDir] [-q queue] [-l long queue] [-s site] [-k key] [-a about] 
 [-R end] [-S stage] [-p chrom.sizes path] [-z reference sequence path] [-h]
   genomeID must be defined in the script, e.g. "hg19" or "mm10" (default "hg19")   
   alternatively, it can be defined using the -z command
   [topDir] is the top level directory (default current directory)
     [topDir]/fastq must contain the fastq files
     [topDir]/splits will be created to contain the temporary split files
     [topDir]/aligned will be created for the final alignment
   [queue] is the queue for running alignments (default "short")
   [long queue] is the queue for running longer jobs such as the hic file creation (default "long")
   [site] must be defined in the script, e.g.  "HindIII" or "MboI" (default "DpnII"); alternatively, this can be the restriction site file
   [end]: use the short read aligner on read end, must be one of 1 or 2 
   -r: use the short read version of the aligner, bwa aln (default: long read, bwa mem)
   [stage]: must be one of "merge", "dedup", "final", or "early". 
      Use "merge" when alignment has finished but the merged_sort file has not yet been created. 
      Use "dedup" when the files have been merged into merged_sort but merged_nodups has not yet been created. 
      Use "final" when the reads have been deduped into merged_nodups but the final stats and hic files have not yet been created
      Use "early" for an early exit, before the final creation of the stats and hic files
   -a: enter description of experiment, enclosed in single quotes
   -p: enter path for chrom.sizes file
   -z: enter path for reference sequence file, BWA index file must be in same directory
   -h: print this help and exit
```
#Details
Juicer aligns the reads, handles chimeras, sorts, dedups, then creates .hic files that can be loaded into Juicebox.  The pipeline now runs slightly differently in that the launching waits for the splits to be done.  This is so it's consistent with the SLURM pipeline (and it's a little more elegant from a coding standpoint).  But because of this, the terminal can hang, so I suggest you run it in a screen.  

#New capabilities

- If you have a new genome, you can either directly modify the script (as usual) or provide the script with the files needed.  You would provide the script with the files needed via "-z reference_sequence_path" (needs to have the BWA index files in same directory), "-p chrom_sizes_path" (these are the chromosomes you want included in .hic file), and "-s site_file" (this is the listing of all the restriction site locations, one line per chromosome).  Note that ligation junction won't be defined in this case.
- Relaunch via the same script.  Type "juicer.sh [options] -S stage" where "stage" is one of merge, dedup, final, early.  "merge" is for when alignment has finished but merged_sort hasn't been created; "dedup" is for when merged_sort is there but not merged_nodups (this will relaunch all dedup jobs); "final" is for when merged_nodups is there and you want the stats and .hic files; and "early" is for early exit.  If your jobs failed at the alignment stage, run "relaunch_prep" and then run juicer.sh.    
- For small files (QC runs, or technically when your data is less than 2GB), you can keep the files gzipped.  Just copy them to the fastq folder as usual but no need to unzip. If you do unzip, it will still work as usual.

#Steps for launching

1. Log into server 
2. Navigate to your directory with your fastqs
Your fastqs should be in a directory titled "fastq" underneath top level directory.
They can be in gzipped form if it's a small QC run; deep sequencing should be unzipped

3. Type screen
 screen
4. Run Juicer with options juicer.sh [options]
The options are things like "-g hg19 -s MboI" etc.

To see the options type "juicer.sh -h"

5. To get back to your terminal, type 
 Ctrl-A Ctrl-D
If you forget and then close your computer, no problem, the jobs will still launch.

6. To check back later that the jobs have been launched, log back into the original server and type 
 screen -r
7. To check on your jobs, type 
 qstat
8. If there are no jobs in your queue, most likely Juicer has successfully completed.  Type
 tail uger.out
You should see
 (-: Pipeline successfully completed (-:
 Run cleanup.sh to remove the splits directory
9. Go ahead and do as the script says:
 [Univa directory]/scripts/cleanup.sh

#Sample output
The "exit code 0" statement means that the split successfully completed.

 (-: Looking for fastq files...fastq files exist
 
 Prepending: UGER (already loaded)
 
 (-: Aligning files matching HIC001/fastq/*_R*.fastq
  in queue short to genome hg19
 
 (-: Created HIC001/splits and HIC001/aligned.  Splitting files
 
 Your job 95416 ("a1439405283split0") has been submitted
 
 Job 95416 exited with exit code 0. 
 
 Your job 95419 ("a1439405283split1") has been submitted
 
 Job 95419 exited with exit code 0.
 
 (-: Starting job to launch other jobs once splitting is complete
 
 Your job 95421 ("a1439405283_001000.fastqcountligations") has been submitted
 
 Your job 95422 ("a1439405283_align1_001000.fastq") has been submitted
 
 Your job 95423 ("a1439405283_align2_001000.fastq") has been submitted
 
 Your job 95424 ("a1439405283_merge_001000.fastq") has been submitted
 
 Your job 95425 ("a1439405283_fragmerge") has been submitted
 
 Your job 95426 ("a1439405283_osplit") has been submitted
 
 Your job 95427 ("a1439405283_finallaunch") has been submitted
 
 Your job 95428 ("a1439405283_done") has been submitted
 
 (-: Finished adding all jobs... please wait while processing.

If you don't see "Pipeline successfully completed" and a reminder to "Run cleanup.sh" at the end of the uger.out file, something probably went wrong. 

#Using a new genome
To use a new genome with Juicer, you will need to index it using BWA, create the restriction site file, and create the chromosome sizes file that will dictate what can be viewed in Juicebox.  Then Juicer can be run using the flags above.

