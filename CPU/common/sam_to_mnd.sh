# How to create old style merged_nodups.txt for 3D-DNA
# File should be dedupped SAM
saamtools view -O SAM -F 1024 $1 | awk -v mnd=1 -f sam_to_pre.awk 
