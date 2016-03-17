$1=="Sequenced"{
gsub(/,/,"",$4);
 total=total+$4;
}
$1=="Unmapped:"{
gsub(/,/,"",$2);
 unmapped=unmapped+$2
}
$1=="Normal"{
gsub(/,/,"",$3);
 regular=regular+$3}
$1=="Chimeric" && $2=="Paired:"{
gsub(/,/,"",$3);
 normal=normal+$3}
$2=="Ambiguous:"{
gsub(/,/,"",$3);
 abnorm=abnorm+$3}
$1=="Alignable"{
gsub(/,/,"",$4);
 alignable=alignable+$4}
$1=="Unique" {
gsub(/,/,"",$3);
 dedup=dedup+$3}
$1=="PCR"{
gsub(/,/,"",$3);
 dups=dups+$3}
$1=="Optical"{
gsub(/,/,"",$3);
 optdups=optdups+$3}
END{
 printf("%s %'d\n", "Sequenced Read Pairs:", total);
 printf(" %s %'d (%0.2f%)\n", "Normal Paired:", regular, regular*100/total);
 printf(" %s %'d (%0.2f%)\n", "Chimeric Paired:", normal, normal*100/total);
 printf(" %s %'d (%0.2f%)\n", "Chimeric Ambiguous:", abnorm, abnorm*100/total);
 printf(" %s %'d (%0.2f%)\n", "Unmapped:", unmapped, unmapped*100/total);
 printf(" %s %'d (%0.2f%)\n", "Alignable (Normal+Chimeric Paired):", alignable, alignable*100/total);
 printf("%s %'d\n", "Unique Reads:", dedup);
 printf("%s %'d\n", "PCR Duplicates:", dups);
  printf("%s %'d\n", "Optical Duplicates:", optdups);
}
