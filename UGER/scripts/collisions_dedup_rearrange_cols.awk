# Rearrange columns
# Two-pass algorithm; could have made it one pass but the behavior would
# have been confusing
# First find unique chromosomes represented in file and give them an index 
# number
# Then multiple chromosome index by 1B and add it to position; sort based
# on this number
BEGIN {
  while ((getline < fname) > 0){
    for (i=3; i<=NF; i+=7){
      chrs[$i]
    }
  }
  ind=1;
  for (i in chrs) {
    a[ind++]=i;
  }
  delete chrs;
  asort(a,b,"@val_num_asc");
  for (i=1; i<=length(b); i++) {
    chrs[b[i]]=i;
  }
  delete a;
  delete b;
}
{
  str=""; 
  for (i=3; i<=NF; i+=7){
    # chr index acts as a hash; same chromosome sorted by position;
    # different chromosome sorted by chromosome index
    num=chrs[$i]*1000000000+$(i+1); 
    str=str" "num;
  }
  split(str,a);
  n=asorti(a,b,"@val_num_asc");  
  # here is where we rearrange the columns according to sort position
  for (j=1; j<=n; j++){ 
    start=((b[j]-1)*7 + 1); 
    printf("%s %s %s %s %s %s %s ",$start,$(start+1),$(start+2),$(start+3),$(start+4),$(start+5),$(start+6));
  } 
  printf("\n");
}
