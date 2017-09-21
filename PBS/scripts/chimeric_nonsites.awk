#!/bin/awk -f
BEGIN{
}
{
printf("%s %s %s %d %s %s %s %d", \$1, \$2, \$3, 0, \$4, \$5, \$6, 1); 
for (i=7; i<=NF; i++) {
    printf(" %s",\$i);
    };
printf("\n");
}
