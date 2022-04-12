cat $1 | awk '{if(NR%4==1){print substr($0, 1, length($0)-21)}else{print $0}}' > $2
