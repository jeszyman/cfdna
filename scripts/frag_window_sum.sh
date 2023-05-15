#!/usr/bin/env bash
input_frag="$1"
output_short="$2"
output_long="$3"

# Functions
make_short(){
    cat $1 | awk '{if ($4 >= 100 && $5 <= 150) print $0}' > $2
}

make_long(){
    cat $1 | awk '{if ($4 >= 151 && $5 <= 220) print $0}' > $2
}

# Run command
make_short $input_frag $output_short
make_long $input_frag $output_long
