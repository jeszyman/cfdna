#!/usr/bin/env bash
trimmomatic PE \
            -threads $1 \
            -trimlog $2 \
            $3 $4 \
            $5 $6 \
            $7 $8 \
            ILLUMINACLIP:"$9":2:30:10 \
            LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:20
