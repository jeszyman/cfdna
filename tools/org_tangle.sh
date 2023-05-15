#!/usr/bin/env bash

#################################################################
###   Tangle An Org-Mode File Through Non-Interactive Emacs   ###
#################################################################

org_file="${1}"


/usr/bin/emacs --batch -l org --eval "(progn (find-file \"$org_file\") (org-babel-tangle))"
