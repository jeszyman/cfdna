# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Repository%20setup%20and%20administration][Repository setup and administration:2]]
#!/usr/bin/env bash

#################################################################
###   Tangle An Org-Mode File Through Non-Interactive Emacs   ###
#################################################################

org_file="${1}"


/usr/bin/emacs --batch -l org --eval "(progn (find-file \"$org_file\") (org-babel-tangle))"
# Repository setup and administration:2 ends here
