#!/bin/sh
# C++ rivals on the held-out CEC 2014 suite, 25 runs, D = 30 and 50.
set -e
for D in 30 50; do
  # jSO (official CEC 2017 code, linked to CEC 2014 functions), one process per function, fixed seeds
  mkdir -p /home/claude/heldout14/jso_d$D
  for f in $(seq 1 30); do
    o=/home/claude/heldout14/jso_d$D/F$f.txt
    [ -s $o ] && grep -c "error value" $o | grep -q 25 && continue
    /home/claude/rivals14/jso/jso14 $D $f 25 $((14000+f)) > $o
  done
done
