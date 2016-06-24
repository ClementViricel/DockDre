#!/bin/bash
source /home/dsimoncini/softs/ccp4-6.4.0/setup-scripts/ccp4.setup-sh

pdbset xyzin $1 xyzout $2 1>/dev/null <<EOF
renumber 1 
exclude HOH 
exclude HEADER 
EOF

pdbcur xyzin $2 xyzout $2 1>/dev/null <<EOF
renchain /*/A 'A'
renchain /*/B 'A'
renchain /*/C 'B'
renchain /*/D 'B'
EOF
