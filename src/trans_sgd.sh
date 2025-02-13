#/bin/bash
grep Result |cut -d' ' -f2-|tr '. ' ',;'|sed 's/-nan//'
