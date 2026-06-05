#!/usr/bin/bash
sed 's/.*\(Iteration.*\)$/\1/' | sed 's/[a-z]*=//ig' | cut -d' ' -f2-6,8- | tr '.' ','
