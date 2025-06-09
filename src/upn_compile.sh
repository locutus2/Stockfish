#!/bin/bash
make clean
touch search.cpp && make build fail_low=no research=no n_conditions=100 && mv stockfish stockfish_upn_fh_100
touch search.cpp && make build fail_low=yes research=no n_conditions=100 && mv stockfish stockfish_upn_fl_100
touch search.cpp && make build fail_low=no research=no n_conditions=2500 && mv stockfish stockfish_upn_fh_2500
touch search.cpp && make build fail_low=yes research=no n_conditions=2500 && mv stockfish stockfish_upn_fl_2500
touch search.cpp && make build fail_low=no research=no n_conditions=1000 && mv stockfish stockfish_upn_fh_1000
touch search.cpp && make build fail_low=yes research=no n_conditions=1000 && mv stockfish stockfish_upn_fl_1000
touch search.cpp && make build fail_low=no research=no n_conditions=200 && mv stockfish stockfish_upn_fh_200
touch search.cpp && make build fail_low=yes research=no n_conditions=200 && mv stockfish stockfish_upn_fl_200
touch search.cpp && make build fail_low=no research=no n_conditions=400 && mv stockfish stockfish_upn_fh_400
touch search.cpp && make build fail_low=yes research=no n_conditions=400 && mv stockfish stockfish_upn_fl_400
