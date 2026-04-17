python3 balmforth_base_state_profile_txt_mod.py --kappa 1e-4 --a 0.2 --Gamma 8 --n-rows 128 --output profile.txt
# make thixo_periodic_rollwave.tst
qcc -O2 -Wall -Wdimensions -o thixo_periodic_rollwave thixo_periodic_rollwave.c -lm
./thixo_periodic_rollwave