# python3 balmforth_base_state_profile_txt_mod.py --kappa 1e-4 --a 0.2 --Gamma 8 --n-rows 128 --output profile.txt
# # make thixo_periodic_rollwave.tst
# qcc -O2 -Wall -Wdimensions -o thixo_periodic_rollwave thixo_periodic_rollwave_consistent.c -lm
# ./thixo_periodic_rollwave

python3 thixo_base_state_explicitFrV_fixed.py --kappa 1e-4 --a 0.2 --Gamma 8 --FrV 4.00 --n-rows 128 --output profile.txt
qcc -O2 -Wall -Wdimensions -o thixo_periodic_rollwave thixo_periodic_rollwave_explicitFrV_fixed.c -lm
./thixo_periodic_rollwave