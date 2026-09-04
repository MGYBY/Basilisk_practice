The code for Kao (1968)'s Frc.

```
python highprec_newtonian_kao_sweep.py \
    --sweep depth_ratio \
    --mode both \
    --minimum 0.01 \
    --maximum 100 \
    --points 81 \
    --spacing log \
    --density-ratio 0.8 \
    --viscosity-ratio 1.0 \
    --output results_depth_both.txt
```

```
python highprec_newtonian_kao_sweep.py \
    --sweep density_ratio \
    --mode surface \
    --minimum 0.10 \
    --maximum 0.99 \
    --points 90 \
    --spacing linear \
    --depth-ratio 1.0 \
    --viscosity-ratio 1.0 \
    --output results_density_surface.txt
```
