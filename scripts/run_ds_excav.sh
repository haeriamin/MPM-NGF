#!/bin/bash

scale_dim=0.25
scale_vel=$(echo "scale=10; 2*$scale_dim" | bc -l) 

c=1
# Depths [m]
for d in $(echo "scale=10; 0.02*$scale_dim" | bc -l) \
         $(echo "scale=10; 0.04*$scale_dim" | bc -l) \
         $(echo "scale=10; 0.05*$scale_dim" | bc -l) \
         $(echo "scale=10; 0.08*$scale_dim" | bc -l) \
         $(echo "scale=10; 0.10*$scale_dim" | bc -l)
do
    # Horizontal speeds [m/s]
    for s in $(echo "scale=10; 0.01*$scale_vel" | bc -l) \
             $(echo "scale=10; 0.04*$scale_vel" | bc -l) \
             $(echo "scale=10; 0.08*$scale_vel" | bc -l) \
             $(echo "scale=10; 0.10*$scale_vel" | bc -l) \
             $(echo "scale=10; 0.15*$scale_vel" | bc -l)
    do
        # Angles [deg]
        for a in 0 3.8 10 30.8 45
        do
            # Motion types
            for m in 1 2
            do
                # echo $c $d $s $a $m
                python3 ds_excav.py $c $d $s $a $m
                c=`expr $c + 1`
            done
        done
    done
done

# Instructions:
# Mark it executable using:
# $ chmod +x run_ds_excav.sh
# Then run it:
# $ ./run_ds_excav.sh