#!/bin/bash

c=1
# Gravity [m/s2]
for g in 1.62 3.72 9.81 
do
    # Slip [%]
    for s in 20 40 70
    do
        # Wheel load [N]
        for l in 100 164 225
        do
            # Wheel diameter [cm]
            for d in 5 15 30
            do
                # Soil internal friction angles [deg]
                for a in 30 37 43
                do
                    # echo $c $g $s $l $d $a
                    python3 ds_wheel.py $c $g $s $l $d $a
                    c=`expr $c + 1`
                done
            done
        done
    done
done

# Instructions:
# Mark it executable using:
# $ chmod +x run_ds_wheel.sh
# Then run it:
# $ ./run_ds_wheel.sh
# Examples:
# python3 ds_wheel.py 1 9.81 20 100 0.30 30
