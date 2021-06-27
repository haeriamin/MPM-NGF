#!/bin/bash

set -e  # Fail on any error
set -x  # Display commands being run

cd ~
sudo apt install gcc python3-pip
sudo apt-get install -y python3-dev git build-essential cmake make g++ libx11-dev python3-pyqt5 libtbb2 libtbb-dev

wget https://raw.githubusercontent.com/yuanming-hu/taichi/legacy/install.py
sudo -H python3 install.py

source ~/.bashrc

cd ~/taichi/python/taichi/
wget -O https://github.com/haeriamin/MPM-NGF/blob/master/modifications/main.py

cd ~/taichi/include/taichi/dynamics/
wget -O https://github.com/haeriamin/MPM-NGF/blob/master/modifications/rigid_body.py

cd ~/taichi/python/taichi/dynamics/
wget -O https://github.com/haeriamin/MPM-NGF/blob/master/modifications/mpm.py

cd ~
source ~/.bashrc
ti install mpm
ti build