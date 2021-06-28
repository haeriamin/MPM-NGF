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
sudo wget -O main.py https://github.com/haeriamin/MPM-NGF/raw/master/modifications/main.py

cd ~/taichi/include/taichi/dynamics/
sudo wget -O rigid_body.h https://github.com/haeriamin/MPM-NGF/raw/master/modifications/rigid_body.h

cd ~/taichi/python/taichi/dynamics/
sudo wget -O mpm.py https://github.com/haeriamin/MPM-NGF/raw/master/modifications/mpm.py

cd ~
source ~/.bashrc
ti install mpm
ti build
