# Modeling rigid body-driven granular flows <br /> via Material Point Method with Nonlocal Granular Fluidity
#### [[Introductory Video](https://youtu.be/_iNQWDR8nNA)] [[Paper 1](https://www.sciencedirect.com/science/article/pii/S0045782522001876)] [[Paper 2](https://www.iaarc.org/publications/2020_proceedings_of_the_37th_isarc/efficient_numerical_methods_for_accurate_modeling_of_soil_cutting_operations.html)] [[Paper 3](https://link.springer.com/article/10.1007/s10035-022-01228-4)]

This is a modified version of Taichi MPM written by [Y. Hu (MIT)](https://github.com/yuanming-hu/taichi_mpm) for modeling granular flows via the Nonlocal Granular Fluidity constitutive model developed by Prof. [K. Kamrin (MIT)](http://web.mit.edu/kkamrin/www/fluid_solid.html) and their interactions with rigid bodies. The numerical approach is developed under the supervision of Prof. [K. Skonieczny (Concordia)](http://users.encs.concordia.ca/~kskoniec/).

<img src="https://github.com/haeriamin/files/blob/master/indExcav.gif" alt="drawing" width="400"> <img src="https://github.com/haeriamin/files/blob/master/wheelSoil.gif" alt="drawing" width="400"> <img src="https://github.com/haeriamin/files/blob/master/explate1.gif" alt="drawing" width="400"> <img src="https://github.com/haeriamin/files/blob/master/explate2exp.gif" alt="drawing" width="400"> <img src="https://github.com/haeriamin/files/blob/master/silo.gif" alt="drawing" width="400"> <img src="https://github.com/haeriamin/files/blob/master/tcFlow.gif" alt="drawing" width="400">


# Install and Run Demo
Tested on Ubuntu 16.04 and 18.04:

1. Run `cd ~ ; wget https://github.com/haeriamin/MPM-NGF/raw/master/modifications/install.sh && chmod +x ./install.sh ; ./install.sh`

2. Run `python3 ~/taichi/projects/mpm/scripts/excav.py`

3. Download and install [Houdini Apprentice](https://www.sidefx.com/products/houdini-apprentice) for visualization (optional)

More explanations can be found [here](https://github.com/yuanming-hu/taichi_mpm#particle-attributes).


# Bibtex
Please cite our papers [[1](https://www.sciencedirect.com/science/article/pii/S0045782522001876), [2](https://www.iaarc.org/publications/fulltext/ISARC_2020_Paper_60.pdf), [3](https://link.springer.com/article/10.1007/s10035-022-01228-4)] if you use this code for your research: 
```
@article{HAERI2022,
	title = {Three-dimensionsal granular flow continuum modeling via material point method with hyperelastic nonlocal granular fluidity},
	journal = {Computer Methods in Applied Mechanics and Engineering},
	author = {Amin Haeri and Krzysztof Skonieczny},
	volume = {394},
	pages = {114904},
	year = {2022},
	issn = {0045-7825},
	doi = {https://doi.org/10.1016/j.cma.2022.114904},
	url = {https://www.sciencedirect.com/science/article/pii/S0045782522001876}
}
```
and/or
```
@inproceedings{HAERI2020,
	title = {Efficient Numerical Methods for Accurate Modeling of Soil Cutting Operations},
	author = {Haeri, Amin and Tremblay, Dominique and Skonieczny, Krzysztof and Holz, Daniel and Teichmann, Marek},
	booktitle = {Proceedings of the 37th International Symposium on Automation and Robotics in Construction (ISARC)},
	isbn = {978-952-94-3634-7},
	publisher = {International Association for Automation and Robotics in Construction (IAARC)},
	editor = {"Osumi, Hisashi", "Furuya, Hiroshi", "Tateyama, Kazuyoshi"},
	doi = {10.22260/ISARC2020/0085},
	year = {2020},
	month = {October},
	pages = {608-615},
	address = {Kitakyushu, Japan}
}
```
and/or
```
@article{HAERI2022Gravity,
	title = {Gravity sensitivity of continuum numerical solvers for granular flow modeling},
	journal = {Granular Matter},
	author = {Amin Haeri and Krzysztof Skonieczny},
	volume = {24},
	year = {2022},
	issue = {2},
	doi = {https://doi.org/10.1007/s10035-022-01228-4},
	url = {https://link.springer.com/article/10.1007/s10035-022-01228-4}
}
```