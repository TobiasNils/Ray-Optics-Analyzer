# Ray-optics Analyzer

Integration of ray-optics simulations for optical systems design in [Blender](https://www.blender.org/).
The raytracing algorithms implemented here are adaptations of the [pyOpTools](https://github.com/cihologramas/pyoptools) v0.1.1 code-base, the addon structure and visualization callblacks are based on [Math Vis (Console)](https://wiki.blender.org/wiki/Extensions:2.6/Py/Scripts/3D_interaction/Math_Viz). 

---

>## pyOpTools

>pyOpTools is a set of packages that allow the simulation of optical systems by raytracing as well as some calculations involving wavefronts, currently under development. It is written in Python and Cython, and is being developed by the technological development group of [Combustión Ingenieros S.A.S](http://www.cihologramas.com), and the applied optics group of the Universidad Nacional de Colombia.

>The pyOpTools is divided in several packages that contain the different library's functionalities:

>    pyoptools.raytrace

>This package contains the classes and functions used to perform the simulation of optical systems using 3D non-sequential raytracing algorithms.

>    pyoptools.misc

>This package contains miscellaneous classes and functions used by the other packages, but that can not be classified in any of them.

---

## Requirements
(Based on experience on Linux Fedora Workstation 31)

* `pip3 install Cython --user`
* `sudo dnf install python3-devel` (`python3-dev` in debian)
* `cd path/to/this/repository/ && python3 setup.py install --user`
* `pip3 install scipy --user` <- check if really necessary for this implementation
* `pip3 install chaospy --user`


## Todo

- [ ] implement and benchmark multiprocessing for System.propagate()
- [ ] use node-editor for material setup?
