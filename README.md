# AMIGOS III #
Pseudo-torsion angle visualization and motif-based structure comparison of nucleic acids

## Dependency ##
AMIGOS III is a plugin to [PyMOL](https://pymol.org/) version 2.5 or later.

## Installation ##
AMIGOS III includes a C++ program ``NaTorsion``, which can be compiled by
```
make
```
or
```
g++ -O3 NaTorsion.cpp -o NaTorsion -static
```
The "-static" flag should be removed on Mac OS, which does not support static executable.

After the compilation is finished, compress the AMIGOIII folder as a zip file, and
install the zip file as through the [PyMOL Plugin Manager](https://pymolwiki.org/index.php/Plugins)

## Reference ##
Morgan Shine, Chengxin Zhang, Anna Marie Pyle (2021)
"AMIGOS III: Pseudo-torsion angle visualization and motif-based structure comparison of nucleic acids"
