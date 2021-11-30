# AMIGOS III #
Pseudo-torsion angle visualization and motif-based structure comparison of
nucleic acids

## Dependency ##
AMIGOS III is a plugin to [PyMOL](https://pymol.org/) version 2.5 or later.
It requires matplotlib, pandas and python3.

## Installation ##
#### Official PyMOL build ####
If you use an official PyMOL build for Shrondinger, you will need to first
check if matplotlib and pandas are available. Launch PyMOL and type in the
following command in the PyMOL command prompt:
```python
import matplotlib
import pandas
```
If the above command complaint that
``ModuleNotFoundError: No module named 'matplotlib'`` or
``ModuleNotFoundError: No module named 'pandas'`` or
use the following command to install matplotlib and pandas respectively:
```python
conda install matplotlib
conda install pandas
```
You can then download the zip file of AMIGOS III plugin for your operating
system (Mac, Windows or Linux) from 
[Releases](https://github.com/pylelab/AMIGOSIII/releases)
and install the zip file through the
[PyMOL Plugin Manager](https://pymolwiki.org/index.php/Plugins)

#### Open-Source PyMOL on Windows ####
If you use Christoph Gohlke's
[Open-Source Pymol on Windows](https://pymolwiki.org/index.php/Windows_Install),
you will need to install download pre-compiled wheel files for
[matplotlib](https://www.lfd.uci.edu/~gohlke/pythonlibs/#matplotlib) and
[pandas](https://www.lfd.uci.edu/~gohlke/pythonlibs/#pandas)
using ``pip install``. After matplotlib and pymol are install, you can install
the zip file for AMIGOS III plugin for Windows through the
[PyMOL Plugin Manager](https://pymolwiki.org/index.php/Plugins)
as mentioned above.

#### PyMOl for a Linux distribution ####
You can also install AMIGOS III to the PyMOL distributed by your Linux distro.
Install ``pymol``, ``pandas`` and ``matplotlib`` for python3 using the package
manager of your distro . Download the zip file for Linux from
[Releases](https://github.com/pylelab/AMIGOSIII/releases)
Launch pymol using root:
```bash
sudo pymol
```
Install the zip file through the
[PyMOL Plugin Manager](https://pymolwiki.org/index.php/Plugins).
After the installation is completed, you can run PyMOL as a normal user
without ``sudo``.

### Other systems ###
If you are using an operating system other than Windows, Mac and Linux, or if
your operating system is not 64bit, you will need to manually install AMIGOS III.
AMIGOS III includes a C++ program ``NaTorsion``, which can be compiled by
```
make
```
or
```
g++ -O3 NaTorsion.cpp -o NaTorsion -static
```
The "-static" flag should be removed on Mac OS, which does not support static
executable. After the compilation is finished, manually move the ``AMIGOIII/``
folder (i.e., the folder containing this README.md file) to
``pmg_tk/AMIGOSIII`` at your pymol installation path.

## Reference ##
Morgan Shine, Chengxin Zhang, Anna Marie Pyle (2021)
"AMIGOS III: Pseudo-torsion angle visualization and motif-based structure
comparison of nucleic acids"
