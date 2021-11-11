# AMIGOS III #
Pseudo-torsion angle visualization and motif-based structure comparison of
nucleic acids

## Dependency ##
AMIGOS III is a plugin to [PyMOL](https://pymol.org/) version 2.5 or later.
It requires matplotlib and python3.
* Not all official Windows build of PyMOL by Shrodinger include matplotlib.
  If AMIGOS III cannot be launched after installation due to lack of
  matplotlib, you re-install 
  [Open-Source Pymol on Windows](https://pymolwiki.org/index.php/Windows_Install)
  where you can install matplotlib.

## Automated Installation ##
If you use 64bit Windows, Linux or Mac, download the zip files for your
operating system at 
[Releases](https://github.com/pylelab/AMIGOSIII/releases)
and install the zip file through the
[PyMOL Plugin Manager](https://pymolwiki.org/index.php/Plugins)

## Manual Installation ##
AMIGOS III includes a C++ program ``NaTorsion``, which can be compiled by
```
make
```
or
```
g++ -O3 NaTorsion.cpp -o NaTorsion -static
```
The "-static" flag should be removed on Mac OS, which does not support static
executable.

After the compilation is finished, compress the ``AMIGOIII/`` folder (i.e.,
this folder) as a zip file, and install the zip file through the
[PyMOL Plugin Manager](https://pymolwiki.org/index.php/Plugins).
Alternatively, you can manually move the ``AMIGOIII/`` folder to
``pmg_tk/AMIGOSIII`` at your pymol installation path.

## Reference ##
Morgan Shine, Chengxin Zhang, Anna Marie Pyle (2021)
"AMIGOS III: Pseudo-torsion angle visualization and motif-based structure
comparison of nucleic acids"
