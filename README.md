# Using *PFC3D* v700 with *OpenFOAM* v2112 for fluid-particle interaction modeling


The repository contains information about solving fluid-particle
interaction problems by coupling *OpenFOAM* based CFD solvers with
*PFC3D*. *OpenFOAM* is an open source C++ framework for numerical
analysis of continuum mechanics problems. Particle Flow Code in Three
Dimensions (*PFC3D*) is a discrete element code produced by Itasca
Consulting Group.

OPENFOAM® is a registered trade mark of OpenCFD Limited, producer and
distributor of the OpenFOAM software. This offering is not approved or
endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM
software and owner of the OPENFOAM® and OpenCFD® trade marks.

Itasca Consulting Group does not provide support for *OpenFOAM*.

The two-way coarse-grid coupling between *PFC3D* and *OpenFOAM*
follows the approach of Tsuji. Porosity and body force fields are
included in the Navier-Stokes equations to account for the presence of
particles in the flow and these terms are linked to drag forces on the
*PFC3D* particles.

A Python module `pyDemFoam` is included which contains modified
versions of the *OpenFOAM* `icoFoam` solver. These solvers are
modified to account for the presence of solid particles. This package
works for only the Ubuntu 20.04 LTS Linux distribution.

# Installation

These instructions work for the x86_64 version of Ubuntu 20.04 LTS.

Make sure your Ubuntu system is updated:
```bash
sudo apt-get update && sudo apt-get upgrade
```

Update to the latest Linux version of *PFC3D* 7.00: https://www.itascacg.com/software/downloads/itasca-linux-software-7-0-update

```bash
wget https://itasca-software.s3.amazonaws.com/itasca-software/v700/itascasoftware_700.145.deb
sudo DEBIAN_FRONTEND=noninteractive apt-get -y install -f ./itascasoftware_700.146.deb
mkdir -p ~/.config/Itasca
touch ~/.config/Itasca/wad700.conf
echo "[weblicense]" >> ~/.config/Itasca/wad700.conf
echo "email=your email here" >> ~/.config/Itasca/wad700.conf
echo "password=your web license password here" >> ~/.config/Itasca/wad700.conf
```

Run the command line version of PFC3D v700
```bash
pfc3d700_console
```

Give the command `license list web` to make sure your web license is
working. Type `exit` to leave *PFC3D*. At this point, *PFC3D* should
be installed and working and we move into building OpenFOAM and the
coupling module.

Install some prerequisite packages:
```bash
sudo apt-get install python-is-python3 python3-numpy ipython3 cython3 python3-dev python3-pip git-core build-essential cmake libfl-dev bison zlib1g-dev qttools5-dev qtbase5-dev libqt5x11extras5-dev gnuplot libreadline-dev libncurses-dev libxt-dev libopenmpi-dev openmpi-bin libboost-system-dev libboost-thread-dev libgmp-dev libmpfr-dev libcgal-dev curl libglu1-mesa-dev
```

Install the `itasca` Python module with `pip install itasca`

## Building *OpenFOAM* v2112

This package and these instructions only work with OpenFOAM v2112.

What follows here is a distillation of the instructions from here:

https://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-7/Ubuntu/20.04

consult this link if you have problems.

First, download and unzip OpenFOAM v2112 with the following commands:

```bash
cd ~
mkdir OpenFOAM
cd OpenFOAM
wget https://dl.openfoam.com/source/v2112/OpenFOAM-v2112.tgz
tar -xzf OpenFOAM-v2112.tgz
ln -s /usr/bin/mpicc.openmpi OpenFOAM-v2112/bin/mpicc
ln -s /usr/bin/mpirun.openmpi OpenFOAM-v2112/bin/mpirun
```

These commands add the OpenFOAM shell environment settings needed to
build and use OpenFOAM. If you have a problem later make sure you did
this step correctly.

```bash
echo "source ~/OpenFOAM/OpenFOAM-v2112/etc/bashrc WM_LABEL_SIZE=64 FOAMY_HEX_MESH=yes" >> ~/.bashrc
. ~/.bashrc
```


```bash
cd $WM_PROJECT_DIR
export QT_SELECT=qt5

./Allwmake -j 12 > log.make 2>&1
```

At this point *OpenFOAM* should be built. Try `icoFoam -help` to
confirm.

## Building `pyDemFoam`

`pyDemFoam` is a Python module which contains modified versions of
*OpenFOAM* solvers. It includes porosity and body force terms to
account for the presence of the particles. Two versions `icoFoam` are
provided. `pyDemIcoFoam` uses the explicit formulation for drag
described in the PFC3D manual in the cfd module section and
`pyDemIcoFoamSemiImplicitDrag` uses the semi implicit drag treatment
described in this paper: Xiao, H., Sun, J., 2010. Algorithms in a
Robust Hybrid CFD-DEM Solver for Particle-Laden Flows. Commun. Comp.
Phys.

```bash
cd
mkdir src
cd src
git clone https://github.com/jkfurtney/PFC3D_OpenFOAM.git

cd PFC3D_OpenFOAM/pyDemFoam/
python setup.py build
```
You will see hundreds of compiler warnings after this step, as long as
there are no error messages everything should work fine. Install the
pyDemFoam Python module:
```bash
python setup.py install --user
```
Test that this worked:

```bash
cd ~
python -c "import pyDemFoam; print(pyDemFoam.__version__)"
```

A version number number like: `2022.01.09` should be shown if the
installation worked correctly.

# Running coupled problems

A demonstration and verification problem is included in the
`dropTest1/` folder.

Open two terminals on the same computer. Change directories to the
`dropTest1/` folder in each terminal. In the first, give the command:

`pfc3d700_console call pfc_dropTest1.py`

and in the other give the command:

`python cfd_dropTest1.py`

The coupled problems should run.

![alt text](dropTest1/dropTest1.png "Model Results")

Only the dropTest1/ and porous1/ examples are working.

# Limitations

This work is intended as a demonstration of how to connect *PFC3D* to
a CFD solver. The implementation given here is limited in the
following ways:

- No linear relaxation is used to stabilize the equations. Numerical
  instabilities are likely to occur.

- No turbulence model is included in the analysis.

- Time derivatives of porosity are not included in the momentum or continuity equations.


# More Information

Documentation for *OpenFOAM* can be found here:
https://www.openfoam.com/documentation/overview

# Building `pyDemFoam` for Python3.6
PFC3D v700 has an internal Python 3.6 interpreter, the pyDemFoam
module can be compiled against Python 3.6 so it can be loaded directly
into PFC3D 700.


```bash
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.6
sudo apt install python3.6-dev

cd pyDemFoam
python3.6 setup.py build
python3.6 setup.py install --user
# manually install these
sudo cp ~/.local/lib/python3.6/site-packages/_pyDemFoam.cpython-36m-x86_64-linux-gnu.so /opt/itascasoftware/v700/python36/lib/python3.6/site-packages/
sudo cp -r ~/.local/lib/python3.6/site-packages/pyDemFoam/ /opt/itascasoftware/v700/python36/lib/python3.6/site-packages/

```
