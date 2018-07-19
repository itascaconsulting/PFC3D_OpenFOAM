# Using *PFC3D* with *OpenFOAM* for fluid-particle interaction modeling.

The repository contains information about solving fluid-particle
interaction problems by coupling *OpenFOAM* based CFD solvers with
*PFC3D*. *OpenFOAM* is an open source C++ framework for numerical analysis
of continuum mechanics problems. Particle Flow Code in Three
Dimensions (*PFC3D*) is a discrete element code produced by Itasca
Consulting Group.

OPENFOAM® is a registered trade mark of OpenCFD Limited, producer and
distributor of the OpenFOAM software. This offering is not approved or
endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM
software and owner of the OPENFOAM® and OpenCFD® trade marks.

Itasca Consulting Group does not provide support for *OpenFOAM*.

The two-way coarse-grid coupling between *PFC3D* and *OpenFOAM* follows
the approach of Tsuji. Porosity and body force fields are included in
the Navier-Stokes equation to account for the presence of particles in
the flow and these terms are linked to drag forces on the *PFC3D* particles.

A Python module `pyDemFoam` is included which contains modified
versions of the *OpenFOAM* `icoFoam` and `simpleFoam` solvers. These
solvers are modified to account for the presence of solid particles.
Currently, *OpenFOAM* only runs on Linux systems and *PFC3D* only runs
on Windows systems. As a work around *OpenFOAM* is run inside of a
*VirtualBox* Ubuntu guest. Python and TCP sockets are used to link
`demIcoFoam` to *PFC3D*.

The following diagram gives an overview of the system.

![alt text](diagram.png "system schematic")


# Installation

Update to the latest *PFC3D* version: http://www.itascacg.com

## Setting up *VirtualBox*

Install VirtualBox 5.0.20 or newer http://download.virtualbox.org/virtualbox/5.0.20/VirtualBox-5.0.20-106931-Win.exe

Install Ubuntu 16.04 into VirtualBox http://www.ubuntu.com/download/desktop

Install the VirtualBox Linux Guest Additions: In the Virtual Box
window select Devices > Insert the Guest Additions CD Image. Follow
the on screen instructions.

Make sure you Ubuntu environment is current with:

`sudo apt-get update`

and

`sudo apt-get upgrade`


```bash
sudo apt-get install python-pip emacs24 git gitk build-essential cmake flex bison zlib1g-dev qt4-dev-tools libqt4-dev libqtwebkit-dev gnuplot libreadline-dev libncurses5-dev libxt-dev libopenmpi-dev openmpi-bin libboost-system-dev libboost-thread-dev libgmp-dev libmpfr-dev python python-dev libcgal-dev python-numpy ipython python-scipy cython
```

Install the `itasca` Python module with `pip install itasca`


## Building *OpenFOAM* v3.0+

What follows here is a distillation of the instructions from here:

https://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM%2B-3.0%2B/Ubuntu

consult this link if you have problems.

First, download and unzip OpenFOAM v3.0+ and the third party
dependencies with the following commands.

```bash
cd ~
mkdir OpenFOAM
cd OpenFOAM
wget "http://downloads.sourceforge.net/openfoamplus/files/OpenFOAM-v3.0%2B.tgz?use_mirror=mesh" -O OpenFOAM-v3.0+.tgz
wget "http://downloads.sourceforge.net/openfoamplus/files/ThirdParty-v3.0%2B.tgz?use_mirror=mesh" -O ThirdParty-v3.0+.tgz
tar -xzf OpenFOAM-v3.0+.tgz
tar -xzf ThirdParty-v3.0+.tgz
```

Add the following line to the end of: `~/.bashrc`

`source $HOME/OpenFOAM/OpenFOAM-v3.0+/etc/bashrc WM_NCOMPPROCS=4 WM_MPLIB=SYSTEMOPENMPI`

To load these changes:
`source ~/.bashrc`


### Building *ParaView*

*ParaView* is a post-processing program which can be used to visualize
*OpenFOAM* results. It it typically used via the `paraFoam` wrapper.

```bash
cd $WM_THIRD_PARTY_DIR
export QT_SELECT=qt4
./Allwmake > log.make 2>&1
wmSET $FOAM_SETTINGS

./makeParaView4 -python -mpi -python-lib /usr/lib/x86_64-linux-gnu/libpython2.7.so.1.0 > log.makepv4_2

```

If you have problems with paraFoam try: `export LIBGL_ALWAYS_SOFTWARE=1`

```bash
wmSET $FOAM_SETTINGS
cd $WM_PROJECT_DIR
export QT_SELECT=qt4

find src applications -name "*.L" -type f | xargs sed -i -e 's=\(YY\_FLEX\_SUBMINOR\_VERSION\)=YY_FLEX_MINOR_VERSION < 6 \&\& \1='

./Allwmake > log.make 2>&1
```

At this point *OpenFOAM* should be built. Try `icoFoam -help` to
confirm.

## Building `pyDemFoam`

`pyDemFoam` is a Python module which contains modified versions of
*OpenFOAM* solvers. It includes porosity and body force terms to
account for the presence of the particles. Modified version of both
`icoFoam` and `simpleFoam` are provided.

```bash
cd
mkdir src
cd src
git clone https://github.com/jkfurtney/PFC3D_OpenFOAM.git

cd PFC3D_OpenFOAM/pyDemFoam/
python setup.py install --user
```

Test that this worked:

```bash
cd ~
python -c "import pyDemFoam; print pyDemFoam.__version__"
```

A version number number like: `2016.06.09` should be shown if the
installation worked correctly.

# Running coupled problems

Under Windows make a clone of this repository.

A demonstration and verification problem is included in the
`dropTest1/` folder.

#### In *PFC3D*

- Start *PFC3D* and open `dropTest1/dropTest1.p3prj`

- Run the file `pfc_dropTest1.py`

#### In Ubuntu:

```bash
cd ~/src/PFC3D_OpenFOAM/dropTest1/
blockMesh
python cfd_dropTest1.py
```

This should launch the coupled calculation.

![alt text](dropTest1/dropTest1.png "Model Results")

# Limitations

This work is intended as a demonstration of how to connect *PFC3D* to
a CFD solver. The implementation given here is limited in the
following ways:

- No linear relaxation is used to stabilize the equations. Numerical
  instabilities are likely to occur.

- No turbulence model is included in the analysis.

- Time derivatives of porosity are not included in the momentum or continuity equations.

- Further verification of the coupled system is underway.

# More Information

Documentation for *OpenFOAM* can be found here:
http://cfd.direct/openfoam/documentation/

## Building a local copy of the *OpenFOAM* documentation

If you are exploring the *OpenFOAM* C++ code it is helpful to have a
local copy of the OpenFOAM C++ API documentation. GNU Global can help
navigating large source codes.

```bash
sudo apt-get install doxygen graphviz global
cd $WM_PROJECT_DIR/doc
./Allwmake
```

The Doxygen reference should be in `$WM_PROJECT_DIR/doc/Doxygen/html/index.html`

python setup.py install --user > log.txt 2>&1 && tail log.txt

## Other projects

* OpenFOAM LAMMPS coupling: https://github.com/xiaoh/sediFoam

* OpenFOAM YADE coupling: http://trace.tennessee.edu/utk_graddiss/21/
