Installation
============

Stable version
--------------

Stable versions of precellar are published on PyPI.
Precompiled binaries are available for x86_64 Linux systems.
So installing it is as simple as running:

```
pip install precellar
```

You may need the [Rust](https://www.rust-lang.org/tools/install) compiler during
the installation process. The Rust compiler can be installed using:

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Nightly build
-------------

The nightly build is the build from the latest source codes, which includes the
latest features, enhancements, and bug fixes that haven't been released. 
The nightly build can be unstable and include some untested features.

The [nightly release](https://github.com/regulatory-genomics/precellar/releases/tag/nightly) page
contains wheel files for the nightly build.
Please download the corresponding wheel file for your platform and use `pip install` to install it.
For example, if you are using a Linux system with Python 3.8, you can use the following command to install it:

```
pip install precellar-x.x.x-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
```

Build from the latest source code 
---------------------------------

Building the precellar library requires the [Rust](https://www.rust-lang.org/tools/install) compiler.
The Rust compiler can be installed using:

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Once you have the Rust compiler properly installed, you can use pip to install the precellar library:

```
pip install 'git+https://github.com/regulatory-genomics/precellar.git#egg=precellar&subdirectory=python'
```