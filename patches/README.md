# Patches

[msgsteiner](http://areeweb.polito.it/ricerca/cmp/code/bpsteiner) does not compile on many 
systems with the default Makefile provided in its source distribution.

In this directory we have prepared patches for various systems to resolve common compilation issues.

If compilation of msgsteiner fails on your system, download the patch file associated with 
your operating system, place it in the directory with the msgsteiner source code, run the `patch` command,
and try the build again with `make`:

```bash
cd msgsteiner-1.3
wget "https://raw.githubusercontent.com/fraenkel-lab/OmicsIntegrator/master/patches/Makefile.linux.patch"
patch Makefile Makefile.linux.patch
make
```
