---
layout: post
title:  "A Hitchiker's Guide to Spack"
date:   2023-01-30 05:30:00 -0500
categories: "software" 
author_link: "/about.html"
author: Sai Vijay Mocherla
--- 

```
> !!! NOTE: This article is still a WORK UNDER PROGRESS. 
> It mainly started out as a note to myself, to document my learning.  
> I am not sure if I'd setup a comments section to this blog, so until then, 
> please feel free to drop any comments or suggestions over an email.
```
Spack is a package manager for HPC clusters and supercomputers. It makes life easy by letting you install a software packages on your machine without the hassel of conflicting dependencies. With spack you can compile multiple version of a package with different configurations and compilers. Yes you've heard it right, all of this can co-exist in one place!

Here's a simple hitchhiker's guide, mostly for researchers/scientists who want to multiple software packages under one roof but sans the headache of dependency nightmares. Skimming through this tutorial will let you use spack as package manager for your work stations (or even personal machines). I assume that you are running some flavor of linux, but for the most part things from this guide can be extended to MacOS. 

Note: For people who have used conda/homebrew, most of the stuff here will be intuitive but I have tried explain the details so that you can use spack's versatile features.

## 1. Installation 
Preferrably do this installation process as a user and NOT as root. Clone the spack git repository and then install libelf, which is need for spack.
```
$ git clone -c feature.manyFiles=true https://github.com/spack/spack.git
$ cd spack/bin
$ ./spack install libelf
$ . spack/share/spack/setup-env.sh
```
libelf is a C-library to handle ELF objects, its a format for executable files(.x), object code(.o), shared libraries(.so), and core dumps.

## 2. Configuring Compilers
Help spack find compilers already installed on your machine.
```
$ spack compiler find
``` 
This should show something like this :
```
==> Added 2 new compilers to $HOME/user/.spack/linux/compilers.yaml
    gcc@10.3.0  gcc@7.3.0
==> Compilers are defined in the following files:
    $HOME/.spack/linux/compilers.yaml
```
So, I already have two versions of GNU's compilers(this includes gcc, g++ and gfortran) `10.3.0` and `7.3.0`. To check you spack's compiler configurations use:
```
$ spack compilers
``` 
This should show : 
```
==> Available compilers
-- gcc ubuntu21.04-x86_64 ---------------------------------------
gcc@10.3.0  gcc@7.3.0

```
You can see that spack identifies my OS as well as the CPU architecture. One can also manually add the compiler information by editing `$HOME/.spack/linux/compilers.yaml`. The file looks like this:
```
compilers:
- compiler:
    spec: gcc@10.3.0
    paths:
      cc: /usr/bin/gcc
      cxx: /usr/bin/g++
      f77: /usr/bin/gfortran
      fc: /usr/bin/gfortran
    flags: {}
    operating_system: ubuntu21.04
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
- compiler:
    spec: gcc@7.3.0
    paths:
      cc: /home/user/anaconda3/bin/gcc
      cxx: /home/user/anaconda3/bin/g++
      f77: /home/user/anaconda3/bin/gfortran
      fc: /home/user/anaconda3/bin/gfortran
    flags: {}
    operating_system: ubuntu21.04
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []

```
The same can also be accessed using `$ spack config edit compilers` (This should open up  in vi/vim on your terminal.)

## 3. Installing a package 
Lets try to install something lite first. For some reason, let's say I need `python 3.6` an older version of CPython, I can install this using spack without creating problems for any other packages, that are already on my machine.  

Here's the command for that:
```
$ spack install python@3.5 +tkinter ~ssl~sqlite3~dbm target="skylake" %gcc@10.3.0 ^zlib@1.2.12
``` 
A spack package is an installation script, a recipe for compiling/building a particular software. Most scientific software and related packages are readily available in spack, and you can find the package list [here](https://spack.readthedocs.io/en/latest/package_list). If search the package list you'll find that Python has been packaged with details about versions and dependecies, [available here](https://spack.readthedocs.io/en/latest/package_list.html#python). Further, if you check the file [python/package.py](https://github.com/spack/spack/blob/develop/var/spack/repos/builtin/packages/python/package.py) you can find information about different variants. 

Now, lets deconstruct the above command. 

- `@` is used to specify of a package.
- `+` can be used to provide variant information, more on this below.
- `~` will let you generate a minimal build, without some of the optional modules. 
- `%` indicates spack which compilers you want to use. 
- `^` to indicate the particular dependency to be used.

Additional arguments:
- `-d` before install is used run in debug mode.
- `-v` after install for verbose output
- `-j 4` to specify how many processors can be used for building the package.


For example, the `+tkinter` includes tkinter module for python. Though a slightly more time consuming option, including `+optimizations` [as described here](https://github.com/spack/spack/blob/78364a6fe48305c026cdd304cdc28c603d58b54c/var/spack/repos/builtin/packages/python/package.py#L170), will enable build-time optimization for `python =< 3.7`,  where the binaries compiled are upto x10 faster than the normal([see here for more info](https://github.com/docker-library/python/issues/160)). Further, to restrict the build to use `zlib 1.2.12`, a particular version of a dependency that's used for compression of files for python, I add `^zlib@1.2.12` to my installation command.   


### Bonus:
Sine I have an intel i7-7700K CPU on my machine, as another example lets install [Intel's oneAPI compilers](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html). In case, you have an AMD chip on your machine, you can try [AMD's AOCC compilers](https://developer.amd.com/amd-aocc/) that are based of the llvm project.
```
$ spack install intel-oneapi-compilers%gcc@10.3.0
```
Note: spack will now fetch the installation file and try to build from source. In case your internet connection is spotty, the fetching process may timeout. Just rerunning the previous command will get you back on track.  
```
==> Installing intel-oneapi-compilers-2022.1.0-wgdwvl37drmmgi2pz2fpwe3msyldjceh
==> No binary for intel-oneapi-compilers-2022.1.0-wgdwvl37drmmgi2pz2fpwe3msyldjceh found: installing from source
==> Fetching https://registrationcenter-download.intel.com/akdlm/irc_nas/18717/l_dpcpp-cpp-compiler_p_2022.1.0.137_offline.sh
==> Fetching https://registrationcenter-download.intel.com/akdlm/irc_nas/18703/l_fortran-compiler_p_2022.1.0.134_offline.sh
==> Error: timeout: The read operation timed out
```

Once, a package is installed you should see something like this:
```
==> intel-oneapi-compilers: Successfully installed intel-oneapi-compilers-2022.1.0-wgdwvl37drmmgi2pz2fpwe3msyldjceh
  Fetch: 5m 9.89s.  Build: 3m 12.56s.  Total: 8m 22.44s.
[+] /home/dwave/spack/opt/spack/linux-ubuntu21.04-skylake/gcc-10.3.0/intel-oneapi-compilers-2022.1.0-wgdwvl37drmmgi2pz2fpwe3msyldjceh

```
Notice there's some mishmash like `wgdwvl37drmmgi2pz2fpwe3msyldjceh` at the end of your package@version that can easily mistaken to be some garble. But, it is not, its a hash used by spack to keep track of different packages. More on this in the next section.

# 4. Maintaining many different Packages
Manier times in research involving scientific computing, software packages used on a regular basis depend on differing versions of common dependencies. Usually, the developers/maintainers of a software offer support from time-to-time by upgrading their code to be compatible with newer stable versions of their upstream dependencies.  

### Uninstall
To uninstall just do `$ spack uninstall package@version %compiler@version`. There can be more than one package of same version that was compiled with a particular compiler version. 

