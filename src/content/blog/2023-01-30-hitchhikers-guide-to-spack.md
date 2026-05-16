---
title: "A Hitchhiker's Guide to Spack"
description: "A short, practical guide to Spack — a package manager for HPC clusters and scientific workstations — that I wrote as a note to myself while learning it."
date: 2023-01-30
tags: ["software", "hpc", "tutorial"]
---

> **Note:** This article started as a note-to-self while I was learning Spack.
> It is still a work in progress. Comments and corrections welcome by
> [email](mailto:vijaysai.mocherla@gmail.com).

[Spack](https://spack.io/) is a package manager for HPC clusters and supercomputers.
It makes life easier by letting you install software on your machine without the
hassle of conflicting dependencies. With Spack you can compile multiple versions
of a package with different configurations and compilers — and all of it can
co-exist in one place.

This guide is mostly aimed at researchers and scientists who want to keep many
software packages under one roof but without the dependency nightmares. Skimming
through this will let you use Spack as a package manager for your workstations
(or even personal machines). I assume you are running some flavour of Linux,
but most of this carries over to macOS. If you have used conda or Homebrew, most
of the concepts will be intuitive.

## 1. Installation

Preferably install Spack as a regular user, **not** as root. Clone the Spack git
repository and install `libelf`, which Spack needs:

```bash
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
cd spack/bin
./spack install libelf
. spack/share/spack/setup-env.sh
```

`libelf` is a C library for handling ELF objects — the format used for
executables (`.x`), object code (`.o`), shared libraries (`.so`), and core dumps.

## 2. Configuring compilers

Help Spack find compilers already installed on your machine:

```bash
spack compiler find
```

You should see something like:

```
==> Added 2 new compilers to $HOME/user/.spack/linux/compilers.yaml
    gcc@10.3.0  gcc@7.3.0
==> Compilers are defined in the following files:
    $HOME/.spack/linux/compilers.yaml
```

So I already have two versions of GCC (10.3.0 and 7.3.0). To list Spack's
compiler configurations:

```bash
spack compilers
```

```
==> Available compilers
-- gcc ubuntu21.04-x86_64 ---------------------------------------
gcc@10.3.0  gcc@7.3.0
```

Spack identifies the OS and CPU architecture for you. You can also add compilers
manually by editing `$HOME/.spack/linux/compilers.yaml`:

```yaml
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
```

The same file can be opened with `spack config edit compilers`.

## 3. Installing a package

Let us install something light first. Say I need Python 3.5 (an older version of
CPython) — I can install it with Spack without disturbing anything else on the
machine:

```bash
spack install python@3.5 +tkinter ~ssl~sqlite3~dbm target="skylake" %gcc@10.3.0 ^zlib@1.2.12
```

A Spack package is an installation recipe for building a particular piece of
software. Most scientific software is already packaged — see the
[package list](https://spack.readthedocs.io/en/latest/package_list.html).
The [`python/package.py`](https://github.com/spack/spack/blob/develop/var/spack/repos/builtin/packages/python/package.py)
file is a good example of variants and dependencies.

Deconstructing the command:

- `@` specifies the version
- `+` enables a variant
- `~` disables a variant (minimal build)
- `%` selects which compiler to use
- `^` pins a dependency to a particular spec

Useful extra flags:

- `-d` before `install` enables debug mode
- `-v` after `install` for verbose output
- `-j 4` controls the parallelism of the build

For example, `+tkinter` includes the `tkinter` module for Python. The
`+optimizations` variant
([source](https://github.com/spack/spack/blob/78364a6fe48305c026cdd304cdc28c603d58b54c/var/spack/repos/builtin/packages/python/package.py#L170))
enables build-time optimization for Python ≤ 3.7 — the resulting binaries are
up to ~10× faster
([context](https://github.com/docker-library/python/issues/160)). To pin a
particular `zlib` version as a dependency, I add `^zlib@1.2.12` to the command.

### Bonus: Intel oneAPI

Since I have an Intel i7-7700K, I can install
[Intel's oneAPI compilers](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)
with:

```bash
spack install intel-oneapi-compilers%gcc@10.3.0
```

(If you are on AMD, [AOCC](https://developer.amd.com/amd-aocc/) is the
equivalent.)

Spack will fetch the installer and build from source. On a flaky connection the
fetch may time out — just re-run the command and Spack will resume:

```
==> Installing intel-oneapi-compilers-2022.1.0-...
==> No binary for intel-oneapi-compilers-... found: installing from source
==> Fetching https://registrationcenter-download.intel.com/akdlm/irc_nas/18717/l_dpcpp-cpp-compiler_p_2022.1.0.137_offline.sh
==> Error: timeout: The read operation timed out
```

A successful install ends with:

```
==> intel-oneapi-compilers: Successfully installed intel-oneapi-compilers-2022.1.0-wgdwvl37drmmgi2pz2fpwe3msyldjceh
  Fetch: 5m 9.89s.  Build: 3m 12.56s.  Total: 8m 22.44s.
[+] /home/user/spack/opt/spack/linux-ubuntu21.04-skylake/gcc-10.3.0/intel-oneapi-compilers-2022.1.0-wgdwvl37drmmgi2pz2fpwe3msyldjceh
```

The garbled-looking suffix (`wgdwvl37...`) is not random — it is the Spack hash
used to uniquely identify the build.

## 4. Maintaining many different packages

In scientific computing, the software packages you use day-to-day often depend
on slightly different versions of common libraries. Spack handles this by
keeping each build addressable by its hash, so you can have many parallel
installs without conflicts.

### Uninstall

To uninstall a build:

```bash
spack uninstall package@version %compiler@version
```

There can be more than one package of the same version compiled with the same
compiler — in that case Spack will ask you to be more specific by hash.

---

*This post is a work in progress and will keep growing as I cover environments,
modules, and recipes for building common scientific stacks.*
