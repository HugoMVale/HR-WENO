---
project: hrweno
license: mit
summary: hrweno is a modern-Fortran implementation of selected WENO and TVD integration schemes.
src_dir: ./src
         ./example
output_dir: _site
page_dir: ./doc
source: true
proc_internals: true
graph: true
coloured_edges: true
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
project_github: https://github.com/HugoMVale/HR-WENO
author: HugoMVale
github: https://github.com/HugoMVale/
email: 57530119+HugoMVale@users.noreply.github.com
dbg: true
predocmark: >
docmark_alt: #
predocmark_alt: <
md_extensions: markdown.extensions.toc
---

About
=====

This package is a modern-Fortran implementation of selected high-resolution [weighted essentially non-oscillatory (WENO)](https://en.wikipedia.org/wiki/WENO_methods) schemes and [total variation diminishing (TVD)](https://en.wikipedia.org/wiki/Total_variation_diminishing) integration methods for solving [hyperbolic conservation equations](https://en.wikipedia.org/wiki/Hyperbolic_partial_differential_equation).

In particular, the package includes WENO schemes up to 5th order, explicit Runge-Kutta TVD methods up to 3rd order, and an explicit 3rd order multi-step TVD method.  

All numerical methods are described in detail by [Shu (1997)](HR-WENO/doc/Shu-WENO-notes.pdf).
