project: Musubi
summary: A multi-level parallel lattice Boltzmann solver within the APES suite.
src_dir: source
src_dir: build/ford
exclude_dir: build/ford/treelm
external: aoturl = https://geb.inf.tu-dresden.de/doxy/aotus
external: temurl = https://geb.inf.tu-dresden.de/doxy/treelm
output_dir: docu
page_dir: doc_pages
media_dir: examples/tutorials/images
copy_subdir: media
project_website: https://geb.inf.tu-dresden.de/pages/musubi.html
project_github: https://github.com/apes-suite/musubi.git
graph: true
graph_maxdepth: 4
graph_maxnodes: 32
display: public
display: protected
display: private
sort: permission
source: false
externalize: true
author: University of Siegen
title: Mainpage
print_creation_date: True
md_extensions: markdown.extensions.toc

Introduction to Musubi
======================

![Apes Logo](|page|/apes_sub_small.png)


Musubi is the multi-level parallel lattice Boltzmann solver within the APES
suite.
It is working on a linearized octree and uses efficient data structures
allowing adaptive parallel simulations.
Musubi offers several collision kernels and is designed in a way to deal with
huge meshes efficiently.

See the [documentation](|page|) for details, tutorials and examples.
