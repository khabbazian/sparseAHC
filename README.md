
An implementation of the agglomerative hierarchical clustering method that accepts sparse similarity matrices. 
A similarity matrix can be interpreted as the adjacency matrix of a weighted undirected graph. 
In this graph, the edge weights represent the similarity of the two endpoints. 
Furthermore, absence of an edge is weak signal of similarity. 
This implementation builds the output dendrograms without expanding the matrix into the memory. 
The memory complexity is linear in terms of the number of non-zero entries in the input. 

[SparseAHC Manual](http://homepages.cae.wisc.edu/~khabbazian/pdfs/sparseAHC.pdf)

### Install using the devtools package.
```
install.packages("devtools")
require(devtools)
install_github("khabbazian/sparseAHC")
require(sparseAHC)
```

### License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<http://www.r-project.org/Licenses/GPL-3>
