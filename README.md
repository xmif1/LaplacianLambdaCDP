# Laplacian (Optimal) Lambda CDP Finder

---
Xandru Mifsud (2021)

_Please cite as_: Mifsud X., "\[Lambda]--Core distance partitions", Linear Algebra Appl. (2021), https://doi.org/10.1016/j.laa.2020.12.012

---

## Requirements

This program is intended to run on Linux platforms, in particular on Debian-based systems such as Ubuntu and
Lubuntu, on which we have tested the implementation with relative success.

For installation, ```cmake``` version 3.17+ is required, as well as ```Eigen``` version 3.3 and above.

## Installation Instructions

Clone the repository, and ```cd``` into the project directory. Then run:

1. ```cmake .```
2. ```make```

## Execution Instructions

Simply ```cd``` into the directory containing the compiled executable, and run ```./LaplacianLambdaCDP <N> <path/to/adjacency_matrix.csv>```,
where ```N``` is a required argument specifying the dimension of the adjacency matrix, represented in ```CSV``` format and
specified by the required filepath ```path/to/adjacency_matrix.csv```.