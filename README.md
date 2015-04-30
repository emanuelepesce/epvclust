epvclust
===========================

##About
"epvclust" is an optimized version of R "pvclust" package (http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/pvclust).
It performs the same work of pvclust for assessing the uncertainty in hierarchical cluster analysis, but with better performance (mostly on time performance).

It's available an offline installable R package.
It's also available the C code for doing multiscale bootstrap.


##About Pvclust
- [Official page] (http://www.sigmath.es.osaka-u.ac.jp/shimo-lab/prog/pvclust/)
- [Theory about multiscale bootstrap] (http://www.is.titech.ac.jp/~shimo/multiboot.html)

##Contents

### R
This folder contains the R epvclust package. For install it execute the file "install.R" with RStudio ("Devtools" R package is a dependency).

"install.R" contains also a runnable example.

### C 
This folder contains the C code for computing multiscale bootstrap.
Main file: "pvclust.c".

### exampleData
This folder contains various dataset with different sizes.

##Authors
[Emanuele Pesce](https://github.com/emanuelepesce) - [University of Salerno] (http://www.unisa.it)

For each suggestion or contribution don't hesitate to contact me.
