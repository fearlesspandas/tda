This is some Python research code which works with SimPy (in Conda3)
to construct a small library for Topological Data Analysis. It's main purpose
is to illustrate the ease of describing simplicial sets using Generating Sets (Described here:
https://arxiv.org/pdf/1703.01547.pdf), and how this approach can be made algorithmically tractible.

It uses the language of Generating Sets to combinatorially generate the boundary matrix of a
simplicial complex. Besides this, the library also provides some standard TDA tools, such as
a function to calculate Betti numbers/homology groups. Further, there is a built in conversion
from incidence matrices to Generating Sets, as well as their associated boundary maps (algorithms
for this also described in the above paper). The incidence matrix may be a more natural first step
for constructing the Generating Set of a particular data set.

Complex.py is a standalone library that provides a class structure for simplicial complexes constructed from Generating Sets. tda.py and BettiCalc.py may share some of the same functionality but also have other helpful functions (Complex.py does not yet have conversion from incidence matrix... coming soon). 
