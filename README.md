The Glasgow Subgraph Solver for Solving Crosswords
==================================================

This is a solver for subgraph isomorphism (induced and non-induced) problems, based upon a series of
papers by subsets of Blair Archibald, Ciaran McCreesh, Patrick Prosser and James Trimble at the
University of Glasgow, and Fraser Dunlop and Ruth Hoffmann at the University of St Andrews. A clique
decision / maximum clique solver is also included.

If you use this software for research, please cite [icgt/McCreeshPT20]. If you use this solver in a
non-research setting, please get in touch if you can. This software is an output of taxpayer funded
research, and it is very helpful for us if we can demonstrate real-world impact when we write grant
applications.

Please contact [Ciaran McCreesh](mailto:ciaran.mccreesh@glasgow.ac.uk) with any queries.

Compiling
---------

To build, type 'make'. You will need a C++17 compiler (we test with GCC 7.3, GCC 8.3, and Clang
7.0.1 on Linux, and Xcode 10.2 on Mac OS X) and Boost (we use 1.65.1 or later, built with threads
enabled).

Running
-------

To run:

```shell session
$ ./glasgow_subgraph_solver --crossword pattern-file target-file
```

The default mode is to display the first found solution, or to prove unsatisfiability if no solution
exists. To count or print all solutions, use one of:

```shell session
$ ./glasgow_subgraph_solver [ --count-solutions | --print-all-solutions ] pattern-file target-file
```

Note that printing all solutions can be exponentially slower than counting solutions.

The solver supports parallel search. Usually you should enable this, as follows:

```shell session
$ ./glasgow_subgraph_solver --parallel ...
```

Note that parallel search, in its default configuration, is non-deterministic.

File Formats
------------

We expect the files to csvs.

In particular, note that auto-detection can easily fail if, for example, the first vertex in the
graph has no neighbours.  We can read LAD, Labelled LAD (labels on vertices, and optionally also on
edges), CSV, and DIMACS 2 formatted graphs. [The LAD
format](https://perso.liris.cnrs.fr/christine.solnon/SIP.html) is a nice simple choice. If you need
to support named vertices, labels on vertices and / or edges, or directed edges, consider using the
CSV format. To specify a directed edge, use a greater-than sign rather than a comma as the delimiter
between the first two columns.  To specify an edge label, include a third column in the file. To
specify a vertex label, leave the second column empty and use the third column for the label. For
example, the following describes a graph with four vertices, with colours for edge labels and shapes
for vertex labels.

```
first>second,red
second>first,blue
first,third,purple
first>first,green
first,,circle
second,,circle
third,,square
fourth,,square
```


Funding Acknowledgements
------------------------

This work was supported by the Engineering and Physical Sciences Research Council (grant numbers
EP/P026842/1, EP/M508056/1, and EP/N007565). This work used the Cirrus UK National Tier-2 HPC
Service at EPCC (http://www.cirrus.ac.uk) funded by the University of Edinburgh and EPSRC
(EP/P020267/1).

References
----------

* [walcom/TomitaSHTW10]: https://dblp.org/rec/html/conf/walcom/TomitaSHTW10
  **walcom/TomitaSHTW10**:
  Etsuji Tomita, Yoichi Sutani, Takanori Higashi, Shinya Takahashi, Mitsuo Wakatsuki:
  A Simple and Faster Branch-and-Bound Algorithm for Finding a Maximum Clique. WALCOM 2010: 191-203.
  DBLP: [walcom/TomitaSHTW10]

* [algorithms/Prosser12]: https://dblp.org/rec/html/journals/algorithms/Prosser12
  **algorithms/Prosser12**:
  Patrick Prosser: Exact Algorithms for Maximum Clique: A Computational Study. Algorithms 5(4):
  545-587 (2012). DBLP: [algorithms/Prosser12].

* [ol/SegundoMRH13]: https://dblp.org/rec/html/journals/ol/SegundoMRH13
  **ol/SegundoMRH13**:
  Pablo San Segundo, Fernando Matía, Diego Rodríguez-Losada, Miguel Hernando: An improved bit
  parallel exact maximum clique algorithm. Optimization Letters 7(3): 467-479 (2013). DBLP:
  [ol/SegundoMRH13].

* [cp/McCreeshP14]: https://dblp.org/rec/html/conf/cp/McCreeshP14
  **cp/McCreeshP14**:
  Ciaran McCreesh, Patrick Prosser: Reducing the Branching in a Branch and Bound Algorithm for the
  Maximum Clique Problem. CP 2014: 549-563. DBLP: [cp/McCreeshP14].

* [jco/BatsynGMP14]: https://dblp.org/rec/html/journals/jco/BatsynGMP14
  **jco/BatsynGMP14**:
  Improvements to MCS algorithm for the maximum clique problem. J. Comb. Optim. 27(2): 397-416
  (2014). DBLP: [jco/BatsynGMP14]

* [cp/McCreeshP15]: https://dblp.org/rec/html/conf/cp/McCreeshP15
  **cp/McCreeshP15**:
  Ciaran McCreesh, Patrick Prosser: A Parallel, Backjumping Subgraph Isomorphism Algorithm Using
  Supplemental Graphs. CP 2015: 295-312. DBLP: [cp/McCreeshP15].

* [lion/KotthoffMS16]: https://dblp.org/rec/html/conf/lion/KotthoffMS16
  **lion/KotthoffMS16**:
  Lars Kotthoff, Ciaran McCreesh, Christine Solnon: Portfolios of Subgraph Isomorphism Algorithms.
  LION 2016: 107-122. DBLP: [lion/KotthoffMS16].

* [aaai/HoffmannMR17]: https://dblp.org/rec/html/conf/aaai/HoffmannMR17
  **aaai/HoffmannMR17**:
  Ruth Hoffmann, Ciaran McCreesh, Craig Reilly: Between Subgraph Isomorphism and Maximum Common
  Subgraph. AAAI 2017: 3907-3914. DBLP: [aaai/HoffmannMR17].

* [jair/McCreeshPST18]: https://dblp.org/rec/html/journals/jair/McCreeshPST18
  **jair/McCreeshPST18**:
  Ciaran McCreesh, Patrick Prosser, Christine Solnon, James Trimble: When Subgraph Isomorphism is
  Really Hard, and Why This Matters for Graph Databases. J. Artif. Intell. Res. 61: 723-759 (2018).
  DBLP: [jair/McCreeshPST18].

* [cpaior/ArchibaldDHMPT19]: http://dblp.org/rec/html/conf/cpaior/ArchibaldDHMP019
  **cpaior/ArchibaldDHMPT19**:
  Blair Archibald, Fraser Dunlop, Ruth Hoffmann, Ciaran McCreesh, Patrick Prosser and James Trimble:
  Sequential and Parallel Solution-Biased Search for Subgraph Algorithms. CPAIOR 2019: 20-38.
  DBLP: [cpaior/ArchibaldDHMPT19].

* [icgt/McCreeshPT20]: http://dblp.org/rec/html/conf/gg/McCreeshP020
  **icgt/McCreeshP020**:
  Ciaran McCreesh, Patrick Prosser, James Trimble:
  The Glasgow Subgraph Solver: Using Constraint Programming to Tackle Hard Subgraph Isomorphism
  Problem Variants. ICGT 2020: 316-324.
  DBLP: [icgt/McCreeshPT20].

<!-- vim: set tw=100 spell spelllang=en : -->
