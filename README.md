# Rank Cover Tree (RCT) Library

> *Important* The document here is being migrated from the original inlined
> documentation in the code and is likely neither complete, nor in in sync 
> with the overall state of the refactoring I'm currently doing. Stay tuned.

![Master Build Status](https://travis-ci.org/mnett82/rct.svg?branch=master)

## Introduction

This package provides a C++ implementation of the rank cover tree (see [1]),
a probabilistic algorithm for finding k-nearest neighbors of point sets in
general metric spaces. The rank cover tree reinterprets the cover set
analysis of the cover tree search structure in terms of neighbor ranks as
measured from the query point, rather than explicit distances. The direct
analysis results in a construction and query time complexity with a far
smaller dependence on the measure of implicit dimension used in the cover
tree analysis. The rank cover tree can also be viewed as a variant of the
generic, scalable SASH heuristic for similarity search, where the rank
cover tree parameter governing the execution have been set to a value
insufficiently large to guarantee correctness with high probability. In this
sense, this implementation and its paper constitute the first formal
analysis of the accuracy of the SASH.

## Quick Start

The following examples assume that the array `DistData** data` contains (at 
least) `size` data items. In order to create an empty RCT use the following 
code.

```
RCT* rct = new RCT();
```

Alternatively, a specific seed value for the random number generator can be
passed as an argument to the constructor:

```
RCT* rct = new RCT(seed);
```

Following the creation of an RCT, the build parameters need to be specified.
More precisely, you may want to provide non-default values for

```
rct->setCoverageParameter(8.0);
rct->setSampleRate(pow(size, 1.0 / 3.0));
```

The RCT can now be initialized by calling 

```
rct->build(data, size);
```

After the construction finishes we can find approximate <em>k</em>
nearest-neighbors to a query object <code>DistData* query</code> by
calling:

```
int num = rct->findNear(query, k);
```

In the following we can access the indices and query-to-item distances
of the query results.

```
int* indices = new int [num];
rct->getResultIndices(indices, num);
float* distances = new float [num];
rct->getResultDists(distances, num);
```
