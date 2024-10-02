[![DOI](https://zenodo.org/badge/866575154.svg)](https://doi.org/10.5281/zenodo.13881896)


# getGenDist
This Python script computes the genetic distance between taxa using a distance matrix obtained from Phylip. It then performs k-means clustering to group pairs of taxa based on their distances from a reference taxon (usually the root of a tree).

# Usage

`python getGenDist.py -i <infile> -a <taxon1> -b <taxon2> -r <ref> -k <clusters> -s <seed>`

```
Options
-i, --infile: Input genetic distance matrix (required) 
-a, --strain1: Name of the first taxon (required)
-b, --strain2: Name of the second taxon (required)
-r, --ref: Name of the reference taxon (required)
-k, --clusters: Number of bins for clustering
-s, --seed: Random seed for sampling
-h, --help: Display help information
```

# Description

The script reads a Phylip genetic distance matrix, retrieves the distances between the specified taxa, and utilizes Jenks natural breaks to cluster taxa based on their genetic distance from the reference taxon.

# Dependencies

pandas
jenkspy

Ensure you have the required libraries installed before running the script.



