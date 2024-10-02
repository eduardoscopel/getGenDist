[![DOI](https://zenodo.org/badge/866575154.svg)](https://doi.org/10.5281/zenodo.13881896)


# getGenDist
This Python script computes the genetic distance between two strains using a distance matrix formatted for Phylip. It also performs k-means clustering to group strains based on their distances from a reference strain.

# Genetic Distance Calculator
This Python script computes the genetic distance between two strains using a distance matrix formatted for Phylip. It also performs k-means clustering to group strains based on their distances from a reference strain.

# Usage

`python getGenDist.py -i <infile> -a <strain1> -b <strain2> -r <ref> -k <clusters> -s <seed>`

```
Options
-i, --infile: Input genetic distance matrix (required) 
-a, --strain1: Name of the first strain (required)
-b, --strain2: Name of the second strain (required)
-r, --ref: Name of the reference strain (required)
-k, --clusters: Number of bins for clustering
-s, --seed: Random seed for sampling
-h, --help: Display help information
```

# Description

The script reads a genetic distance matrix, retrieves the distances between the specified strains, and utilizes Jenks natural breaks to cluster strains based on their genetic distance from the reference strain.

# Dependencies

pandas
jenkspy

Ensure you have the required libraries installed before running the script.



