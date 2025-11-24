# cluster_metrics: Calculates clustering metric scores

Uses the ClusterSim package to calculate metrics for cluster
consistency. These metrics do not require a priori knowledge of the
truth set, which makes them ideal for scRNA clustering.

## Usage

``` r
cluster_metrics(
  so,
  cluster_list,
  dims = 1:20,
  reduction = "pca",
  silhouette = TRUE
)
```

## Arguments

- so:

  A Seurat single cell RNA (scRNA) object

- cluster_list:

  A vector of cluster resolution names (e.g. "SCT_snn_res.0.2")

- dims:

  An integer vector of the principal components to be used in
  calculating the cluster metrics. May not exceed the number of
  principal components available

- reduction:

  A string for the reduction to be used in calculating scores

- silhouette:

  A Boolean for whether to include silhouette scoring

## Value

Returns a table of clustering metrics for each cluster resolution
selected

## Details

Calculates three consistency scores to evaluate clustering
effectiveness: Calinski-Harabasz: Lower scores preferred.
Davies-Bouldin: Higher scores preferred. Silhouette score: Higher scores
preferred. Recommended to only run when data has fewer than 90,000 data
points
