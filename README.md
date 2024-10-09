
# Population Diversity Analysis Tool

This Python script calculates the diversity of multiple populations (based on Heterozygosity and SSD) and generates correlation plots for allele frequency data. It includes functions for calculating population-level diversity measures: pairwise differencing, pooling, averaging, and fixing.

## Features

- **Reading Data:** Load population frequency data from a CSV file. Each row corresponds to a population and each column is the allele frequency of a single locus.
- **Correlation Plots:** Generate scatter plots showing the correlation of all pairs of measures for population sets of specified sizes k.
- **Diversity Measures:**
  - Het (Heterozygosity) and SSD (Split System Diversity) calculated through various methods:
    - Pooling: all populations of a given set of size k are pooled together and treated as a single population.
    - Averaging: the diversity is calculated within each population in the set and then the average diversity is calculated.
    - Pairwise Differencing: the average pairwise differences of allele frequencies are calculated for all pairs of populations in the set.
    - Fixing: pooled diversity is calculated after fixation.
- **Brute Force Analysis:** Calculate diversity measures for all possible combinations of population sets of a specified size k.

## Requirements

- Python 3.x
- Libraries:
  - numpy
  - pandas
  - matplotlib

## Usage

1. Clone the repository or download the script file.
2. Make sure you have Python and the required libraries installed.
3. Run the script with appropriate parameters.

## How to Use

### 1. Reading Data

Ensure your population frequency data is in CSV format and located in the appropriate directory. By default, the script looks for a file named `test_example.csv`. You can specify a different file path if needed. 

```python
# Example usage
d = load_frequencies()
```

### 2. Correlation Plots

Generate scatter plots showing the correlation of all pairs of measures for population sets of specified sizes k. The script looks for a file named 'test_example__k2_output.csv' which is the output of the brute-force() of all measures. Adjust the parameters as needed.

```python
# Example usage
scatter_plots_all_pairs(k=2, path_df='test_example__k2_output.csv', col_index='Unnamed: 0')
```

### 3. Diversity Measures

Use the provided functions to calculate diversity measures. Adjust parameters as needed.

```python
# Example usage
het_score = Het_Pooling(freqs)
ssd_score = SSD_Pooling(freqs)
```

### 4. Brute Force Analysis

Perform a brute force analysis to calculate diversity measures for all possible combinations of population sets of a specified size k. Results are saved to a CSV file.

```python
# Example usage
brute_force_all_HETandSSD(pop_freqs=d, k=2, save_file_as='output.csv')
```

## Data Source

The original VCF files used in this project are from the paper:

- Moore et al., "Conservation genomics of anadromous Atlantic salmon across its North American range: outlier loci identify the same patterns of population structure as neutral loci
," *Molecular Ecology*, 2014. [Link to Paper]( https://doi.org/10.1111/mec.12972)

The analyzed data includes the allele frequency of multiple loci in 50 Atlantic salmon populations.

## License

This project is licensed under the Apache License 2.0. Users of this code must cite both the repository and the paper as follows:

- Repository: [population-diversity-tool](https://github.com/nabhari/population-diversity-tool)
- Paper: Abhari et al., Measuring  Genetic Diversity Across Populations, 2024 [under review].

See the [LICENSE](http://www.apache.org/licenses/LICENSE-2.0) file for more details.

