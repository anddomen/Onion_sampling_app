# Product Sampling App

**A link to the app is [here.](https://ahrouj-andrea0domen.shinyapps.io/product_sampling_app/)**

A little Shiny app used to simulate finished product sampling to detect contaminated lots. The app uses Monte Carlo simulations to evaluate the effectiveness of different sample sizes taken from a lot of a given size and contamination prevalence. It can compare two different scenarios side-by-side. The results include the probability of detecting positive lots and the distribution of positive samples within those lots. This was created for use as an education tool.

## Features

The following input parameters are user-adjustable:

-   Lot size (1,000 to 30,000,000)

-   Number of samples taken (1 to 500 samples)

-   Number of contaminated units

Visual results:

-   Stacked bar chart showing detection rate by lot

-   Distribution histogram of positive samples within detected lots

## How it Works

The simulations rely on a Zero-Adjusted Gamma (ZAGA) distribution to simulate contamination levels across samples. The prevalence parameter $\nu$ is calculated based on user input. The location and scale parameters ($\mu$ and $\sigma$ respectively) were determined from combined experimental data. Rather than using a binomial distribution, the ZAGA distribution allows us to examine the severity of the contamination. Currently the app does not utilize this, but we hope to in future versions. The user-specified number of samples are drawn from this distribution and aggregated. This is repeated for 10,000 simulated lots.

## Acknowledgements

Special thanks to Dr. James Molyneux for coding help and optimization. Claude was used in the creation of this app to address errors, however, the bulk of the code remains human written.
