# multiDEG
R package to perform DEG analysis using Wilcox rank test, DEGseq2, and edgeR methods

# Installing

```
devtools::install_github("nicolau/multiDEG")
```


# Loading dependencies packages

```
library(DESeq2)
library(edgeR)
library(tidyverse)
```

# Loading example expression and phenodata files

```
data(raw.exp)
data(phenodata)
```

# Get Venn Diagram plot for overlap between three methods
## For Up-regulated genes
```
plot(listDEGs = results, type = "up")
```

## For Down-regulated genes
```
plot(listDEGs = results, type = "down")
```
