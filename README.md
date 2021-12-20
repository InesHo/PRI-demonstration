# PRI-demonstration

The purpose of this tutorial is to serve as a brief demonstration for ***PRI - Pattern Recognition of Immune Cells***. This is a novel approach in analyszing and visualizing mass flow cytometry data developed in the research group [Baumgrass of the German Rheumatism Research Centre Berlin (DRFZ)](https://www.drfz.de/en/forschung/pb1/ag/signaltransduktion/).

Here we present a simplified excerpt from PRI to ensure comprehensible use and demonstration of the core functionality. We do not include database connection and GUI application for data analysis here. 

PRI uses a bin-based apporach in plotting three marker combinations enabling us to recognize activation patterns in the underlying data. We call these plots ***binplots***. In a binplot, the area of the markers X versus Y is categorized into bins of size 0.2 x 0.2. By binning, we always combine several cells into one bin for which we can then calculate certain metrics. Using this approach, we can visualize a third parameter in these plots. This third feature can depict a variety of options, for example the cell density or a statistical measure calculated for a third marker, referred to as marker Z. 

![plot](https://github.com/InesHo/PRI-demonstration/blob/main/images/general_example.png)

### Prerequisites
We use the library 'flowCore' v1.52.1 to read the FCS-FIles. If this library is not yet installed, it can be done via [BioConductor](https://www.bioconductor.org/packages/release/bioc/html/flowCore.html) in the ***R terminal***. This must be done before you can use the code examples from this tutorial. 
To do so run the following code in your ***R terminal***:

```
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
```

After you installed 'flowCore' you can either follow this link to get to the [demonstration](http://htmlpreview.github.io/?https://github.com/InesHo/Pri-demonstration/blob/master/code/demonstration_code.html) or you can source a common R file from your ***systems terminal*** that recreates the plots from the demonstration code and saves them in the 'images' folder. To do so open a terminal window and use cd (change directory) command line tool to change tothe 'PRI-demonstration' directory.

```
cd [path to PRI-demonstration]
```

From this folder, please execute the command
```
Rscript 'example/example.R'
```

to create the result plots as there are presented in the [demonstration code](http://htmlpreview.github.io/?https://github.com/InesHo/Pri-demonstration/blob/master/code/demonstration_code.html). 
