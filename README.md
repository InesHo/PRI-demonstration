# PRI-demonstration

The purpose of this repo is to serve as a breif demonstration and tutorial for PRI - Pattern Recognition of Immune Cells. This is a novel approach in analyszing and visualizing mass flow cytometry data developed in our group (https://www.drfz.de/en/forschung/pb1/ag/signaltransduktion/).

Here we present a simplified excerpt from PRI to ensure simple, comprehensible use and demonstration of the core functionality. We do not include a database connection and the GUI application for data analysis here. 

PRI uses a bin-based apporach in plotting three marker combinations enabling us to recognize activation patterns in the underlying data. We call these plots 'binplots'. In a binplot, the area of the markers X versus Y is categorized into bins of size 0.2 x 0.2. By binning, we always combine several cells into one bin for which we can then calculate certain metrics. Using this approach, we can visualize a third parameter in these plots. This third feature can depict a variety of options, for example the cell density or a statistical measure calculated for a third marker, referred to as marker Z. 

After you installed flowCore (see below) you can follow this Link to get to the demonstration: http://htmlpreview.github.io/?https://github.com/InesHo/Pri-demonstration/blob/master/code/demonstration_code.html

To run this code in R you additionally need to install the flowCore package using via BiocManager. This can be done b running the following code in your R-Terminal:
```
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
```

