# PRI-demonstration

The purpose of this tutorial is to serve as a brief demonstration for ***PRI - Pattern Recognition of Immune Cells***. This is a novel approach in analysing and visualizing mass flow cytometry data developed in the research group [Baumgrass of the German Rheumatism Research Centre Berlin (DRFZ)](https://www.drfz.de/en/forschung/pb1/ag/signaltransduktion/).

Here we present a simplified excerpt from PRI to ensure comprehensible use and demonstration of the core functionality. 

PRI uses a bin-based apporach in plotting three marker combinations enabling the user to recognize activation patterns in the underlying data. We call these plots ***binplots***. In a binplot, the area of the markers X versus Y is categorized into bins of size 0.2 x 0.2. By binning, there are always several cells combined into one bin for which certain metrics can calculated. Using this approach, it is possible to visualize a third parameter in these plots. This third feature can depict a variety of options, for example the cell density or a statistical measure calculated for a third marker, referred to as marker Z. 

![plot](https://github.com/InesHo/PRI-demonstration/blob/main/images/general_example.png)



### Prerequisites
We use the library 'flowCore' v1.52.1 to read the FCS-FIles. If this library is not yet installed, it can be installed via [BioConductor](https://www.bioconductor.org/packages/release/bioc/html/flowCore.html) in the ***R terminal***. This must be done before the code examples from this tutorial can be used.
To do so run the following code in your ***R terminal***:

```
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
```

After you installed 'flowCore' you can either follow this [link to the demonstration](http://htmlpreview.github.io/?https://github.com/InesHo/Pri-demonstration/blob/master/code/demonstration_code.html) or you can source a simple R file from your ***systems terminal*** that recreates the plots from the demonstration code and saves them in the 'images' folder. 

Start by opening a terminal window and use the cd (change directory) commandline tool to navigate to the location where you want to save the files for this tutorial. Then clone the repository using 'git clone' . 

```
cd [path to destination for PRI-demonstration]
git clone https://github.com/InesHo/PRI-demonstration.git
```

Change into the 'PRI-demonstration' folder and execute the R-script.

```
cd PRI-demonstration
Rscript 'example/example.R'
```

This creates the result plots as they are presented in the [demonstration code](http://htmlpreview.github.io/?https://github.com/InesHo/Pri-demonstration/blob/master/code/demonstration_code.html). 



### Publications featuring PRI
- Gryzik, S., Hoang, Y., Lischke, T., Mohr, E., Venzke, M., Kadner, I., Poetzsch, J., Groth, D., Radbruch, A., Hutloff, A., & Baumgrass, R. (2020). Identification of a super-functional Tfh-like subpopulation in murine lupus by pattern perception. eLife, 9, e53226. https://doi.org/10.7554/eLife.53226
- tbc ...

