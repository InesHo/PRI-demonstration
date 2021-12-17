# PRI-demonstration

This repo contains a simplified excerpt from PRI - Pattern Recognition of Immune Cells as we use it in our group (https://www.drfz.de/en/forschung/pb1/ag/signaltransduktion/) to analyze mass flow cytometry data. 
To ensure simple, comprehensible use and demonstration of the core functionality, we do not include a database connection and the GUI application for data analysis here. 



To run this code in R you additionally need to install the flowCore package using via BiocManager. This can be done b running the following code in your R-Terminal:
```
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
```

After you installed flowCore you can follow this Link to get to the demonstration: http://htmlpreview.github.io/?https://github.com/InesHo/Pri-demonstration/blob/master/code/demonstration_code.html


