# CellQuantAndLoc
CellQuantAndLoc quantifies yeast cell protein locations at the cell membrane and vacuolar level.  

## Pipeline Overview

1.  **Preprocessing**:  conversion to grayscale and normalization of all channels.

![Preprocessing](/images/1-Orig-Grayscale-Normalized.png)


2.  **Membrane localizations**.  Image processing is applied to the GFP channel to detect PM localizations.  Detected segments are subjected to exclusionary criteria based on location, circularity, size, etc.  Excluded segments are removed.


![PM localizations](/images/2-Membranes-Removed-Inner-Membrane.png)


3.  **Vacuolar localizations**.  A simple threshold is applied to the CMAC channel to get a first pass on vacuolar localizations.  As the image is a 2D slice all detected localizations may not be contiguous.  Further exclusionary filters are applied.

4.  **Computed Features**. Each cell with its associated PM localization, vacuolar localizations, and computed features are sorted into a dataframe.  Further exclusionary criteria applies.

5.  **Output**.  Final output is exported to a .csv file.

6.  **Visualization**.  Pipeline can be run entirely in the CellQuantAndLoc Visualization Tool.  Results can also be saved as a .RDS file and read into the visualization tool for quick processing.



![Summarized Result](/images/output1.png)

![Table](/images/table.png)

![Output log](/images/output-log.png)


### Requirements
1.  **Input Format**: CMAC, GFP, and DIC .tif files. Working on .nd2 compatibility.

2.  **R 3.6.2** or greater

3.  **Packages**: Bioconductor::EBImage, R-Shiny
