# Reformat colData

Arbitrary colData can be appended to the data based on the results of exploratory data analysis. colData addition can be executed by uploading a csv with rownames matching cell ids with new variables as columns.

Clustering is done at different resolutions, ranging from 0.2-2. Additionally, two or more columns of the colData can be combined into one single column by simply selecting the columns of interest in SingleCellExperiment colData Object Columns and by providing a new name for the merged column.  
