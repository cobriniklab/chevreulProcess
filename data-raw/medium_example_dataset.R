# ## code to prepare `medium_example_dataset` goes here
#
# ## Not run:
# ## The code to prepare the .Rda file file from the larger dataset:
# medium_example_dataset <- chevreuldata::human_gene_transcript_sce()
# 
# altExp(medium_example_dataset, "transcript") <- NULL
# medium_example_dataset <- swapAltExp(medium_example_dataset, "gene") 
# altExp(medium_example_dataset, "integrated") <- NULL
# 
# myhvgs <- getTopHVGs(medium_example_dataset, n = 200)
# medium_example_dataset <- medium_example_dataset[myhvgs,]
# 
# assay(medium_example_dataset, "merged") <- NULL
# assay(medium_example_dataset, "reconstructed") <- NULL
# 
# colData(medium_example_dataset) <- colData(medium_example_dataset)[,duplicated(colnames(colData(medium_example_dataset)))]
#
# usethis::use_data(medium_example_dataset, overwrite = TRUE)
