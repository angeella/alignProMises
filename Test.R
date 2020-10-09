
require(ARIpermutation)
require(vMFPmodel)

mask_STG <- RNifti::readNifti("C:/Users/Angela Andreella/Documents/Thesis_Doc/Hyperaligment/Computation/AuditoryData/Pvalues/STG/mask_Superior_Temporal_Gyrus.nii.gz")
sub_ids <- sapply(c(21:24),function(x) paste0(0,x))
img <- array(NA, c(sum(mask_STG==1), 310, length(sub_ids)))

for (sid in 1:length(sub_ids)) {  
  wb <- RNifti::readNifti(paste0("C:/Users/Angela Andreella/Documents/Thesis_Doc/Hyperaligment/Computation/AuditoryData/Data_preprocess_FSL/sub-", sub_ids[sid],".feat.nii.gz"))
  voxel <- dim(wb)[1] * dim(wb)[2] * dim(wb)[3]
  img[,,sid] <- array(wb, dim = c(voxel, 310))[which(mask_STG==1),]
}

data <- img
data<- array(rnorm(23*23*5), dim = c(23,23,5))

system.time(out <-vMFPmodel(data, maxIt = 1, t = 1, k = 1, scaling = TRUE, reflection = TRUE, subj = FALSE, centered = FALSE))
