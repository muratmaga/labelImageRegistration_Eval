library(ANTsR)
log = read.csv(file="/home2/rachel/P01/MagaLab/Data/Genotypes/embryos_2025-04-24_scanned.csv")
f=dir(path='/home2/rachel/P01/MagaLab/Data/RigidAligned/fullRes_simpleLabels/') #files with corrected labels
sample.name = substr(f, 1, 9)
log = log[which(log$ScanID %in% sample.name),c(2,6, 16)]

ko = grep("MUT", log$Genotype)
wt = c(4,  1, 26, 38, 17, 15, 25,  8, 28, 14)

log = log[c(wt, ko), ]
resample=c(0.054, 0.054, 0.054)

ref.image = antsImageRead('/home2/rachel/P01/MagaLab/AtlasBuilding/Atlas_02/Masked_NRRDs/maskedtemplate0.nii.gz')
ref.image = resampleImage(ref.image, resample)

ref.label = antsImageRead('/home2/rachel/P01/MagaLab/Data/RigidAligned/labatlas-labels-17_resampled-simple-label.nii.gz')
ref.label = resampleImage(ref.label, resample)

ref.mask  = antsImageRead('/home2/rachel/P01/MagaLab/Data/RigidAligned/maskedtemplate0__lowRes-wholebody-label.nii.gz')

n.labels = 16

n=nrow(log)
results=NULL



for (i in 1:n) {
  mov.image = antsImageRead(dir(patt=log$ScanID[i], path='/home2/rachel/P01/MagaLab/Data/RigidAligned/fullRes_Volumes/', full.names = TRUE) )
  mov.image = resampleImage(mov.image, resample)
  
  mov.label = antsImageRead(dir(patt=log$ScanID[i], path='/home2/rachel/P01/MagaLab/Data/RigidAligned/fullRes_simpleLabels/', full.names = TRUE) )
  mov.label = resampleImage(mov.label, resample)
  
  mov.mask = antsImageRead(dir(patt=log$ScanID[i], path='/home2/rachel/P01/MagaLab/Data/RigidAligned/lowRes_Masks-wholebody/', full.names = TRUE) )
  
  sim.reg = antsRegistration(fixed = ref.image, 
                             moving = mov.image, 
                             mask = ref.mask, 
                             movingMask = mov.mask, 
                             typeofTransform = 'Similarity')
  
  syn.reg = antsRegistration(fixed = ref.image, 
                             moving = mov.image, 
                             mask = ref.mask, 
                             movingMask = mov.mask, 
                             initialTransform = sim.reg$fwdtransforms,
                             typeofTransform = 'antsRegistrationSyNQuick[so]')
  
  ref.labels.in.subject.blind = antsApplyTransforms(fixed=mov.image, moving = ref.label, transformlist = syn.reg$invtransforms, interpolator = 'genericLabel')
  antsImageWrite(ref.labels.in.subject.blind, paste0("~/Desktop/validation_output/SyN.Labels/", log$ScanID[i], "-label.nii.gz"))
  
  
  # done with blind intensity registration. 
  
  # begin intensity registration leave one out 
  
  for (label in 1:16) {
    ref.label.temp = antsImageClone( ref.label)
    ref.label.temp[ref.label==label] = 0 #remove the label being held.
    mov.label.temp = antsImageClone( mov.label)
    mov.label.temp[mov.label==label] = 0 #remove the label being held.
    
    label.reg = labelImageRegistration(fixedLabelImages = ref.label.temp, 
                                       movingLabelImages = mov.label.temp, 
                                       fixedIntensityImages = ref.image, 
                                       movingIntensityImages = mov.image, 
                                       fixedMask = ref.mask, 
                                       movingMask = mov.mask, 
                                       initialTransforms = 'similarity', 
                                       typeOfDeformableTransform = 'antsRegistrationSyNQuick[so]')
    
    ref.labels.in.subject.label = antsApplyTransforms(fixed=mov.image, moving = ref.label, transformlist = label.reg$invtransforms, interpolator = 'genericLabel')
    antsImageWrite(ref.labels.in.subject.label, paste0("~/Desktop/validation_output/Label.Labels/", label,"/", log$ScanID[i], "-label.nii.gz"))
    
    # make new images with just the held label in them
    GS.label = antsImageClone(mov.label)
    GS.label[mov.label!=label] = 0
    GS.label[mov.label==label] = 1
    
    blind.label = antsImageClone(ref.labels.in.subject.blind)
    blind.label[ref.labels.in.subject.blind != label] = 0
    blind.label[ref.labels.in.subject.blind == label] = 1
    
    label.label = antsImageClone(ref.labels.in.subject.label)
    label.label[ref.labels.in.subject.label != label] = 0
    label.label[ref.labels.in.subject.label == label] = 1
    
    overlap.blind = labelOverlapMeasures(blind.label , GS.label)[1,]
    overlap.label = labelOverlapMeasures(label.label , GS.label)[1,]
    
    results=rbind(results, cbind(f[i], label, overlap.blind, overlap.label))
    
    remove(ref.label.temp, mov.label.temp, label.reg, ref.labels.in.subject.label, GS.label, blind.label, label.label, overlap.label, overlap.blind)
    
  }
}
write.csv(file="~/Desktop/validation_output/overlaps.csv", results)


