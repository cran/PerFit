# Define print() function for class "PerFit":
print.PerFit <- function(x, #x = an object from 'PerFit' class
                         cutoff.obj=NULL, #cutoff.obj = an object from 'PerFit.cutoff' class
                         ModelFit="NonParametric", Nreps=1000, 
                         IP=x$IP, IRT.PModel=x$IRT.PModel, Ability=x$Ability, Ability.PModel=x$Ability.PModel, mu=0, sigma=1, 
                         Blvl = 0.05, Breps = 1000, CIlvl = 0.95, 
                         UDlvl = NA, ...)
{
  N <- dim(x$Matrix)[1]; I <- dim(x$Matrix)[2]
  # Sanity check - Class PerFit:
  Sanity.cls(x)  
  # 
  dico.PFS <- c("Cstar", "C.Sato", "U3", "ZU3", "G", "Gnormed", "D.KB", "r.pbis", "NCI", "Ht", "A.KB", "E.KB", "lz", "lzstar")
  poly.PFS <- c("Gpoly", "Gnormed.poly", "U3poly", "lzpoly")
  
  # Compute cutoff:
  if (is.null(cutoff.obj))
  {
    cutoff.res <- cutoff(x, ModelFit, Nreps, IP, IRT.PModel, Ability, Ability.PModel, mu, sigma, Blvl, Breps, CIlvl, UDlvl)
  } else
  {
    Sanity.clsPO(cutoff.obj)
    cutoff.res <- cutoff.obj
  }
  
  # Compute flagged:
  flagged.res <- flagged.resp(x, cutoff.res, scores=FALSE)[[1]][,1]
  
  # Summarize results:
  flagged.bin <- rep("", N)
  flagged.bin[flagged.res] <- "*"
  all.PFS <- data.frame(PerFit.SE(x), Flagged=flagged.bin)
  print(all.PFS)
  # 
  cat(paste0("\nPFS = ", x$PFStatistic, "\n"))
  cat(paste0("Cutoff = ", cutoff.res$Cutoff, " (SE = ", cutoff.res$Cutoff.SE, ").\n"))
  cat(paste0("Tail = ", cutoff.res$Tail, ".\n"))
  cat(paste0("Proportion of flagged respondents = ", cutoff.res$Prop.flagged, ".\n"))
  cat("(N.B.: The cutoff varies each time cutoff() is run due to bootstrapping.)\n\n")
  # 
  cat(paste0("Identified respondents - ", length(flagged.res), " in total:\n"))
  cat("   ", flagged.res, "\n\n")
}
