load("~/Dropbox/monkeybris/data/SNP.Rdat")
source('~/Dropbox/monkeybris/rscript/pedigree.preproc.batch.R')
source("~/Dropbox/tools/cheetahmm/cheetahmm.R")
source('~/Dropbox/monkeybris/rscript/setup.R')
reped <- ped.matchnames(as.character(fv$MON$Animal.ID),pedigree$id)
pedigree <- ped.replace(pedigree,reped$oldID,reped$ID)
redped <- ped.trace(fv$MON$Animal.ID,pedigree)
n <- nrow(redped)

numped <- data.frame(id=1:n,sire=vector("numeric",n),dam=vector("numeric",n))
numped$sire <- match(redped$sire,redped$id,nomatch = 0)
numped$dam <- match(redped$dam,redped$id,nomatch = 0)
numped <- as.matrix(numped)

X <- array(-1,dim = c(n,ncol(fv$X)))
X[match(rownames(fv$X),redped$id),] = as.matrix(fv$X)
X[is.na(X)] <- -1

shitlist <- c()
for (i in 1:length(unique(fv$SNP$chrom)))
{
  ic <- unique(fv$SNP$chrom)[i]
  loci <- which(fv$SNP$chrom == ic)
  loc <- fv$SNP$loc[loci]
  nl <- length(loc)

  D <- abs(kronecker(t(rep(1,times=nl)),loc) - kronecker(rep(1,times=nl),t(loc)))
  foo <- which(D<3e5 & D>0,arr=T)
  dump <- vector("numeric")
  
  moo <- foo[!(foo[,1] %in% dump | foo[,2] %in% dump),]
  while (length(moo) > 0)
  {
    dump <- c(dump, as.numeric(names(which.max(table(moo[,1])))))
    moo <- foo[!(foo[,1] %in% dump | foo[,2] %in% dump),]
  }
  
  shitlist <- c(shitlist,loci[dump])
}

X <- X[,-shitlist]
