source('~/Dropbox/monkeybris/rscript/pedigree.preproc.batch.R')
source('~/code/snparray/SNPpreproc.R')
dat <- fread("~/analysis/SNPannotation/SNPmaster_qc.csv")
SNPdat <- fread("~/analysis/SNPannotation/SNPdat.csv")
#source("~/code/snparray/processraw.R")
reped <- ped.matchnames(dat$ID,pedigree$id)
pedigree <- ped.replace(pedigree,reped$oldID,reped$ID)
redped <- ped.trace(dat$ID,pedigree)
n <- nrow(redped)

numped <- data.frame(id=1:n,sire=vector("numeric",n),dam=vector("numeric",n))
numped$sire <- match(redped$sire,redped$id,nomatch = 0)
numped$dam <- match(redped$dam,redped$id,nomatch = 0)
numped <- as.matrix(numped)

#fv <- cullSNP(as.data.frame(dat[,-1,with=F]),SNPdat,mthresh=0.25,lthresh = 0.1)
X <- array(-1,dim = c(n,ncol(dat)-1))
X[match(dat$ID,redped$id),] <- as.matrix(dat[,-1,with=F])
X[is.na(X)] <- -1

if (!exists("type")) type = "none"

if (type == "eigenanalysis")
{
  shitlist <- SNPdelink(SNPdat)
  X <- X[,-shitlist]
  Xdat <- SNPdat[-shitlist]
  Ximp <- as.matrix(fread("~/analysis/SNPannotation/SNPmaster_qc_imputed.csv")[,-1,with=F])[,-shitlist]
} else if (type=="LD")
{
  m <- nrow(fv$SNP)
  fv$SNP$chrom <- as.character(fv$SNP$chrom)
  chrom_match <- matrix(fv$SNP$chrom,m,m) == matrix(fv$SNP$chrom,m,m,byrow = T)
  dist <- abs(matrix(fv$SNP$loc,m,m) - matrix(fv$SNP$loc,m,m,byrow=T))
  dist[!chrom_match] <- Inf
  shitlist <- which(colSums(dist<3e5)==1)
  X <- X[,-shitlist]
  dist <- dist[-shitlist,-shitlist]
}