# DEPENDENCIES ------------------------------------------

packages <- c("ggpubr",
              "reshape",
              "tidyverse",
              "tidyr",
              "vegetarian",
              "vegan")

is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}

for(package in packages){
  # check if package is installed
  if (!is.installed(package)){
    install.packages(package)
  }
  else {print(paste(package, " library already installed"))}
}

for(package in packages){
  eval(bquote(library(.(package))))
}


# # ARCHIVED DEPENDENCY ---------------------------------------------------

# Download package vegetarian from CRAN archive

url <- "https://cran.r-project.org/src/contrib/Archive/vegetarian/vegetarian_1.2.tar.gz"
pkgFile <- "vegetarian_1.2.tar.gz"
download.file(url = url, destfile = pkgFile)

# Expand the zip file using whatever system functions are preferred
# look at the DESCRIPTION file in the expanded package directory
# Install dependencies list in the DESCRIPTION file

# install.packages(c("ada", "ipred", "evd"))

# Install package
install.packages(pkgs=pkgFile, type="source", repos=NULL)

# Delete package tarball
unlink(pkgFile)



# RANGE - ALPHA DIV (t-tests) ---------------------------------------------

#FILE SETUP
dat <- read.csv("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/YST1-metadata-allSoil.csv", row.names = 1, header = TRUE)
tab <- read.csv("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/feature-table.csv", row.names = 1, header = TRUE)
# check row names align, should be 0
# sum(as.numeric(rownames(dat)) - as.numeric(rownames(tab)))

#CALCULATE RICHNESS
q0a <- apply(tab,1,d,q=0, lev="alpha")
# add to dataframe
dat$q0 <- q0a

# TTESTS 
# Q0 native alpha diversity vs invaded alpha diversity
t.test(q0a[which(dat$range=="Native")],
       q0a[which(dat$range=="Invaded")], 
       paired = FALSE, 
       alternative = "two.sided", 
       var.equal = FALSE)
# shannon native alpha diversity vs invade alpha diversity
t.test(dat$shannon_r1000[which(dat$range=="Native")],
       dat$shannon_r1000[which(dat$range=="Invaded")], 
       paired = FALSE, 
       alternative = "two.sided", 
       var.equal = FALSE)




# EXPT - ALPHA DIV (linear models) ----------------------------------------

#FILE SETUP
dat <- read.csv("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/YST1-metadata-rootwash.csv", row.names = 1, header = TRUE)
dat <- dat[(which(dat$feature.count>0)),]

tab <- read.csv("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/feature-table-rootwash-t.csv", row.names = 1, header = TRUE)

# check files in same order
cbind(rownames(dat),rownames(tab))

#CALCULATE RICHNESS
q0a <- apply(tab,1,d,q=0, lev="alpha")
# add to dataframe
dat$q0 <- q0a

# remove control samples
tab <- tab[(which(dat$soilRange!="Control")),]
dat <- dat[(which(dat$soilRange!="Control")),]


# Model with interactions
q.aov2 <- aov(q0 ~ seedRange + soilRange + JulGerm + StatsBlock + seedRange*soilRange, data = dat)
summary(q.aov2)
s.aov2 <- aov(shannon_r1000 ~ seedRange + soilRange + JulGerm + StatsBlock + seedRange*soilRange, data = dat)
summary(s.aov2)

TukeyHSD(x=q.aov2, 'seedRange', conf.level=0.95)
TukeyHSD(x=s.aov2, 'seedRange', conf.level=0.95)

# BETA DIV (permanova/adonis) ---------------------------------------------

dat <- read.csv("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/YST1-metadata-allSoil.csv", row.names = 1, header = TRUE)

# remove samples with <1000 seqs
dat <- subset(dat, feature.count > 1000)
tab<-tab[(which(rownames(dat) %in% rownames(tab))),]

# Distance matricies
uwu_dm <- read.table("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/diversity-fieldSoil/uwu-distance-matrix.tsv", header = TRUE, row.names = 1)
uwu_dm <- as.dist(uwu_dm)
wu_dm <- read.table("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/diversity-fieldSoil/wu-distance-matrix.tsv", header = TRUE, row.names = 1)
wu_dm <- as.dist(wu_dm)

# ADONIS - nested
set.seed(1990)
adonis(formula = uwu_dm ~ range + range/Soil.Site,
       data = dat,
       permutations = 9999,
       by = margin) #type3 sum of sq

set.seed(1990)
adonis(formula = wu_dm ~ range + range/Soil.Site,
       data = dat,
       permutations = 9999,
       by = margin) #type3 sum of sq







# INOC COMPARISONS --------------------------------------------------------
dat <- read.csv("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/YST1-metadata-allSoil.csv", row.names = 1, header = TRUE)

# remove samples with <1000 seqs
dat <- subset(dat, feature.count > 1000)

uwU = read.table("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/diversity-fieldSoil/uwu-distance-matrix.tsv", row.names = 1)
wU = read.table("~/Dropbox/projects/Starthistle1/dataAnalysis_YST1/diversity-fieldSoil/wu-distance-matrix.tsv", row.names = 1)

# specifying which site combos for functions (all with sequenced innoculum bulk)
bulk_sites <- c("CAZ_YST", "GIL_nonYST", "GIL_YST", "GRA_nonYST", "GRA_YST","LEB_nonYST", "LEB_YST", "TRI_nonYST")


#dataframe where all the results will go
ttestResults <- (matrix(c("group1", "group2", "distance", "t", "df", "p"), nrow = 1))


# function to pull correct pairwise distances - WITHIN SITE VS BETWEEN SITE
# (average of innoc sample compared to all that site's other samples,
# vs innoc sample compared to average of all other sites' (non innoc) samples)


prws_dsts <- function(site = NA, mat = NA, dat = NA, reftyp = NA, alttyp = NA){
  # site = location and whether it's in YST patch or not
  # mat = distance matrix to use
  # dat = metadata file 
  # reftyp & alttyp = whether it's a sample that was used for innoc or not 
  #                   (yes/no as coded in metadata)
  
  # pull the parts of the matrix we want
  # focal site reference sample type
  focalsamp <- which(dat$sourceCode==site & dat$used.as.inoc==reftyp)
  # focal site alternate sample type
  altsamesite <- which(dat$used.as.inoc==alttyp & dat$sourceCode==site)
  # alternate sample type, all other sites
  othersamp <- which(dat$sourceCode!=site)
  
  # pulling all rows/columns from dm for among and averaging, excluding focal site's alttyp
  # (distance of among comparison for focal reftyp vs all other site's alttyp)
  # averaging first row without focal lgint
  self <- mean(as.matrix(mat[c(focalsamp,altsamesite),c(focalsamp,altsamesite)])[1,-1])
  other <- mean(as.matrix(mat[c(focalsamp,othersamp),c(focalsamp,othersamp)])[1,-1])
  
  return(c(self,other))
}

#Within vs between - paired t-tests - unweighted uniFrac

# make array of all within (column 1) and all between (column 2) using function
inoc_compar_uw <- array(dim=c(length(bulk_sites),2))

for(i in (1:(length(bulk_sites)))){
  inoc_compar_uw[i,] <- prws_dsts(bulk_sites[i], mat = uwU, dat = dat, reftyp = "yes", alttyp = "no")
}
# paired t-test of column 1 and column 2 of inoc_compar_uw
tRes <- t.test(inoc_compar_uw[,1],inoc_compar_uw[,2],paired=TRUE)


ttestResults <- rbind(ttestResults,
                      c("inoc-other within site", 
                        "inoc-other outside site", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))



#Within vs between - paired t-tests - weighted uniFrac

# make array of all within (column 1) and all between (column 2) using function
inoc_compar_w <- array(dim=c(length(bulk_sites),2))

for(i in (1:(length(bulk_sites)))){
  inoc_compar_w[i,] <- prws_dsts(bulk_sites[i], mat = wU, dat = dat, reftyp = "yes", alttyp = "no")
}
# paired t-test of column 1 and column 2 of large_f_uw
tRes <- t.test(inoc_compar_w[,1],inoc_compar_w[,2],paired=TRUE)


ttestResults <- rbind(ttestResults,
                      c("inoc-other within site", 
                        "inoc-other outside site", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))



# INOC COMPARISONS (plotting) ---------------------------------------------



# same as above just adding site names as a column
prws_dsts <- function(site = NA, mat = NA, dat = NA, reftyp = NA, alttyp = NA){
  # site = location and whether it's in YST patch or not
  # mat = distance matrix to use
  # dat = metadata file 
  # reftyp & alttyp = whether it's a sample that was used for innoc or not 
  #                   (yes/no as coded in metadata)
  
  # pull the parts of the matrix we want
  # focal site reference sample type
  focalsamp <- which(dat$sourceCode==site & dat$used.as.inoc==reftyp)
  # focal site alternate sample type
  altsamesite <- which(dat$used.as.inoc==alttyp & dat$sourceCode==site)
  # alternate sample type, all other sites
  othersamp <- which(dat$sourceCode!=site)
  
  # pulling all rows/columns from dm for among and averaging, excluding focal site's alttyp
  # (distance of among comparison for focal reftyp vs all other site's alttyp)
  # averaging first row without focal lgint
  self <- mean(as.matrix(mat[c(focalsamp,altsamesite),c(focalsamp,altsamesite)])[1,-1])
  other <- mean(as.matrix(mat[c(focalsamp,othersamp),c(focalsamp,othersamp)])[1,-1])
  
  return(c(site,self,other))
}

# make array of all within (column 1) and all outside (column 2) using function
inoc_compar_uw <- array(dim=c(length(bulk_sites), 3))
colnames(inoc_compar_uw) <- c("site", "within", "outside")
for(i in (1:(length(bulk_sites)))){
  inoc_compar_uw[i,] <- prws_dsts(bulk_sites[i], mat = uwU, dat = dat, reftyp = "yes", alttyp = "no")
}

# changing from array to dataframe so I can melt
inoc_compar_uw <- as.data.frame(inoc_compar_uw)
# melting
mdat <- melt(inoc_compar_uw, id="site")
# separating first column for better plotting
mdat <- separate(data = mdat, col = site, into = c("site", "patch"), sep = "_")
mdat$value <- as.numeric(mdat$value)


uwU <- ggplot(mdat, 
              aes(x = variable, y = value)) + 
  geom_boxplot() +  
  xlab("") +
  ylab("unweighted unifrac distance") +
  # ylim(0,1.25) +
  theme_classic() 


# weighted unifrac same as above
inoc_compar_w <- as.data.frame(inoc_compar_uw)
mdat <- melt(inoc_compar_w, id="site")
mdat <- separate(data = mdat, col = site, into = c("site", "patch"), sep = "_")
mdat$value <- as.numeric(mdat$value)

wU <- ggplot(mdat, 
             aes(x = variable, y = value)) + 
  geom_boxplot() +  
  xlab("") +
  ylab("weighted unifrac distance") +
  # ylim(0,1.25) +
  theme_classic() 

ggarrange(wU, uwU,
          labels = c("a)","b)"))
