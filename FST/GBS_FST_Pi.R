### R modules for calculating pairwise FST and nucleotide diversity in a population with 2 subpopulations. 
### Specifically for GBS data with a high percentage of missing data. 
### Requires .geno and .ind files in Eigensoft format. Includes ability to set a minimum percentage of sample coverage.

## Function to read in a .geno file. Returns a list with a matrix of genotypes, 
## a vector with population assignments for each individual (0 or 1), and
## a vector with randomized population assignments for sampling null.
makeGenoMatrix <- function(genofile, indfile){
    geno_data = scan(genofile,what = "character")
    L = length(geno_data)
    S = nchar(geno_data[[1]])
    geno_matrix = matrix(data = NA, nrow = S, ncol = L)
    for (i in 1:L){
        genostring = strsplit(geno_data[[i]], split="")
        for (j in 1:S){
            geno_matrix[j,i] = as.numeric(genostring[[1]][j])}
    }
    geno_matrix[geno_matrix == 9] = NA
    ##Read in .ind file in order to get population assignment for genotypes.
    testpop = read.table(indfile,header=FALSE)
    test_vector = testpop$V2
    test_sub = {}
    for (i in 1:S){
        if (test_vector[i] == 1){
            test_sub = c(test_sub,i)}
    }
    ## Creates a null vector, where individuals are randomly assigned to the 2 populations.
    null_vector = sample(c(rep(1,sum(test_vector)),rep(0,S-sum(test_vector))),size=S)
    null_sub = {}
    for (i in 1:S){
        if (null_vector[i] == 1){
            null_sub = c(null_sub,i)}
    }
    return(list(geno_matrix,test_sub,null_sub))
}

##THIS MODULE CALCULATES PI (NUCLEOTIDE DIVERSITY) FOR EACH LOCUS AND WRITES TO A FILE
## Note that this is pi when considering individuals from just 2 populations. Will need to adapt to get pi from 
## entire dataset.
getPi <- function(geno_matrix, suf,mincov){
    pi = {}
    L <- ncol(geno_matrix)
    S <- nrow(geno_matrix)
    for (i in 1:L){
        locus = {}
        locus = geno_matrix[,i]
        locus = locus[!is.na(locus)]
        if(length(locus) >= S*mincov){
            locus_sum = 0
            for (j in 2:length(locus)){
                jcount = j - 1
                for (k in 1:jcount){
                    diff = 0
                    diff = abs(locus[j]-locus[k])
                    locus_sum = locus_sum + diff}
            }
            locus_sum = (2*locus_sum)/(length(locus)^2) 
            pi[i] = locus_sum}
        else{pi[i] = NA}
    }
    write(pi,file=paste(suf,"pi",sep="."),ncol = 1)
    meanPi = mean(pi, na.rm = TRUE)
    return(meanPi)
}

##THIS MODULE CALCULATES FREQUENCIES AND FST FOR TEST
#Counts and mean frequencies of alleles    
getCoverFst <- function(geno_matrix, suf,test_sub,null_sub,mincov){
    L <- ncol(geno_matrix)
    S <- nrow(geno_matrix)
    geno_means = vector(mode = "numeric", length = L)
    geno_cover = vector(mode = "numeric", length = L)
    geno_cover_p1 = vector(mode = "numeric", length = length(test_sub))
    geno_cover_p2 = vector(mode = "numeric", length = L - length(test_sub))
    ##THIS MODULE CALCULATES major allele frequency (AF) FOR EACH LOCUS
    af = {}
    for (i in 1:L){
        L_maf = sum(!is.na(geno_matrix[,i]))
        if(L_maf >= S*mincov){
            L_sum = sum(geno_matrix[,i],na.rm=TRUE)
            af[i] = L_sum/(2*L_maf)
        }
        else{ af[i] = NA}
    }
    #Converts AF to minor allele frequency (MAF)
    maf = {}
    for (i in 1:L){
        if (!is.na(af[i])){
            if (af[i] > 0.5){
                maf[i] = 1 - af[i]} 
            else{ maf[i] = af[i] }
        }
        else{ maf[i] = NA}
    }
    #WRITE MAF and AF TO FILE
    write(maf,file=paste(suf,"maf",sep="."), ncolumns=1)
    write(af, file=paste(suf,"af",sep="."), ncolumns=1)
    #NUMBER OF INDIVIDUALS WITH SNP (FOR REMOVING SINGLETONS)
    allele_class = {}
    for (i in 1:L){
        if(!is.na(maf[i])){
            allele_class[i] = maf[i]*geno_cover[i]}
        else{ allele_class[i] = NA} }
    write(allele_class,file=paste(suf,"AC",sep="."),ncolumns=1)
    # FST calculation
    for (i in 1:L){
        geno_means[i] = mean(geno_matrix[,i],na.rm=TRUE)/2
    }
    for (i in 1:L){
        geno_cover[i] = sum(!is.na(geno_matrix[,i]))
        geno_cover_p1[i] = sum(!is.na(geno_matrix[test_sub,i]))
        geno_cover_p2[i] = sum(!is.na(geno_matrix[-test_sub,i]))
    }
    #Use test_sub to look at Fst values
    pop_1 = {}
    pop_2 = {}
    pop_mean = {}
    FST = {}
    for (i in 1:L){
        if(geno_cover_p1[i] >= 2 && geno_cover_p2[i] >= 2 && geno_cover[i] >= S*mincov){
            if(geno_means[i] != 1 && geno_means[i] != 0){
                pop_1[i] = mean(geno_matrix[test_sub,i],na.rm=TRUE)/2
                pop_2[i] = mean(geno_matrix[-test_sub,i],na.rm=TRUE)/2
                #  pop_mean[i] = (pop_1[i] + pop_2[i])/2
                FST[i] = (geno_cover_p1[i]*(abs(geno_means[i] - pop_1[i])^2)+ geno_cover_p2[i]
                          *abs(geno_means[i] - pop_2[i])^2)/(geno_cover[i]*geno_means[i]*(1 - geno_means[i]))
            }
            else{ FST[i] = 0 }
        }
        else{ FST[i] = NA }
    }
    FST = abs(FST)
    meanFST = mean(FST, na.rm = TRUE)
    ##WRITE FST TO FILE
    write(FST,file=paste(suf,"fst",sep="."),ncolumns=1)
    write(geno_cover,file=paste(suf,"genocover",sep="."),ncolumns=1)
    
    #NULL FST
    geno_cover_p1_N = vector(mode = "numeric", length = length(null_sub))
    geno_cover_p2_N = vector(mode = "numeric", length = L - length(null_sub))
    for (i in 1:L){
        geno_cover[i] = sum(!is.na(geno_matrix[,i]))
        geno_cover_p1_N[i] = sum(!is.na(geno_matrix[null_sub,i]))
        geno_cover_p2_N[i] = sum(!is.na(geno_matrix[-null_sub,i]))
    }
    #Use null_sub to look at Fst values
    pop_1_N = {}
    pop_2_N = {}
    pop_mean = {}
    FST_N = {}
    for (i in 1:L){
        if(geno_cover_p1[i] >= 2 && geno_cover_p2[i] >= 2 && geno_cover[i] >= S*mincov){
            if(geno_means[i] != 1 && geno_means[i] != 0){
                pop_1_N[i] = mean(geno_matrix[null_sub,i],na.rm=TRUE)/2
                pop_2_N[i] = mean(geno_matrix[-null_sub,i],na.rm=TRUE)/2
                #  pop_mean[i] = (pop_1[i] + pop_2[i])/2
                FST_N[i] = (geno_cover_p1_N[i]*(abs(geno_means[i] - pop_1_N[i])^2) + geno_cover_p2_N[i]*abs(geno_means[i] - pop_2_N[i])^2)/(geno_cover[i]*geno_means[i]*(1 - geno_means[i]))
            }
            else{FST_N[i] = 0 }
        }
        else{ FST_N[i] = NA}
        FST_N = abs(FST_N)}
    #WRITE FST NULL TO FILE
    write(FST_N,file=paste(suf,"null_fst",sep="."),ncolumns=1)
    return(meanFST)
}
##THIS MODULE CALCULATES H-W EQUILIBRIUM FOR EACH LOCUS
getHW <- function(geno_matrix, suf,mincov){
    af.table <- read.table(file=paste(suf,"af",sep="."), header=FALSE)
    af <- as.numeric(af.table$V1)
    HW_chi = {}
    L <- ncol(geno_matrix)
    S <- nrow(geno_matrix)
    for (i in 1:L){
        L_maf = sum(!is.na(geno_matrix[,i]))
        if(L_maf >= S*mincov){
            p = af[i]
            q = 1 - p
            exp_pp = L_maf*(p^2)
            exp_pq = L_maf*(2*p*q)
            exp_qq = L_maf*(q^2)
            act_pp = sum(geno_matrix[,i] == 2,na.rm=TRUE)
            act_pq = sum(geno_matrix[,i] == 1,na.rm=TRUE)
            act_qq = sum(geno_matrix[,i] == 0,na.rm=TRUE)
            HW_chi[i] = ((act_pp - exp_pp)^2/exp_pp) + ((act_pq - exp_pq)^2/exp_pq) + ((act_qq - exp_qq)^2/exp_qq)
        }
        else{HW_chi = NA}
    }
    #WRITE HARDY-WEINBERG CHI-SQUARED STAT (CAN DERIVE p-values FROM THIS)
    write(HW_chi,file=paste(suf,"HW",sep="."),ncolumns=1)
}

##Code to iterate through all pairs of populations in a given pop file 
##to calculate pi and fst for each locus. Outputs to ".pi" and ".fst" files
poptable <- read.table("../../Pop2Int.txt", header = FALSE)
pops <- unique(as.vector(poptable$V3))
suf <- "over10k-min75H32u"
mincov = 0.75

fst.df <- data.frame(matrix(data=NA,nrow=length(pops),ncol=length(pops)))
rownames(fst.df) <- pops
colnames(fst.df) <- pops
pi.df <- data.frame(matrix(data=NA,nrow=length(pops),ncol=length(pops)))
rownames(pi.df) <- pops
colnames(pi.df) <- pops
popcombos <- unlist(combn(pops,2))
for(i in seq(1,length(popcombos), by = 2)){
    popsuf = paste(suf,popcombos[i],popcombos[i+1],sep="_")
    if(!file.exists(paste(popsuf,"geno",sep = "."))){
        popsuf = paste(suf,popcombos[i+1],popcombos[i],sep="_")}
    glist = makeGenoMatrix(paste(popsuf,"geno",sep = "."),paste(popsuf,"ind",sep="."))
    geno_matrix = glist[[1]]
    test_sub = glist[[2]]
    null_sub = glist[[3]]
    meanPi = getPi(geno_matrix,popsuf,mincov)
    meanFST = getCoverFst(geno_matrix, popsuf,test_sub, null_sub,mincov)
    fst.df[rownames(fst.df)==popcombos[i],popcombos[i+1]] <- meanFST
    pi.df[rownames(pi.df)==popcombos[i],popcombos[i+1]] <- meanPi
    getHW(geno_matrix,popsuf,mincov)
}
print("done")











## Read in .geno file
geno_data = scan("over10k-min75H32u_Barkeley_Sound_Conchaphila.geno",what = "character")
L = length(geno_data)
S = nchar(geno_data[[1]])
geno_matrix = matrix(data = NA, nrow = S, ncol = L)
for (i in 1:L){
    genostring = strsplit(geno_data[[i]], split="")
    for (j in 1:S){
        geno_matrix[j,i] = as.numeric(genostring[[1]][j])
    }
}
geno_matrix[geno_matrix == 9] = NA

##Read in .ind file in order to get population assignment for genotypes.
testpop = read.table("over10k-min75H32u_Barkeley_Sound_Conchaphila.ind",header=FALSE)
test_vector = testpop$V2

test_sub = {}
for (i in 1:S){
    if (test_vector[i] == 1){
        test_sub = c(test_sub,i)
    }
}
## Creates a null vector, where individuals are randomly assigned to the 2 populations.
null_vector = sample(c(rep(1,sum(test_vector)),rep(0,S-sum(test_vector))),size=S)
null_sub = {}
for (i in 1:S){
    if (null_vector[i] == 1){
        null_sub = c(null_sub,i)
    }
}

##THIS MODULE CALCULATES PI (NUCLEOTIDE DIVERSITY) FOR EACH LOCUS AND WRITES TO A FILE
## Note that this is pi when considering individuals from just 2 populations.
## Will need to adapt for entire dataset.
pi = {}
for (i in 1:L){
    locus = {}
    locus = geno_matrix[,i]
    locus = locus[!is.na(locus)]
    if(length(locus) >= 2){
        locus_sum = 0
        for (j in 2:length(locus)){
            jcount = j - 1
            for (k in 1:jcount){
                diff = 0
                diff = abs(locus[j]-locus[k])
                locus_sum = locus_sum + diff
            }
        }
        locus_sum = (2*locus_sum)/(length(locus)^2) 
        pi[i] = locus_sum
    }
    else{
        pi[i] = NA
    }
}

write(pi,file="over10k-min75H32u_Barkeley_Sound_Conchaphila.pi",ncol = 1)

##THIS MODULE CALCULATES FREQUENCIES AND FST FOR TEST
#Counts and mean frequencies of alleles
geno_means = vector(mode = "numeric", length = L)
geno_cover = vector(mode = "numeric", length = L)
geno_cover_p1 = vector(mode = "numeric", length = length(test_sub))
geno_cover_p2 = vector(mode = "numeric", length = L - length(test_sub))
for (i in 1:L){
    geno_means[i] = mean(geno_matrix[,i],na.rm=TRUE)/2
}
for (i in 1:L){
    geno_cover[i] = sum(!is.na(geno_matrix[,i]))
    geno_cover_p1[i] = sum(!is.na(geno_matrix[test_sub,i]))
    geno_cover_p2[i] = sum(!is.na(geno_matrix[-test_sub,i]))
}

#Use test_sub to look at Fst values
pop_1 = {}
pop_2 = {}
pop_mean = {}
FST = {}
for (i in 1:L){
    pop_1[i] = mean(geno_matrix[test_sub,i],na.rm=TRUE)/2
    pop_2[i] = mean(geno_matrix[-test_sub,i],na.rm=TRUE)/2
    #  pop_mean[i] = (pop_1[i] + pop_2[i])/2
    FST[i] = (geno_cover_p1[i]*(abs(geno_means[i] - pop_1[i])^2) + 
    geno_cover_p2[i]*abs(geno_means[i] - pop_2[i])^2)/(geno_cover[i]*geno_means[i]*(1 - geno_means[i]))
}

FST = abs(FST)

##WRITE FST TO FILE
write(FST,file="over10k-min75H32u_Barkeley_Sound_Conchaphila.fst",ncolumns=1)
write(geno_cover,file="over10k-min75H32u_Barkeley_Sound_Conchaphila.geno_cover",ncolumns=1)

##THIS MODULE CALCULATES FREQUENCIES AND FST FOR NULL
#Counts and mean frequencies of alleles
geno_cover_p1_N = vector(mode = "numeric", length = length(null_sub))
geno_cover_p2_N = vector(mode = "numeric", length = L - length(null_sub))

for (i in 1:L){
    geno_cover[i] = sum(!is.na(geno_matrix[,i]))
    geno_cover_p1_N[i] = sum(!is.na(geno_matrix[null_sub,i]))
    geno_cover_p2_N[i] = sum(!is.na(geno_matrix[-null_sub,i]))
}

#Use null_sub to look at Fst values
pop_1_N = {}
pop_2_N = {}
pop_mean = {}
FST_N = {}
for (i in 1:L){
    pop_1_N[i] = mean(geno_matrix[null_sub,i],na.rm=TRUE)/2
    pop_2_N[i] = mean(geno_matrix[-null_sub,i],na.rm=TRUE)/2
    #  pop_mean[i] = (pop_1[i] + pop_2[i])/2
    FST_N[i] = (geno_cover_p1_N[i]*(abs(geno_means[i] - pop_1_N[i])^2) + geno_cover_p2_N[i]*abs(geno_means[i] - pop_2_N[i])^2)/(geno_cover[i]*geno_means[i]*(1 - geno_means[i]))
}

FST_N = abs(FST_N)

#WRITE FST TO FILE
write(FST_N,file="over10k-min75H32u_Barkeley_Sound_Conchaphila.fst_null",ncolumns=1)

##THIS MODULE CALCULATES major allele frequency (AF) FOR EACH LOCUS
af = {}
for (i in 1:L){
    L_maf = sum(!is.na(geno_matrix[,i]))
    L_sum = sum(geno_matrix[,i],na.rm=TRUE)
    af[i] = L_sum/(2*L_maf)
}

#Converts AF to minor allele frequency (MAF)
maf = {}
for (i in 1:L){
    if (!is.na(af[i])){
        if (af[i] > 0.5){
            maf[i] = 1 - af[i]
        } else{
            maf[i] = af[i]
        }
    }
}

#WRITE MAF and AF TO FILE
write(maf,file="over10k-min75H32u_Barkeley_Sound_Conchaphila.maf", ncolumns=1)
write(af, file="over10k-min75H32u_Barkeley_Sound_Conchaphila.af", ncolumns=1)

##THIS MODULE CALCULATES H-W EQUILIBRIUM FOR EACH LOCUS
HW_chi = {}
for (i in 1:L){
    L_maf = sum(!is.na(geno_matrix[,i]))
    p = af[i]
    q = 1 - p
    exp_pp = L_maf*(p^2)
    exp_pq = L_maf*(2*p*q)
    exp_qq = L_maf*(q^2)
    act_pp = sum(geno_matrix[,i] == 2,na.rm=TRUE)
    act_pq = sum(geno_matrix[,i] == 1,na.rm=TRUE)
    act_qq = sum(geno_matrix[,i] == 0,na.rm=TRUE)
    HW_chi[i] = ((act_pp - exp_pp)^2/exp_pp) + ((act_pq - exp_pq)^2/exp_pq) + ((act_qq - exp_qq)^2/exp_qq)
}
#WRITE HARDY-WEINBERG CHI-SQUARED STAT (CAN DERIVE p-values FROM THIS)
write(HW_chi,file="over10k-min75H32u_Barkeley_Sound_Conchaphila.HW",ncolumns=1)

#NUMBER OF INDIVIDUALS WITH SNP (FOR REMOVING SINGLETONS)
allele_class = {}
for (i in 1:L){
    allele_class[i] = maf[i]*geno_cover[i]
}

write(allele_class,file="over10k-min75H32u_Barkeley_Sound_Conchaphila.AC",ncolumns=1)

## Create matricx with all data
library(ggplot2)
pi = read.table("over10k-min75H32u_Barkeley_Sound_Conchaphila.pi")
f_null = read.table("over10k-min75H32u_Barkeley_Sound_Conchaphila.fst_null")
f_test = read.table("over10k-min75H32u_Barkeley_Sound_Conchaphila.fst")
cover_test = read.table("over10k-min75H32u_Barkeley_Sound_Conchaphila.geno_cover") 
af = read.table("over10k-min75H32u_Barkeley_Sound_Conchaphila.af")
HW_chi = read.table("over10k-min75H32u_Barkeley_Sound_Conchaphila.HW")
allele_class = read.table("over10k-min75H32u_Barkeley_Sound_Conchaphila.AC")

FST_data = cbind(f_null,f_test,cover_test,pi,af,HW_chi,allele_class)
names(FST_data) = c("F_null","F_test","Cover","Pi","AF","HW","AC")

##SUBSETS FST DATA INTO COVERAGE GROUPS
for (i in min(cover_test):max(cover_test)){
    name = paste("FST_data_cover_",i, sep = "")
    assign(name,subset(FST_data, cover_test == i))
}
#SUBSETS FST DATA INTO GROUPS WITH FST = 1.0
FST_data_HT = subset(FST_data, F_test == 1.0)
FST_data_HN = subset(FST_data, F_null == 1.0)

#SUBSETS FST DATA INTO GROUPS WITH FST < 0.1
FST_data_LT = subset(FST_data, F_test < 0.1)
FST_data_LN = subset(FST_data, F_null < 0.1)

#SUBSETS FST DATA INTO GROUPS WITH Pi > 0.8
Hi_Pi = subset(FST_data, Pi > 0.8)
Lo_Pi = subset(FST_data, Pi < 0.2)

Pi_cat = cut(FST_data$Pi,breaks=10)

#SUBSETS FST DATA TO EXCLUDE SINGLETONS (AC == 1)
FST_data_AC = subset(FST_data, allele_class > 1)
FST_data_AC_HT = subset(FST_data_AC, F_test == 1)

#Plot densities of these values
plot(density(FST_data$Pi,bw=0.01, na.rm=TRUE),xlim=c(0,1),main="Distribution of Major Allele Frequencies across Loci", xlab="Major Allele Frequency", sub="Species")
plot(density(FST_data$Cover,bw=0.5),xlim=c(4,16),main="Density Distribution of Individuals per Locus: ", xlab="Number of Individuals (Max = )",sub="N = , bw = ")

##PLOT FST AGAINST GENO_COVER
plot(density(FST_data$F_test,na.rm=TRUE,bw=.004),ylim=c(0,60),main="Density Distribution for Fst Values:  (G)",xlab="Fst Values", sub = "N =    bw = ")
lines(density(FST_data$F_null,na.rm=TRUE,bw=.004),col=2)

plot(FST_data$F_test,FST_data$Cover,col=2,main="Fst for all Loci by Number of Samples",xlab="Fst",ylab="Number of Individuals with Locus ",sub="Black = Null, Red = ")
points(FST_data$F_null,FST_data$Cover)

#PLOTS NUCLEOTIDE DIVERSITY FOR HIGH VS. LOW FST LOCI
plot(density(FST_data_LT$Pi),main="High vs. Low Fst Loci Nucleotide Diversity Distributions ()",xlim = c(0,1), xlab="Pi (Normalized)",sub="Black = Low (Fst < 0.1) [N = ] -- Red = High (Fst = 1.0) [N = ]")
lines(density(FST_data_HT$Pi),col=2)

plot(density(FST_data_HN$Cover,bw=1),xlim=c(0,161))
lines(density(FST_data_HT$Cover),col=2)

##PLOTS DENSITY DISTRIBUTION OF NUCLEOTIDE DIVERSITY FOR FST == 1.0
TOP_DENS = ggplot(FST_data_HT, aes(x = Pi))
TOP_DENS + geom_density(fill = "black") + ggtitle(expression(atop("Fst Distributions for Nucleotide Diversity Classes", atop(italic("E. burchellii & E. lucanoides (N = 85, L = 40,722)", ""))))) + labs(x = "Nucleotide Diversity", y = "Density")

##GGPLOT NUCLEOTIDE DIVERSITY CLASSES AGAINST FST FOR NULL AND TEST
p = ggplot(FST_data, aes(factor(Pi_cat), F_null))

p + geom_boxplot() + ggtitle(expression(atop("Fst Distributions for Nucleotide Diversity Classes: Null Hypothesis", 
atop(italic("E. lucanoides (N = 15, L = 63,087", ""))))) + labs(x = "Nucleotide Diversity Class", y = "Fst Value")

