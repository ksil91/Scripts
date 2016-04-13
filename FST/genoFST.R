
getSamples <- function(pop1,pop2,popfile, genofile){
    indtable = read.table(popfile,header=FALSE)
    newindtable <- matrix(data = NA, nrow = nrow(indtable), ncol = 2)
    
    inds <- c()
    x <- 1
    for (i in 1:nrow(indtable)){
        if (indtable[i,2] == pop1){
            inds <- c(inds, i)
            newindtable[x,1] <- indtable[i,1]
            newindtable[x,2] <- 0
            x <- x +1
        }
        if (indtable[i,2] == pop2){
            inds <- c(inds, i)
            newindtable[x,1] <- indtable[i,1]
            newindtable[x,2] <- 1
            x <- x +1
        }
    }
    newindtable <- newindtable[1:x-1,] 
    geno_data = scan(genofile,what = "character")
    L = length(geno_data)
    S = nchar(geno_data[[1]])
    geno_matrix = matrix(data = NA, nrow = S, ncol = L)
    for (i in 1:L) {
        genostring = strsplit(geno_data[[i]], split="")
        for (j in 1:S) {
            if (j %in% inds){
                geno_matrix[j,i] = as.numeric(genostring[[1]][j])
            }
        }
    }
    geno_matrix <- geno_matrix[rowSums(is.na(geno_matrix)) != ncol(geno_matrix),]
    geno_matrix[geno_matrix == 9] = NA
    #geno_matrix <- geno_matrix[, colSums(is.na(geno_matrix)) != nrow(geno_matrix)]
    
    return(list(geno_matrix, newindtable))
}

makeFiles <- function(p1,p2,popfile,genofile,prefix){
    
    samp <- getSamples(p1,p2,popfile,genofile)
    suffix <- paste(p1,p2,".txt",sep="_")
    gm <- samp[[1]]
    nit <- samp[[2]]
    test_vector = as.numeric(nit[,2])
    S <- nrow(gm)
    L <- ncol(gm)
    
    test_sub = {}
    for (i in 1:S){
        if (test_vector[i] == 1){
            test_sub = c(test_sub,i)
        }
    }
    
    null_vector = sample(c(rep(1,sum(test_vector)),rep(0,S-sum(test_vector))),size=S)
    null_sub = {}
    for (i in 1:S){
        if (null_vector[i] == 1){
            null_sub = c(null_sub,i)
        }
    }
    

    
    ##THIS MODULE CALCULATES FREQUENCIES AND FST FOR TEST
    #Counts and mean frequencies of alleles
    geno_means = vector(mode = "numeric", length = L)
    geno_cover = vector(mode = "numeric", length = L)
    geno_cover_p1 = vector(mode = "numeric", length = length(test_sub))
    geno_cover_p2 = vector(mode = "numeric", length = length(test_sub))
    for (i in 1:L){
        geno_means[i] = 1 - mean(gm[,i],na.rm=TRUE)/2
    }
    for (i in 1:L){
        geno_cover[i] = sum(!is.na(gm[,i]))
        geno_cover_p1[i] = sum(!is.na(gm[test_sub,i]))
        geno_cover_p2[i] = sum(!is.na(gm[-test_sub,i]))}
    
    #Use test_sub to look at Fst values
    pop_1 = {}
    pop_2 = {}
    pop_mean = {}
    FST = {}
    for (i in 1:L){
        #if (geno_means[i] != 0 & geno_means[i] != 1){
        pop_1[i] = 1 - mean(gm[test_sub,i],na.rm=TRUE)/2
        pop_2[i] = 1 - mean(gm[-test_sub,i],na.rm=TRUE)/2
        FST[i] = (geno_cover_p1[i]*(abs(pop_1[i] - geno_means[i])^2) + geno_cover_p2[i]*abs(pop_2[i] - geno_means[i])^2)/(geno_cover[i]*geno_means[i]*(1 - geno_means[i]))
        #}
        #else{
        #    FST[i] = 0
        #}
    }
    FST = abs(FST)
    write(FST,file=paste("FST",prefix,suffix,sep="_"),ncol =1)
    
    #Use null_sub to look at Fst values
    geno_cover_p1_N = vector(mode = "numeric", length = length(null_sub))
    geno_cover_p2_N = vector(mode = "numeric", length = L - length(null_sub))
    
    for (i in 1:L){
        geno_cover[i] = sum(!is.na(gm[,i]))
        geno_cover_p1_N[i] = sum(!is.na(gm[null_sub,i]))
        geno_cover_p2_N[i] = sum(!is.na(gm[-null_sub,i]))}
    pop_1_N = {}
    pop_2_N = {}
    pop_mean = {}
    FST_N = {}
    for (i in 1:L){
        pop_1_N[i] = 1 - mean(gm[null_sub,i],na.rm=TRUE)/2
        pop_2_N[i] = 1 - mean(gm[-null_sub,i],na.rm=TRUE)/2
        FST_N[i] = (geno_cover_p1_N[i]*(abs(pop_1_N[i]-geno_means[i])^2) + geno_cover_p2_N[i]*abs(pop_2_N[i]- geno_means[i])^2)/(geno_cover[i]*geno_means[i]*(1 - geno_means[i]))
    }
    FST_N = abs(FST_N)
    write(FST_N,file=paste("FST_null",prefix,suffix,sep="_"),ncol=1)
    write(geno_cover,file=paste("genocover",prefix,suffix,sep="_"),ncol=1)
    return(FST)
}

##THIS MODULE CALCULATES PI (NUCLEOTIDE DIVERSITY)
for (i in c(1:19)){
    pop1 = i
    geno_matrix <- getSamplesInd(pop1,"../best.samples","Am100newbest.geno")[[1]]
    S <- nrow(geno_matrix)
    L <- ncol(geno_matrix)
    pi = {}
    for (i in 1:L){
        locus = {}
        locus = geno_matrix[,i]
        locus = locus[!is.na(locus)]
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
    
    write(pi,file=paste("pi",pop1,"Am100newbest",sep="_"),ncol = 1)
    
}


###READS IN .USNPS.GENO FILE AND SAMPLE FILE 
#Change size of fst matrix relative to the total number of populations in your dataset
fstmeans_2 <- matrix(data=NA,ncol=19,nrow=19)
for (i in c(1:18)){
    n = i +1
    for (x in c(n:19)){
        pop1 = i
        pop2 = x
        fst <- makeFiles(pop1,pop2,"../best.samples","Am100bestulLimp.usnps.geno",prefix = "Am100bestulLimp")
        fstmeans_2[i,x] <- mean(fst,na.rm=TRUE)
    }
}

for (i in c(2:ncol(fstmeans_2))){
    n = i -1
    for (j in c(1:n)){
        fstmeans_2[i,j] <- fstmeans_2[j,i]
    }
}
fstmeans_2[is.na(fstmeans_2)] <- 0
###WRITE TABLE OF MEAN FST VALUES TO FILE
write.table(fstmeans_2, "FSTmean_Am100bestulLimp.txt",row.names=FALSE,col.names=FALSE)

### Mantel test of FST and water distance
waterdist <- read.csv("~/Downloads/Ostrea Phylogeography Samples - Sheet6.csv",row.names=1,header=TRUE)
waterdist[is.na(waterdist)] <- 0
library("ape")
mantel.test(waterdist,fstmeans_2,graph=TRUE,alternative="greater")

fstmeanscorr <- fstmeans_2/(1-fstmeans_2)
plot(as.dist(waterdist),as.dist(fstmeanscorr), xlab = "Coastal distance (km)",ylab="FST/(1-FST)")
abline(lm(as.dist(fstmeanscorr)~as.dist(waterdist)),col="red",lty=2)

#' pretty plot
library(MASS)
dens <- kde2d(as.dist(waterdist),as.dist(fstmeans), n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(as.dist(waterdist), as.dist(fstmeanscorr), pch=20,cex=.5, xlab = "Coastal distance (km)",ylab="FST/(1-FST")
image(dens, col=myPal(300), add=TRUE)
abline(lm(as.dist(fstmeanscorr)~as.dist(waterdist)))
