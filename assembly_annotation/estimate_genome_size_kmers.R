# estimate_genome_size_kmers.R
# Estimate haploid genome size based on kmers 
# Jessie Pelosi 
# Estimating genome size from 10KP reads passed through jellyfish

# 17 mers
spec1_17 <- read.table("10KP_17mer.histo")
plot(spec1_17[3:200,], type="l", xlab = "Coverage", ylab = "Frequency", main= "Jellyfish 17mers")
png(file="17mer_jellyfish_r.png", width = 900, height =750)
dev.off()
points(spec1_17[3:200,])
sum(as.numeric(spec1_17[3:200,1]*spec1_17[3:200,2]))

max(spec1_17[3:200,2])
#12

sum(as.numeric(spec1_17[1:10000,1]*spec1_17[1:10000,2]))/12
# 2232035601 ~2.2Gb if 3:200
# 3839603757 ~3.8 Gb if 3:10000
# 3947768519 ~3.9 Gb if 1:10000

# 21 mers
spec1_21 <- read.table("10KP_21mer_reads.histo")
plot(spec1_21[2:200,], type="l", xlab = "Coverage", ylab = "Frequency", main = "Jellyfish 21mers")
png(file="21mer_jellyfish_r.png", width = 900, height =750)
dev.off()
points(spec1_21[3:200,])
sum(as.numeric(spec1_21[3:200,1]*spec1_21[3:200,2]))

max(spec1_21[3:200,2])
#11

sum(as.numeric(spec1_21[1:10000,1]*spec1_21[1:10000,2]))/11
# 2709696107 ~2.7Gb if 3:200
# 4139095955 ~4.1 Gb if 3:10000
# 4354407842 ~4.4 Gb if 1:10000


# Estimating genome size form 10KP reads passed through kmerfreq

specfreq_17 <- read.table("10KP_kmerfreq_17.kmer.freq.stat")
plot(specfreq_17[3:200,1:2], type = "l", xlab = "Coverage", y = "Frequency", main = "KmerFreq 17mers")
points(specfreq_17[3:200,])
sum(as.numeric(specfreq_17[3:200,1]*specfreq_17[3:200,2]))

max(specfreq_17[3:200,2])
#12

sum(as.numeric(specfreq_17[1:65535,1]*specfreq_17[1:65535,2]))/12
# 2232594150 ~2.2 Gb if 3:200
# 5000615165 ~5 Gb if 3:end 
# 5108996103 ~5.1 Gb if 1:end

# fit poisson distrubtion 
singleC <- sum(as.numeric(specfreq_17[3:30,1]*specfreq_17[3:30,2]))/12
plot(1:200,dpois(1:200,12)*singleC, type = "l", col=3, lty =2, xlab = "Coverage", ylab = "Frequency", main = "KmerFreq 17mers")
lines(specfreq_17[1:200,1:2], type = "l")
png(file="17mer_kmerfreq_r.png", width = 900, height =750)
dev.off()


###############All short read data ######################

# 21 mers
spec1_21 <- read.table("../../Dissertation/CH3 Lygodium_microphyllum_genome/lygodium_genome/reads.histo")
plot(spec1_21[3:200,], type="l", xlab = "Coverage", ylab = "Frequency", main = "Jellyfish 21mers")
png(file="21mer_jellyfish_shortreads.png", width = 900, height =750)
dev.off()
points(spec1_21[3:200,])
sum(as.numeric(spec1_21[3:200,1]*spec1_21[3:200,2]))
max(spec1_21[10:10000,2]) #52x 
sum(as.numeric(spec1_21[1:10000,1]*spec1_21[1:10000,2]))/52
# 21 mer : 3.58 Gb @ 52x
 
(sum(as.numeric(spec1_21[2:75,1]*spec1_21[2:75,2]))) / (sum(as.numeric(spec1_21[2:10000,1]*spec1_21[2:10000,2])))

# 31 mers
spec1_31 <- read.table("../../Dissertation/CH3 Lygodium_microphyllum_genome/lygodium_genome/genome_analyses/shortreads/31mer_reads.histo")
plot(spec1_31[3:200,], type="l", xlab = "Coverage", ylab = "Frequency", main = "Jellyfish 31mers")
png(file="31mer_jellyfish_shortreads.png", width = 900, height =750)
dev.off()
points(spec1_31[3:200,])
sum(as.numeric(spec1_31[3:200,1]*spec1_31[3:200,2]))
max(spec1_31[10:10000,2]) #47x
sum(as.numeric(spec1_31[1:10000,1]*spec1_31[1:10000,2]))/47
# 31mer : 4.1 Gb @ 47x

# 41 mers
spec1_41 <- read.table("../../Dissertation/CH3 Lygodium_microphyllum_genome/lygodium_genome/genome_analyses/shortreads/41mer_reads.histo")
plot(spec1_41[3:200,], type="l", xlab = "Coverage", ylab = "Frequency", main = "Jellyfish 41mers")
png(file="41mer_jellyfish_shortreads.png", width = 900, height =750)
dev.off()
points(spec1_41[3:200,])
sum(as.numeric(spec1_41[3:200,1]*spec1_41[3:200,2]))
max(spec1_41[10:10000,2]) #42x
sum(as.numeric(spec1_41[1:10000,1]*spec1_41[1:10000,2]))/42
# 41mer : 4.4 Gb @ 42x 

#51 mers
spec1_51 <- read.table("shortreads/51mer_reads.histo")
plot(spec1_51[3:200,], type="l", xlab = "Coverage", ylab = "Frequency", main = "Jellyfish 51mers")
png(file="51mer_jellyfish_shortreads.png", width = 900, height =750)
dev.off()
points(spec1_51[3:200,])
sum(as.numeric(spec1_51[3:200,1]*spec1_51[3:200,2]))
max(spec1_51[4:10000,2]) #37x
sum(as.numeric(spec1_51[1:10000,1]*spec1_51[1:10000,2]))/37
# 51mer : 4.6 Gb @ 37x

#61 mer
spec1_61 <- read.table("shortreads/61mer_reads.histo")
plot(spec1_61[3:200,], type="l", xlab = "Coverage", ylab = "Frequency", main = "Jellyfish 61mers")
png(file="61mer_jellyfish_shortreads.png", width = 900, height =750)
dev.off()
points(spec1_61[3:200,])
sum(as.numeric(spec1_61[3:200,1]*spec1_61[3:200,2]))
max(spec1_61[4:10000,2]) #32x
sum(as.numeric(spec1_61[1:10000,1]*spec1_61[1:10000,2]))/32
# 61mer : 4.8 Gb @ 32x
