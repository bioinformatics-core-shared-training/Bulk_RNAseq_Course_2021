#### Data import
## load libraries

library( tximport )
library(DESeq2)
library(tidyverse)

### Data set
# Status : Infected and Uninfected
# Time point : day 11 and day 33
# 3 reps /group
# Total of 12 samples

### Read sample metadata
sampleinfo <- read_tsv( 'data/samplesheet.tsv' )
sampleinfo

# %>% : CMD shift m or ctrl shift m
sampleinfo  %>% 
  arrange( SampleName, Status )

### Read count data
files <- str_c( 'salmon', sampleinfo$SampleName, 'quant.sf', sep='/' )
files
length(files)

# create mapping between file names and sample anmes
files <- set_names( files, sampleinfo$SampleName)
files
# read transcript and gene mapping file
tx2gene <- read_tsv( 'references/tx2gene.tsv' )
head(tx2gene)
dim(tx2gene)

# read counts data
txi <- tximport( files = files, type = 'salmon', tx2gene = tx2gene)

str( txi )
class(txi$counts)
dim(txi$counts)

dir.create( 'salmon_outputs' )

# Save txi object
saveRDS( txi, file='salmon_outputs/txi.rds')

################################################################
# dplyr
# select : for selecting columns

sampleinfo %>% 
  select( SampleName, TimePoint)

# filter : filter rows based on conditions
sampleinfo %>% 
  filter( Status == 'Infected' )

# rename
sampleinfo %>% 
  rename( time_point = TimePoint )

### data exploration
rawCounts <- txi$counts
class(rawCounts)
dim(rawCounts)
head(rawCounts)

rawCounts <- round( rawCounts, digits = 0 )
head(rawCounts)

### Filter genes
dim(rawCounts)

keep <- rowSums( rawCounts ) > 5

table(keep, useNA = 'always')

filtCounts <- rawCounts[ keep, ]
dim(rawCounts)
dim(filtCounts)

head(rawCounts[ keep, c('SRR7657878', 'SRR7657883')])

### Data exploration

summary(filtCounts)

# boxplot : raw counts 
boxplot( filtCounts, main = 'Raw counts', las = 2)

## raw counts : mean var (sd) relationship

plot( rowMeans(filtCounts), rowSds(filtCounts))

# zoom
plot( rowMeans(filtCounts), rowSds(filtCounts),
      xlim=c(0,10000),
      ylim=c(0,5000),
      main = 'Raw counts : mean vs SD'
      )

### Data trasformation

## log2 transformation

logcounts <- log2(filtCounts + 1)

summary(logcounts)

# boxplot : log2 counts
statusCols <- str_replace_all( sampleinfo$Status, c(Infected = 'red', Uninfected = 'orange'))
boxplot(logcounts, col = statusCols, las = 2, main='log2 transformation')
abline( h = median(logcounts), col= 'blue')

# logcounts : mean vs SD

plot( x=rowMeans(logcounts), rowSds( logcounts), main='log transformation : Mean Vs SD' )

# VST
vst_counts <- vst(filtCounts)

# boxplot : VST
boxplot( vst_counts, las=2, main='VST transformation')

# VST : mean vs SD 
plot( x = rowMeans(vst_counts), y = rowSds(vst_counts), ylim=c(0,5), main='VST : mean vs SD')

### Challenge 1
rlogcounts <- rlog(filtCounts)
boxplot(rlogcounts, col=statusCols, las=2, main='rlog transformation')
abline( h = median( rlogcounts), col='blue')

### PCA analysis
dim(filtCounts)

library( ggfortify )

rlogcounts <- rlog( filtCounts)
head(rlogcounts)
pcdat <- prcomp( t(rlogcounts) )

class(pcdat)
autoplot(pcdat)

# color and shape
autoplot( pcdat, 
          data = sampleinfo,
          colour = 'Status',
          shape = 'TimePoint',
          size=5
          
          )

# add the labels on PCA
library( ggrepel )
autoplot( pcdat, 
          data = sampleinfo,
          colour = 'Status',
          shape = 'TimePoint',
          size=5
) +
  geom_text_repel( mapping = aes( x = PC1, y = PC2, label = SampleName ), box.padding = 1 )

## correct sample info
sampleinfo_corrected <- sampleinfo %>% 
  mutate( Status = case_when( SampleName == 'SRR7657873' ~ 'Infected',
                              SampleName == 'SRR7657882' ~ 'Uninfected',
                              TRUE ~ Status
                              ) )
write_tsv( x=sampleinfo_corrected, file='results/Sampleinfo_Corrected.txt' )

# 
autoplot( pcdat, 
          data = sampleinfo_corrected,
          colour = 'Status',
          shape = 'TimePoint',
          size=5
)


