##Sample code for 'Learning To Adapt'
##Prepared 3/29
##Prepared by Lindsay Mico
##lindsaymico@gmail.com

##LOAD THE LIBRARIES
##########################################
library(sp)
library(pwr)
library(ggplot2)
library(spsurvey)
library(maptools)
library(bestglm)
library(Hmisc)

##SAMPLING WITH GRTS
##########################################

## read in a shapefile and plot it
file <-readShapeLines("luck-ash.shp")
plot(file)

# Read in dbf file and look at it
att <- read.dbf("luck-ash")
head(att)
##concert to lower case
names(att) <- tolower(names(att))
head(att)

##### Equal probability GRTS survey design
Equaldsgn <- list(None=list(panel=c(Panel_1=50),
        seltype='Equal'))

##Other random draws are also ok
sample(100000000,1) # run to get random seed
dsgntime <- proc.time()

Equalsites <- grts(design=Equaldsgn, 
        DesignID='EQUAL',
        type.frame='linear',
        src.frame='shapefile',
        in.shape='luck-ash',
        att.frame=att,
        prjfilename='luck-ash',
        out.shape='Equal Sites New')


##Plot the output
output <-readShapePoints("Equal Sites New.shp")
p<-plot(output)
p


##### Stratified GRTS survey design with over sample.
####################################################
# create strata variable
att$PerInt <- as.factor( att$minor2 )
levels(att$PerInt) <- list(Perennial='0', Intermittent='610')

# specify design
Stratdsgn <- list(Perennial=list(panel=c(Panel_1=50),
        seltype='Equal',
        over=50),
        Intermittent=list(panel=c(Panel_1=50),
                seltype='Equal',
                over=50))

sample(100000000,1)
dsgntime <- proc.time()

Stratsites <- grts(design=Stratdsgn,
        DesignID='STRATIFIED',
        type.frame='linear',
        src.frame='shapefile',
        in.shape='luck-ash',
        att.frame=att,
        stratum='PerInt',
        prjfilename='luck-ash',
        out.shape='Stratified Sites')

##### Unequal probability GRTS survey design with over sample and stratification
#################################################################
# create strata variable
att$PerInt <- as.factor( att$minor2 )
levels(att$PerInt) <- list(Perennial='0', Intermittent='610')

# create Strahler categories for unequal probability design
att$strahcat <- as.factor(att$strahler)
levels(att$strahcat) <- list('1st'=c('0','1'), '2nd'='2', '3rd+'=c('3','4','5') )

Unequaldsgn <- list(Perennial=list(panel=c(Panel_1=75),
        seltype='Unequal',
        caty.n=c('1st'=25 , '2nd'=25 ,'3rd+'=25 ),
        over=36),
        Intermittent=list(panel=c(Panel_1=25),							
                seltype='Unequal',
                caty.n=c('1st'=17 , '2nd'=5 ,'3rd+'=3 ),
                over=0))

sample(100000000,1) 
dsgntime <- proc.time()

Unequalsites <- grts(design=Unequaldsgn,
        DesignID='UNEQUAL',
        type.frame='linear',
        src.frame='shapefile',
        in.shape='luck-ash',
        att.frame=att,
        stratum='PerInt',
        mdcaty='strahcat',
        prjfilename='luck-ash',
        out.shape='Unequal Sites')


##### Panels for surveys over time with unequal probability GRTS 
##### (cont) survey design with over sample and stratification
############################################################
# create strata variable
att$PerInt <- as.factor( att$minor2 )
levels(att$PerInt) <- list(Perennial='0', Intermittent='610')

# create Strahler categories for unequal probability design
att$strahcat <- as.factor(att$strahler)
levels(att$strahcat) <- list('1st'=c('0','1'), '2nd'='2', '3rd+'=c('3','4','5') )

Paneldsgn <- list(Perennial=list(panel=c(Year1=17, Year2=17, YearAll=16),
        seltype='Unequal',
        caty.n=c('1st'=15 , '2nd'=15 ,'3rd+'=20 ),
        over=50),
        Intermittent=list(panel=c(YearOnce=25),						          
        seltype='Unequal',										
        caty.n=c('1st'=17 , '2nd'=5 ,'3rd+'=3 ),						
        over=0))

sample(100000000,1) # run to get random seed
dsgntime <- proc.time()

Panelsites <- grts(design=Paneldsgn,
        DesignID='UNEQUAL',
        type.frame='linear',
        src.frame='shapefile',
        in.shape='luck-ash',
        att.frame=att,
        stratum='PerInt',
        mdcaty='strahcat',
        prjfilename='luck-ash',
        out.shape='Panel Sites')

##summary of sites
addmargins(table( Panelsites$panel,  Panelsites$mdcaty, Panelsites$stratum) )

###############################################################################
###############################################################################

## Power Analysis examples
###########################################################
## Drawn from http://www.statmethods.net/stats/power.html

# What is the power of a one-tailed t-test, with a
# significance level of 0.01, 25 people in each group, 
# and an effect size equal to 0.75?

pwr.t.test(n=25,d=0.75,sig.level=.01,alternative="greater")

# Using a two-tailed test proportions, and assuming a
# significance level of 0.01 and a common sample size of 
# 30 for each proportion, what effect size can be detected 
# with a power of .75? 

pwr.2p.test(n=30,sig.level=0.01,power=0.75)

# For a one-way ANOVA comparing 5 groups, calculate the
# sample size needed in each group to obtain a power of
# 0.80, when the effect size is moderate (0.25) and a
# significance level of 0.05 is employed.

pwr.anova.test(k=5,f=.25,sig.level=.05,power=.8)


##############################################################################
##############################################################################

###Build synthetic data

##Make a data frame with wood counts
##Include restoration sites, reference sites, before and after the project (10 per bin across 5 years)
##Lets make this consistent with the stratfied PASA design above

siteID <- rep(seq(from = 1, to = 20, by=1),5)
##Restoration sites are rest = 1
rest <- rep(c(rep(1,10),rep(0,10)),5)
own <- rep(rbinom(20,1,prob=0.5),5)
year <- numeric()
for (i in 2016:2020) {
        temp <- rep(i,20)
        #print(temp)
        year <- append(year,temp)
}

##Initialize a 100 long vector
LWD <- rep(0,100)
##Restoration and Control start as equal (but with some variation around 10)
LWD[1:20] <- rnorm(20, mean = 10, sd = 1)

#Set trend
restTrend = 2 #LWD Units/year
contTrend = .5 #LWD Units/year


##Simulate trend for rest
LWD[21:30] <- LWD[1:10]+restTrend+rnorm(1, mean = 0, sd = .5)
LWD[41:50] <- LWD[21:30]+restTrend+rnorm(1, mean = 0, sd = .5)
LWD[61:70] <- LWD[41:50]+restTrend+rnorm(1, mean = 0, sd = .5)
LWD[81:90] <- LWD[61:70]+restTrend+rnorm(1, mean = 0, sd = .5)
##Simulate trend for control
LWD[31:40] <- LWD[11:20]+contTrend+rnorm(1, mean = 0, sd = .5)
LWD[51:60] <- LWD[31:40]+contTrend+rnorm(1, mean = 0, sd = .5)
LWD[71:80] <- LWD[51:60]+contTrend+rnorm(1, mean = 0, sd = .5)
LWD[91:100] <- LWD[71:80]+contTrend+rnorm(1, mean = 0, sd = .5)

##Make a data frame
df <- data.frame(siteID,rest,year,own,LWD)
df$rest <- as.factor(df$rest)

##Plot the data over the years
qplot(df$year, df$LWD, data=df, color = rest ) 


##Throw some missing data and errors into the data
for (i in 1:nrow(df)){
        if (rnorm(1, 0, 1) > 2) {df$LWD[i] <- NaN}
}

##Turns a random data point into 100
df$LWD[sample(1:nrow(df), 1)] <- 100

##Plot the data over the years
qplot(df$year, df$LWD, data=df, color = rest ) 

##CLEAN YOUR DATA
#############################################

##Run QA/QC on the data and flag the errors
errors <- df$siteID[df$LWD > 50] ; errors

##Find the missing values
df$siteID[is.na(df$LWD)==TRUE]

##Get the means for restoration vs control
##Note the na.rm = FALSE
summarize(df$LWD, rest, mean, na.rm=FALSE)
##Now we change it to TRUE 
summarize(df$LWD, rest, mean, na.rm=TRUE)
##Note that for unequal probability distribution samples you would use the funcntion > weighted.mean()

## Impute the missing values to the mean
## This doesn't change the mean BUT it artifically decreases the variance
## Consider adjusting n if you do this a lot
df[is.nan(df$LWD)==TRUE,]$LWD <- mean(df$LWD, na.rm=TRUE)
summarize(df$LWD, rest, mean, na.rm=FALSE)

##T-Testing
##############################################
##Test restoration vs control
t.test(df$LWD~df$rest)

##Subset the data to look only at 2020
df2 <- df[df$year==2020,]
t.test(df2$LWD~df2$rest)


##ANOVAS
#############################################
df<-as.data.frame(df)
fit <- aov(df$LWD ~ df$rest*df$year, data=df)
summary(fit)

##LINEAR REGRESSION
#######################################################
##Post restoration equals 10% increase per year + random noise
##Post reference equals 5% increase per year + random noise
##Note that rest = 1 is a restoration site.
model <- lm(LWD ~ year + rest + own ,data = df)
summary(model)

modelRest <- lm(LWD ~ year,data = subset(df, rest==1))
summary(modelRest)
modelCont <- lm(LWD ~ year,data = subset(df, rest==0))
summary(modelCont)

##Use this to examine the fit of the regression
plot(model)














