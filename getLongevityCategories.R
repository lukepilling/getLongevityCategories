
#####################################################################################################################
#####################################################################################################################
####
####  Luke C. Pilling | L.Pilling@exeter.ac.uk
####      Epidemiology and Public Health Group, University of Exeter Medical School, U.K.
####
####  Function adapted from analyses by Ambarish Dutta et al. 
####      “Longer Lived Parents: Protective Associations With Cancer Incidence and Overall Mortality” 
####      Journals Gerontol. Ser. A Biol. Sci. Med. Sci., 2013. DOI: 10.1093/gerona/glt061. PMID: 23685624
####
####  In turn based on an analysis by Cheung and Robine, describing the "normal fitting" of age at death distributions
####      "Increase in common longevity and the compression of mortality: the case of Japan"
####      Population studies, 2007. DOI: 10.1080/00324720601103833. PMID: 17365875
####      
####  Version 0.140930
####      
#####################################################################################################################
#####################################################################################################################
####
####  This function [ getLongevityCategories() ] will return a list of values corresponding to the
####      cut-offs that should be applied to the participants (i.e. the "offspring") to categorise them as
####      long-, intermediate-, or short-lived based on the age at death of the parents.
####
####  The function will also create a PDF file displaying the cut-offs graphically on histograms of
####      age-at-death of the mothers and fathers separately.
####
####  Essential inputs;
####      ageMum            ==  vector of age at death of mother for all participants (missings OK)
####      ageDad            ==  vector of age at death of father for all participants (missings OK)
####
####  Optional inputs;
####      doPlots           ==  create a PDF file with plots - default is TRUE
####      pdfFilename       ==  file name for PDF output - default is "Parental_age_at_death.DATE.pdf"
####      returnCategories  ==  default=TRUE -- return vector for all participants categorising them as offspring of;
####                               [0]  Both short-lived
####                               [1]  One short, one intermediate
####                               [2]  Both intermediate
####                               [3]  One intermediate, one long
####                               [4]  Both long-lived
####                               [8]  Discordant (one short, one long)
####                               [9]  Discordant (at least one parent died prematurely)
####                               [NA] Not dead (one or both parents not record as dead)
####                               
#####################################################################################################################
#####################################################################################################################

getLongevityCategories <- function(
    ageMum=stop("Required: [ageMum] a vector of numeric values for mothers age at death (missings OK)"), 
    ageDad=stop("Required: [ageDad] a vector of numeric values for fathers age at death (missings OK)"), 
    doPlots=TRUE,
    pdfFilename=paste("Parental_age_at_death.", Sys.Date(), ".pdf", sep=""),
    returnCategories=TRUE
)
{
    ##############################################################################################################################
    ## check inputs are numeric
    if (!is.numeric(ageMum))  stop("Mothers age at death [ageMum] needs to be a numeric vector")
    if (!is.numeric(ageDad))  stop("Fathers age at death [ageDad] needs to be a numeric vector")
    
    ## save "all" values for later
    ageMum.all <- ageMum
    ageDad.all <- ageDad
    
    ## remove missing values
    ageMum <- ageMum[ !is.na(ageMum)  &  ageMum >= 10 ]
    ageDad <- ageDad[ !is.na(ageDad)  &  ageDad >= 10 ]
    
    ##############################################################################################################################
    ## get summary statistics for most common age at death
    ageMum.median <- median(ageMum)
    ageMum.mode   <- Mode(ageMum)
    ageMum.sd     <- sd(ageMum)
    ageMum.n      <- length(ageMum)
    
    ageDad.median <- median(ageDad)
    ageDad.mode   <- Mode(ageDad)
    ageDad.sd     <- sd(ageDad)    
    ageDad.n      <- length(ageDad)
    
    ## create subset of ages which are just the "right"/"normal" part of the distribution, based on median age at death
    ageMum.old    <- ageMum[ageMum >= ageMum.median]
    ageDad.old    <- ageDad[ageDad >= ageDad.median]
    
    ## get summary statistics for most common age at death -- just old individuals
    ageMum.old.median <- median(ageMum.old)
    ageMum.old.mode   <- Mode(ageMum.old)
    ageMum.old.sd     <- sd(ageMum.old)
    ageMum.old.n      <- length(ageMum.old)
    
    ageDad.old.median <- median(ageDad.old)
    ageDad.old.mode   <- Mode(ageDad.old)
    ageDad.old.sd     <- sd(ageDad.old)    
    ageDad.old.n      <- length(ageDad.old)
    
    ##############################################################################################################################
    ## input for non-linear model (i.e. to fit left side of curve) is a frequency distribution of the age-at-death of the "old" group
    
    ## mothers frequency tables
    freq.table.m <- table(ageMum)
    age.m        <- as.numeric( names( freq.table.m ) ) 
    den.m        <- as.numeric( freq.table.m )
    ageMum.freq  <- data.frame(age=age.m, den=den.m)

    freq.table.m.o  <- table(ageMum.old)
    age.m.o         <- as.numeric( names( freq.table.m.o ) ) 
    den.m.o         <- as.numeric( freq.table.m.o )
    ageMum.old.freq <- data.frame(age=age.m.o, den=den.m.o)
    
    ## fathers frequency tables
    freq.table.d <- table(ageDad)
    age.d        <- as.numeric( names( freq.table.d ) ) 
    den.d        <- as.numeric( freq.table.d )
    ageDad.freq  <- data.frame(age=age.d, den=den.d)

    freq.table.d.o  <- table(ageDad.old)
    age.d.o         <- as.numeric( names( freq.table.d.o ) ) 
    den.d.o         <- as.numeric( freq.table.d.o )
    ageDad.old.freq <- data.frame(age=age.d.o, den=den.d.o)
    
    ##############################################################################################################################
    ## fit non-linear regression curve using "old" values
    ##     i.e. getting the corresponding left-side of the normal age-related death distribution that is masked by the "premature" deaths
    ageMum.old.nls <- nls(den ~ (a*exp(-0.5*((age-x0)/b)^2)),                                  ## fit non-linear regression curve
                          data = ageMum.old.freq,
                          start = list( a=ageMum.old.n, b=ageMum.old.sd, x0=ageMum.old.mode) ) ## use values calculated above
    
    ageDad.old.nls <- nls(den ~ (a*exp(-0.5*((age-x0)/b)^2)),                                  ## fit non-linear regression curve
                          data = ageDad.old.freq,
                          start = list( a=ageDad.old.n, b=ageDad.old.sd, x0=ageDad.old.mode) ) ## use values calculated above
    
    ## This gives the right hand side of the normal distribuion
    	## uncomment if you are interested in seeing the "right"/"normal" part of the distribution only from the [nls] model above
    	##    this is then fitted to the "left" part of the data to create a normal distribution using the [predict] function below
    #par(mfrow=c(2,1))
    #plot(ageMum.old.freq$age, fitted.values(ageMum.old.nls))
    #plot(ageDad.old.freq$age, fitted.values(ageDad.old.nls))
    
    ## Now to get the fitted values for the left hand side of the normal distribution
    ageMum.old.nls.predict <- predict(ageMum.old.nls, ageMum.freq) ## new data frame containing ages of death in all mothers
    ageMum.freq <- cbind(ageMum.freq, ageMum.old.nls.predict)      ## the predicted values are then attached as a column to ageMum.freq dataframe

    ageDad.old.nls.predict <- predict(ageDad.old.nls, ageDad.freq) 
    ageDad.freq <- cbind(ageDad.freq, ageDad.old.nls.predict)      
    
    ## Exclude premature deaths
    ageMum.prematureCutOff <- ( ageMum.mode - ( 2 * ageMum.sd ) )          ## premature definition: ( mode - 2*SD )
    ageMum.freq <- ageMum.freq[ageMum.freq$age > ageMum.prematureCutOff,]

    ageDad.prematureCutOff <- ( ageDad.mode - ( 2 * ageDad.sd ) )         
    ageDad.freq <- ageDad.freq[ageDad.freq$age > ageDad.prematureCutOff,]

    ##############################################################################################################################
    ## Determine cut-offs used to define 3 groups: "children of short-, intermediate-, and long-lived parents"

    ## Mode & SD of normal distribution used to determine longevity cut-offs
    ageMum.norm    <- ageMum[ ageMum > ageMum.prematureCutOff ]
    ageMum.norm.M  <- ageMum.freq$age[ ageMum.freq$ageMum.old.nls.predict == max(ageMum.freq$ageMum.old.nls.predict) ]  ##  not actual mode of the data - this is the "peak" of the computed distribution
    ageMum.norm.SD <- sd( ageMum.norm )

    ageDad.norm    <- ageDad[ ageDad > ageDad.prematureCutOff ]
    ageDad.norm.M  <- ageDad.freq$age[ ageDad.freq$ageDad.old.nls.predict == max(ageDad.freq$ageDad.old.nls.predict) ]  ##  not actual mode of the data 
    ageDad.norm.SD <- sd( ageDad.norm )
    
    ## get 4 boundaries
    ageMum.norm.1  <- min(ageMum.norm)                                  ## below this parents "prematurely" died based on the normal distribution
    ageMum.norm.2  <- round( (ageMum.norm.M - ageMum.norm.SD) - 0.5 )   ## boundary between short/intermediate group  |  -0.5 so that it is always rounded down
    ageMum.norm.3  <- round( (ageMum.norm.M + ageMum.norm.SD) + 0.5 )   ## boundary between intermediate/long group   |  +0.5 so that it is always rounded up
    ageMum.norm.4  <- max(ageMum.norm)                                  ## longest-lived parent in participant group
    ageMum.cutoffs <- c(ageMum.norm.1, ageMum.norm.2, ageMum.norm.3, ageMum.norm.4)

    ageDad.norm.1  <- min(ageDad.norm)             
    ageDad.norm.2  <- round( (ageDad.norm.M - ageDad.norm.SD) - 0.5 )   
    ageDad.norm.3  <- round( (ageDad.norm.M + ageDad.norm.SD) + 0.5 )   
    ageDad.norm.4  <- max(ageDad.norm)             
    ageDad.cutoffs <- c(ageDad.norm.1, ageDad.norm.2, ageDad.norm.3, ageDad.norm.4)
    
    ##############################################################################################################################
    ## create PDF with histograms to display the "normal" age-at-death distribution and cut-offs graphically
    if (doPlots)
    {
        pdf(pdfFilename, width=10, height=7)
    
        doParLongHistogram(ageMum, ageMum.cutoffs, ageMum.freq, ageMum.norm.M, ageMum.norm.SD, "Mothers age at death")
        doParLongHistogram(ageDad, ageDad.cutoffs, ageDad.freq, ageDad.norm.M, ageDad.norm.SD, "Fathers age at death")
    
        dev.off()
    }
    
    ##############################################################################################################################
    ## print values to screen
    cat("\n__Mothers__\n")
    cat( paste( "Premature < ", min(ageMum.norm), "\n", sep="") )
    cat( paste( "Short-lived (", min(ageMum.norm), " : ", ageMum.norm.2-1, ")", "\n",  sep="") )
    cat( paste( "Intermediate (", ageMum.norm.2, " : ", ageMum.norm.3-1, ")", "\n",  sep="") )
    cat( paste( "Long-lived >", ageMum.norm.3, "\n",  sep="") )

    cat("\n__Fathers__\n")
    cat( paste( "Premature < ", min(ageDad.norm), "\n",  sep="") )
    cat( paste( "Short-lived (", min(ageDad.norm), " : ", ageDad.norm.2-1, ")", "\n",  sep="") )
    cat( paste( "Intermediate (", ageDad.norm.2, " : ", ageDad.norm.3-1, ")", "\n",  sep="") )
    cat( paste( "Long-lived >", ageDad.norm.3, "\n\n",  sep="") )
    
    if (doPlots)  cat(paste("PDF output created: ", pdfFilename, "\n", "in directory: ", getwd(), "\n\n", sep=""))

    ##############################################################################################################################
    ## categorise offspring
    if (returnCategories)  
    {
        offspring        <- rep(NA, length(ageMum.all))
        offspring.mother <- rep(NA, length(ageMum.all))
        offspring.father <- rep(NA, length(ageMum.all))
        
        ## offspring categorisation just on mother
        offspring.mother[   ageMum.all >= ageMum.cutoffs[1]  &  ageMum.all < ageMum.cutoffs[2] ]  <- 0  ## mother was short-lived
        offspring.mother[   ageMum.all >= ageMum.cutoffs[2]  &  ageMum.all < ageMum.cutoffs[3] ]  <- 1  ## mother was intermediated-lived
        offspring.mother[   ageMum.all >= ageMum.cutoffs[3]                                    ]  <- 2  ## mother was long-lived
        offspring.mother[   ageMum.all <  ageMum.cutoffs[1]                                    ]  <- 9  ## mother died prematurely
    
        ## offspring categorisation just on father
        offspring.father[   ageDad.all >= ageDad.cutoffs[1]  &  ageDad.all < ageDad.cutoffs[2] ]  <- 0  ## father was short-lived
        offspring.father[   ageDad.all >= ageDad.cutoffs[2]  &  ageDad.all < ageDad.cutoffs[3] ]  <- 1  ## father was intermediated-lived
        offspring.father[   ageDad.all >= ageDad.cutoffs[3]                                    ]  <- 2  ## father was long-lived
        offspring.father[   ageDad.all <  ageDad.cutoffs[1]                                    ]  <- 9  ## father died prematurely

        ## offspring categorised based on mother AND father
        offspring[   ageMum.all >= ageMum.cutoffs[1]  &  ageMum.all < ageMum.cutoffs[2]    ##  both parents in short-lived category
                   & ageDad.all >= ageDad.cutoffs[1]  &  ageDad.all < ageDad.cutoffs[2] ]  <-  0
    
        offspring[   ageMum.all >= ageMum.cutoffs[2]  &  ageMum.all < ageMum.cutoffs[3]    ##  both parents in intermediate-lived category
                   & ageDad.all >= ageDad.cutoffs[2]  &  ageDad.all < ageDad.cutoffs[3] ]  <-  2
    
        offspring[   ageMum.all >= ageMum.cutoffs[3]                                       ##  both parents in long-lived category
                   & ageDad.all >= ageDad.cutoffs[3]                                    ]  <-  4
        
        
        offspring[   ageMum.all >= ageMum.cutoffs[1]  &  ageMum.all < ageMum.cutoffs[2]    ##  mum short, dad intermediate
                   & ageDad.all >= ageDad.cutoffs[2]  &  ageDad.all < ageDad.cutoffs[3] ]  <-  1
        offspring[   ageMum.all >= ageMum.cutoffs[2]  &  ageMum.all < ageMum.cutoffs[3]    ##  dad short, mum intermediate
                   & ageDad.all >= ageDad.cutoffs[1]  &  ageDad.all < ageDad.cutoffs[2] ]  <-  1
        
        offspring[   ageMum.all >= ageMum.cutoffs[2]  &  ageMum.all < ageMum.cutoffs[3]    ##  mum intermediate, dad long
                   & ageDad.all >= ageDad.cutoffs[3]                                    ]  <-  3
        offspring[   ageDad.all >= ageDad.cutoffs[2]  &  ageDad.all < ageDad.cutoffs[3]    ##  dad intermediate, mum long
                   & ageMum.all >= ageMum.cutoffs[3]                                    ]  <-  3
        
        
        offspring[   ageMum.all >= ageMum.cutoffs[1]  &  ageMum.all < ageMum.cutoffs[2]    ##  DISCORDANT: mum short, dad long
                   & ageDad.all >= ageDad.cutoffs[3]                                    ]  <-  8
        offspring[   ageDad.all >= ageDad.cutoffs[1]  &  ageDad.all < ageDad.cutoffs[2]    ##  DISCORDANT: dad short, mum long
                   & ageMum.all >= ageMum.cutoffs[3]                                    ]  <-  8
    
        offspring[   ageMum.all < ageMum.cutoffs[1]                                        ##  DISCORDANT: at least one parent premature
                   | ageDad.all < ageDad.cutoffs[1]                                     ]  <-  9
            
        ## summarise (i.e. N per category) and print to screen
        offspring.mother.freq <- table( offspring.mother )
        offspring.father.freq <- table( offspring.father )
        offspring.freq        <- table( offspring )

        cat("__Offspring categorised by age-of-death of mother__\n")
        cat( paste( "[0]  N=", offspring.mother.freq[[1]], "  Short-lived\n",  sep="") )
        cat( paste( "[1]  N=", offspring.mother.freq[[2]], "  Intermediate-lived\n",  sep="") )
        cat( paste( "[2]  N=", offspring.mother.freq[[3]], "  Long-lived\n",  sep="") )
        cat( paste( "[9]  N=", offspring.mother.freq[[4]], "  Died prematurely\n",  sep="") )
        cat( paste( "[NA] N=", table(is.na(offspring.mother))[[2]], "  Not dead\n\n",  sep="") )

        cat("__Offspring categorised by age-of-death of father__\n")
        cat( paste( "[0]  N=", offspring.father.freq[[1]], "  Short-lived\n",  sep="") )
        cat( paste( "[1]  N=", offspring.father.freq[[2]], "  Intermediate-lived\n",  sep="") )
        cat( paste( "[2]  N=", offspring.father.freq[[3]], "  Long-lived\n",  sep="") )
        cat( paste( "[9]  N=", offspring.father.freq[[4]], "  Died prematurely\n",  sep="") )
        cat( paste( "[NA] N=", table(is.na(offspring.father))[[2]], "  Not dead\n\n",  sep="") )
        
        cat("__Offspring categorised by age-of-death of both parents__\n")
        cat( paste( "[0]  N=", offspring.freq[[1]], "  Both short-lived\n",  sep="") )
        cat( paste( "[1]  N=", offspring.freq[[2]], "  One short, one intermediate\n",  sep="") )
        cat( paste( "[2]  N=", offspring.freq[[3]], "  Both intermediate\n",  sep="") )
        cat( paste( "[3]  N=", offspring.freq[[4]], "  One intermediate, one long\n",  sep="") )
        cat( paste( "[4]  N=", offspring.freq[[5]], "  Both long-lived\n",  sep="") )
        cat( paste( "[8]  N=", offspring.freq[[6]], "  Discordant (one short, one long)\n",  sep="") )
        cat( paste( "[9]  N=", offspring.freq[[7]], "  Discordant (at least one parent died prematurely)\n",  sep="") )
        cat( paste( "[NA] N=", table(is.na(offspring))[[2]], "  Not dead (one or both parents not recorded as dead)\n\n",  sep="") )
        
        return( data.frame( offspring , offspring.mother , offspring.father ) )
    }
}

##############################################################################################################################
### function to get the "mode"
Mode <- function(x)  {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

##############################################################################################################################
### function to create histogram
doParLongHistogram <- function(age, cuts, freq, M, SD, lab)  {
    
    ## Draw histogram
    brks <- length( unique(age) )
    hist(age, breaks=brks, col="lightgrey", border="lightgrey", xlab="", ylab="", main="", axes=T, xaxt='n')
    
    ## Set titles/labels for axes (not axis tickmarks which is set by subsequent axis command). line is for position of the label
    title (main="Categories of parental longevity, using \"normal age-related deaths\"", 
           xlab=lab, ylab="Frequency",  cex.lab=1.3, line=2.5)
    
    ## Set the line for the fitted normal curve from the NLS regression model
    lines(freq[,1], freq[,3], lwd=1.2)  ## 1 = age,  3 = predicted normal curve
    
    ## Change the tickmarks and fonts for the axis and tickmark labels. At argument is for where the tickmarks have to be placed. lwd is
    ##     to change the size of major tickmarks. cex.axis is to change the font of the axis labels
    axis(1, at=seq(10,150,10), label=T, tck=-0.025)
    axis(1, at=seq(10,150,5),  label=F, tck=-0.025)
    axis(1, at=seq(10,150,1),  label=F, tck=-0.007)
    
    ## Add vertical ablines at the cut-offs for the 3 groups
    abline(v=cuts[1], lty=2, col="black", lwd=1)  # min
    abline(v=cuts[2], lty=2, col="black", lwd=1)  # 
    abline(v=cuts[3], lty=2, col="black", lwd=1)  # 
    abline(v=cuts[4], lty=2, col="black", lwd=1)  # max
    
    ## Calculate values for positioning of text
    max.age.cat <- max(table(age))
    y.text <- max.age.cat - ( max.age.cat / 15 )
    y.sub  <- y.text - ( max.age.cat / 25 )
    y.sup  <- y.text + ( max.age.cat / 25 )
    
    x.premature <- cuts[1] - 1
    x.short     <- cuts[2] - ( ( cuts[2] - cuts[1] ) / 2 )
    x.inter     <- cuts[3] - ( ( cuts[3] - cuts[2] ) / 2 )
    x.long      <- cuts[4] - ( ( cuts[4] - cuts[3] ) / 2 )
    
    ## Add text. cex argument is for the times of font size from default
    text (x.premature, y.text, labels="Premature deaths", cex=1, pos=2)
    text (x.short,     y.text, labels="Short",            cex=1)
    text (x.short,     y.sub,  labels="( <M - 1SD )",     cex=0.8)
    text (x.inter,     y.sup,  labels="Mode (M)",         cex=1)
    text (x.inter,     y.text, labels="Intermediate",     cex=1)
    text (x.inter,     y.sub,  labels="( M ± 1SD )",      cex=0.8)
    text (x.long,      y.text, labels="Long",             cex=1)
    text (x.long,      y.sub,  labels="( >M + 1SD )",     cex=0.8)
    
    ## add subtitle with mode and SD
    subtitle <- paste("Normal age-related distribution Mode = ", M, ", SD = ", signif(SD, 3), sep="")
    mtext(subtitle, line=0.8, cex=0.8)

}
