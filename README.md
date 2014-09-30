# getLongevityCategories

R function to bin population study participants into short-, intermediate-, and long-lived categories based on the age of death of their parents

Luke C. Pilling | L.Pilling@exeter.ac.uk   
Epidemiology and Public Health Group, University of Exeter Medical School, U.K.

Function adapted from analyses by Ambarish Dutta et al. 
> “Longer Lived Parents: Protective Associations With Cancer Incidence and Overall Mortality”
> Journals Gerontol. Ser. A Biol. Sci. Med. Sci., 2013. DOI: 10.1093/gerona/glt061. PMID: 23685624

In turn based on an analysis by Cheung and Robine, describing the "normal fitting" of age at death distributions
> "Increase in common longevity and the compression of mortality: the case of Japan" 
> Population studies, 2007. DOI: 10.1080/00324720601103833. PMID: 17365875
 

## Description

This function [ getLongevityCategories() ] will return a list of values corresponding to the cut-offs that should be applied to the participants (i.e. the "offspring") to categorise them as long-, intermediate-, or short-lived based on the age at death of the parents.

The function will also create a PDF file displaying the cut-offs graphically on histograms of age-at-death of the mothers and fathers separately.

#### Essential inputs;
* ageMum == vector of age at death of mother for all participants (missings OK)
* ageDad == vector of age at death of father for all participants (missings OK)

#### Optional inputs;
* doPlots           ==  create a PDF file with plots - default is TRUE  
* pdfFilename       ==  file name for PDF output - default is "Parental_age_at_death.DATE.pdf"  
* returnCategories  ==  default=TRUE -- return vector for all participants categorising them as offspring of;  

> [0]  Both short-lived  
> [1]  One short, one intermediate  
> [2]  Both intermediate  
> [3]  One intermediate, one long  
> [4]  Both long-lived  
> [8]  Discordant (one short, one long)  
> [9]  Discordant (at least one parent died prematurely)  
> [NA] Not dead (one or both parents not record as dead)  

## Outputs

By default text is written to the screen (the summarises the number of participants categorised as premature, short, intermediate, and long-lived), a vector output is returned, and a PDF document is created.

## Example

```R
> source("getLongevityCategories.R")
> load("yourData.Rdat")

> ageMum <- yourData$mothers_age_at_death
> ageDad <- yourData$fathers_age_at_death

> longevityCategories <- getLongevityCategories( ageMum, ageDad )

__Mothers__
Premature < 56
Short-lived (56 : 72)
Intermediate (73 : 92)
Long-lived >93

__Fathers__
Premature < 46
Short-lived (46 : 65)
Intermediate (66 : 89)
Long-lived >90

PDF output created: Parental_age_at_death.2014-09-30.pdf
in directory: U:/

__Offspring categorised by age-of-death of mother__
[0]  N=3540  Short-lived
[1]  N=8547  Intermediate-lived
[2]  N=699  Long-lived
[9]  N=1364  Died prematurely
[NA] N=15850  Not dead

__Offspring categorised by age-of-death of father__
[0]  N=4014  Short-lived
[1]  N=10682  Intermediate-lived
[2]  N=873  Long-lived
[9]  N=1082  Died prematurely
[NA] N=13349  Not dead

__Offspring categorised by age-of-death of both parents__
[0]  N=870  Both short-lived
[1]  N=3916  One short, one intermediate
[2]  N=5261  Both intermediate
[3]  N=933  One intermediate, one long
[4]  N=61  Both long-lived
[8]  N=297  Discordant (one short, one long)
[9]  N=2247  Discordant (at least one parent died prematurely)
[NA] N=16415  Not dead (one or both parents not recorded as dead)
```

#### Parental_age_at_death.2014-09-30.pdf

Histograms showing the age at death of mothers and fathers with the imputed normal distribution of age-at-death and categories based on this distribution.

![Histograms of parental longevity with imputed categories](http://s2.postimg.org/6tf4fyk89/example.png)
