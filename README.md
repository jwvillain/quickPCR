# quickPCR
Collection of functions to expedite qPCR analysis.

## Installation

``` r
# install.packages("devtools") 
devtools::install_github("jwvillain/quickPCR")

# for README examples
library(quickPCR)
library(ggplot2)
```

## Import data and perform initial data processing
``` r
#Import data:
data<-read.table('qPCR data2.txt', na.strings = "",fill = TRUE,header = T)
conditionKey<-read.csv('conditionKey2.csv')

#view data:
head(data)
  Well Sample_Name Target_Name    Task Reporter Quencher       CT
1   B4         S10        RPS9 UNKNOWN      VIC  NFQ-MGB 23.47855
2   B4         S10        LGR5 UNKNOWN      FAM  NFQ-MGB 34.50635
3   B5         S10        RPS9 UNKNOWN      VIC  NFQ-MGB 23.50337
4   B5         S10        LGR5 UNKNOWN      FAM  NFQ-MGB 34.58119
5   B6         S10        RPS9 UNKNOWN      VIC  NFQ-MGB 23.52603
6   B6         S10        LGR5 UNKNOWN      FAM  NFQ-MGB 34.78944

#Note: 'conditionKey' must have the sample name in the first column and condition in the second. Sample names must match those in the 'data' dataframe
head(conditionKey) 
  Sample_name      Condition
1         S32 0.5 ng/mL Drug
2         S33 0.5 ng/mL Drug
3         S34 0.5 ng/mL Drug
4         S35 0.5 ng/mL Drug
5         S36 0.5 ng/mL Drug
6         S37   1 ng/mL Drug

#Initial data processing
processedData<-quickProcess(data_df = data,
                             normalizer_char = "RPS9", #name of your normalizer
                             sampleName_num = 2, #column number with your sample name information in data_df
                             targetGene_num = 3, #column number with your target gene name information in data_df
                             CT_num = 7, #column number with your CT information in data_df
                             conditionKey_df = conditionKey)
```
## Example 1: Process data and then go directly into calculating relative quantitative values (RQVs)
Format your data and calculate 2^-deltaCT values for downstream analyses.
``` r

```
## Example 2: Process data and then go directly into calculating arbitrary values (AUs)
Format your data and calculate 2^-deltaCT values for downstream analyses.
``` r

```
## Example 3: Process data, normalize to a specified gene and then calculate RQVs
Format your data and calculate 2^-deltaCT values for downstream analyses.
``` r

```
## Example 4: Process data, normalize to a specified gene and then calculate AUs
Format your data and calculate 2^-deltaCT values for downstream analyses.
``` r

```
