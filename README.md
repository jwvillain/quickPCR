# quickPCR
Collection of functions to expedite qPCR analysis.

## Installation

``` r
# install.packages("devtools") 
devtools::install_github("jwvillain/quickPCR")
```

## Import data and perform initial data processing

``` r
# Load libraries
library(quickPCR)
library(ggplot2)

#Import data:
data<-read.table('qPCR data2.txt', na.strings = "",fill = TRUE,header = T)
conditionKey<-read.csv('conditionKey2.csv')

#view data:
#Note: imported data can be formated in different ways. Minimium for the imported data is three columns defining the 1) sample names 2) gene names 3) CT values.
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
<details><summary>Expected output</summary>
<p>
For the normalizer gene, an "NA" value will appear for delta CT (dCT) and 2^-dCT across all samples. Additionally, this tool will only support up to 3 technical replicates (CT1, CT2, CT3) for analysis. If fewer than three replicates are used, an "NA" value will be placed in that cell of the dataframe as shown below.
  
``` r
head(processedData)
  Sample_Name Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT
2         S10   Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741
3         S10   Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589
4         S10   Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812
5         S10   Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231
1         S10   Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA
7         S11   Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657
```
</p>
</details>

## Example 1: Go directly into calculating relative quantitative values (RQVs)

``` r
RQV_noNormalize<-quickRQV(data_df = processedData,
               control_char = "Control", #Character specifying the control condition
               RQV_input_num = 10) #Numeric specifying the column with the 2^-dCT values
```

<details><summary>Expected output</summary>
<p>
Average for the control condition should equal 1 for each gene.
  
``` r
head(RQV_noNormalize)
  Sample_Name Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT       RQV
2         S10   Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741 0.3729184
3         S10   Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589 0.6623044
4         S10   Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812 0.4378348
5         S10   Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231 0.6999969
1         S10   Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA        NA
7         S11   Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657 0.2178802
```
</p>
</details>

## Example 2: Go directly into calculating arbitrary values (AUs)

``` r

```

## Example 3: Normalize to a specified gene and then calculate RQVs

``` r

```

## Example 4: Normalize to a specified gene and then calculate AUs

``` r

```

