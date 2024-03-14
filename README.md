# quickPCR
Collection of functions to expedite qPCR analysis.
<br>
<img src="https://github.com/jwvillain/quickPCR/blob/main/Figures/github figure.jpg" width="400" height="400">

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

#Initial data processing
processedData<-quickProcess(data_df = data,
                             normalizer_char = "RPS9", #name of your normalizer
                             sampleName_num = 2, #column number with your sample name information in data_df
                             targetGene_num = 3, #column number with your target gene name information in data_df
                             CT_num = 7, #column number with your CT information in data_df
                             conditionKey_df = conditionKey)
```
<details><summary>Input data</summary>
<p>
  
``` r
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
```

</p>
</details>

<details><summary>Expected output</summary>
<p>
For the normalizer gene, an "NA" value will appear for delta CT (dCT) and 2^-dCT across all samples. Additionally, this tool will only support up to 3 technical replicates (CT1, CT2, CT3) for analysis. If fewer than three replicates are used, an "NA" value will be placed in that cell of the dataframe as shown below.

Non-numeric values, which usually result from lowly expressed genes where the instrument labels those wells "undetermined," are replaced with a numeric value of 40. This replacement value can be adjusted by specifying "undetermined_num = X" when running quickProcess(). X is the numeric values you want to replace non-numeric values with (for example, can set to 50). 
  
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
Calculate RQVs

``` r
RQV_noNormalize<-quickRQV(data_df = processedData,
               control_char = "Control", #Character specifying the control condition
               RQV_input_num = 10) #Numeric specifying the column with the 2^-dCT values
```

<details><summary>Expected output</summary>
<p>
Average RQV for the control condition should equal 1 for each gene.
  
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
<br>

Calculate p-values using RQVs.

``` r

signif<-quickSignif(data_df = RQV_noNormalize,
            reference_condition_char = "Control", #Character specifying the condition you want to compare to
            test_char = "wilcox", #Character specifying statistical test (wilcox or ttest)
            data_input_num = 11) #Numeric with the values you are using to calculate p-values
```

<details><summary>Expected output</summary>
<p>
Reminder that your control condition used to calculate RQVs will have an average value of 1.

``` r

head(signif)
Average RQV for the control condition should equal 1 for each gene.
         Condition1 Condition2 Gene_Target Average_Condition1 SD_Condition1 Average_Condition2 SD_Condition2 wilcox_pvalue
1    4 ng/mL Drug    Control        LGR5          5.0780174    3.40115028                  1     1.0790082    0.01904762
2    4 ng/mL Drug    Control       MKI67          0.7719667    0.47690271                  1     0.2597114    0.25714286
3    4 ng/mL Drug    Control       OLFM4          0.7404012    0.50978977                  1     0.7748819    0.60952381
4    4 ng/mL Drug    Control        CDH1          1.4315459    0.51006291                  1     0.4228664    0.35238095
5  0.5 ng/mL Drug    Control        LGR5          0.4901223    0.36141733                  1     1.0790082    0.73015873
6  0.5 ng/mL Drug    Control       MKI67          0.7035484    0.25707187                  1     0.2597114    0.11111111

```

</p>
</details>
<br>

Calculate z-score values using RQVs

``` r
ZScore_RQV_noNormalize<-quickZScore(data_df = RQV_noNormalize,
                          data_input_num = 11) #Numeric with the values you are using to calculate Z-scores
```

<details><summary>Expected output</summary>
<p>

``` r

head(ZScore_RQV_noNormalize)
   Sample_Name    Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT        RQV    Z.Score
2          S10      Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741  0.3729184 -0.6160301
7          S11      Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657  0.2178802 -0.6780490
12         S12      Control        LGR5 33.43422 32.91566 33.26767   33.20585         0.26474682  9.951401 0.0010100198  0.8400462 -0.4291682
17         S13 4 ng/mL Drug        LGR5 29.03057 29.07294 29.10276   29.06876         0.03627536  6.349655 0.0122620626 10.1985120  3.3144346
21         S14 4 ng/mL Drug        LGR5 30.24761 30.57711 30.56170   30.46214         0.18594941  7.995673 0.0039179838  3.2586365  0.5383239
26         S15 4 ng/mL Drug        LGR5 32.18607 31.83747 32.03858   32.02070         0.17498527  9.590382 0.0012972016  1.0788989 -0.3336216

```

</p>
</details>
<br>

Generate plot of the RQV values to get a quick view of the data.

``` r
qPCR_plot<-quickPlot(data_df = subset(RQV_noNormalize, RQV_noNormalize[,3] == "LGR5"),
          input_num = 11, #Numeric with the values you are wanting to graph
          control_char = "Control") #Character specifying the control condition

ggsave("RQV_noNormalize.pdf",
       plot = qPCR_plot,
       units = 'in',
       width = 8,
       height = 8)

```

<details><summary>Expected output</summary>
<p>
Need to subset to specify one gene to plot. Red dotted line marks the average of the specified control condition.
<br>
<img src="https://github.com/jwvillain/quickPCR/blob/main/Figures/RQV_noNormalize.png" width="400" height="400">

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

