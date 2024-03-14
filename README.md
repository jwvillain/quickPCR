# quickPCR
Collection of functions to expedite qPCR analysis. See image below for workflow options.
<br>
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
<br>

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
Average RQV for the control condition should equal 1 for each gene.
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

Generate plot of the RQV values to get a quick view of the data. Use "?quickPlot" to see more information on how to customize plots.

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
<br>

## Example 2: Go directly into calculating arbitrary units (AUs)

Calculate AUs by multiplying 2^-dCT by 1000.

``` r
AU_noNormalize<-quickAU(data_df = processedData,
                          AU_input_num = 10) #Numeric specifying the column with the 2^-dCT values
```

<details><summary>Expected output</summary>
<p>
  
``` r
head(AU_noNormalize)
  Sample_Name Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT          AU
2         S10   Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741   0.4483741
3         S10   Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589 102.4488589
4         S10   Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812  60.2626812
5         S10   Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231 359.3804231
1         S10   Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA          NA
7         S11   Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657   0.2619657
```

</p>
</details>
<br>

Calculate p-values using AUs.

``` r
signif<-quickSignif(data_df = AU_noNormalize,
            reference_condition_char = "Control", #Character specifying the condition you want to compare to
            test_char = "wilcox", #Character specifying statistical test (wilcox or ttest)
            data_input_num = 11) #Numeric specifying the column with data you want to use to calculate significance
```

<details><summary>Expected output</summary>
<p>

``` r
head(signif)
        Condition1 Condition2 Gene_Target Average_Condition1 SD_Condition1 Average_Condition2 SD_Condition2 wilcox_pvalue
1    4 ng/mL Drug    Control        LGR5          6.1054954     4.0893336           1.202338      1.297333    0.01904762
2    4 ng/mL Drug    Control       OLFM4        101.9072854    70.1664114         137.637934    106.653147    0.60952381
3    4 ng/mL Drug    Control       MKI67        119.4120278    73.7699168         154.685463     40.173571    0.25714286
4    4 ng/mL Drug    Control        CDH1        734.9598500   261.8677806         513.402909    217.100838    0.35238095
5  0.5 ng/mL Drug    Control        LGR5          0.5892929     0.4345459           1.202338      1.297333    0.73015873
6  0.5 ng/mL Drug    Control       MKI67        108.8287045    39.7652814         154.685463     40.173571    0.11111111

```
</p>
</details>
<br>

Calculate z-score values using AUs

``` r
ZScore_AU_noNormalize<-quickZScore(data_df = AU_noNormalize,
                                    data_input_num = 11) #Numeric specifying the column with data you want to use to calculate Z-score values
```

<details><summary>Expected output</summary>
<p>

``` r

head(ZScore_AU_noNormalize)
   Sample_Name    Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT         AU    Z.Score
2          S10      Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741  0.4483741 -0.6160301
7          S11      Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657  0.2619657 -0.6780490
12         S12      Control        LGR5 33.43422 32.91566 33.26767   33.20585         0.26474682  9.951401 0.0010100198  1.0100198 -0.4291682
17         S13 4 ng/mL Drug        LGR5 29.03057 29.07294 29.10276   29.06876         0.03627536  6.349655 0.0122620626 12.2620626  3.3144346
21         S14 4 ng/mL Drug        LGR5 30.24761 30.57711 30.56170   30.46214         0.18594941  7.995673 0.0039179838  3.9179838  0.5383239
26         S15 4 ng/mL Drug        LGR5 32.18607 31.83747 32.03858   32.02070         0.17498527  9.590382 0.0012972016  1.2972016 -0.3336216
> 

```

</p>
</details>
<br>

Generate plot of the AU values to get a quick view of the data. Use "?quickPlot" to see more information on how to customize plots.

``` r
qPCR_plot2<-quickPlot(data_df = subset(ZScore_AU_noNormalize, ZScore_AU_noNormalize[,3] == "LGR5"),
                     input_num = 11, #Numeric specifying the column with data you want to use to generate your plot
                     control_char = "Control") #Character specifying the control condition

ggsave("AU_noNormalize.pdf",
       plot = qPCR_plot2,
       units = 'in',
       width = 8,
       height = 8)

```

<details><summary>Expected output</summary>
<p>
Need to subset to specify one gene to plot. Red dotted line marks the average of the specified control condition.
<br>
<img src="https://github.com/jwvillain/quickPCR/blob/main/Figures/AU_noNormalize.png" width="400" height="400">

</p>
</details>
<br>

## Example 3: Normalize to a specified gene and then calculate RQVs
For certain datasets, you may need to normalize your 2^-dCT values to another gene. If you are looking at changes in a cell population within a heterogenous tissue (like primary intestinal tissue), you may need to normalize to a another gene from your qPCR data. For example, if looking at changes in LGR5 stem cells from the epithelium, you may want to normalize 2^-dCT values to an epithelial marker (CDH1).
<br>

Normalize your 2^-dCT values to another gene from your qPCR data.

``` r
normalized<-quickNormalize(data_df = processedData,
                           normalizer2_char = "CDH1", #character value specifying which gene you want to normalize your 2^-dCT values to
                           twoToNeg_dCT_num = 10) #Numeric specifying the column with the data you want to process
```

<details><summary>Expected output</summary>
<p>

``` r
head(normalized)
  Sample_Name Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized
2         S10   Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309
3         S10   Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589    0.2850707837
4         S10   Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812    0.1676849304
5         S10   Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231    1.0000000000
1         S10   Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA              NA
7         S11   Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836

```

</p>
</details>
<br>

Calculate RQVs

``` r
RQV_Normalize<-quickRQV(data_df = normalized,
               control_char = "Citric Acid Control",
               RQV_input_num = 11)
```

<details><summary>Expected output</summary>
<p>
<br>
  
Average RQV for the control condition should equal 1 for each gene. For the gene you normalized to in the previous step (in this case CDH1), the normalized values will be 1 across all your samples. 
  
``` r
head(RQV_Normalize)
  Sample_Name Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized       RQV
2         S10   Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309 0.6452716
3         S10   Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589    0.2850707837 0.8673992
4         S10   Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812    0.1676849304 0.6668514
5         S10   Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231    1.0000000000 1.0000000
1         S10   Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA              NA        NA
7         S11   Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836 0.4331445
```

</p>
</details>
<br>

Calculate p-values using the gene normalized RQVs.

``` r
signif<-quickSignif(data_df = RQV_Normalize,
            reference_condition_char = "Control", #Character specifying the condition you want to compare to
            test_char = "wilcox", #Character specifying statistical test (wilcox or ttest)
            data_input_num = 12 #Numeric specifying the column with data you want to use to calculate significance
```

<details><summary>Expected output</summary>
<p>

When running this step after quickNormalize(), quickSignif() will automatically identify which gene you normalized to and remove it from the output. 

``` r
head(signif)
       Condition1 Condition2 Gene_Target Average_Condition1 SD_Condition1 Average_Condition2 SD_Condition2 wilcox_pvalue
1    4 ng/mL Drug    Control        LGR5          5.3606959     4.3643295                  1     0.7272935    0.06666667
2    4 ng/mL Drug    Control       MKI67          0.5920915     0.4730187                  1     0.3823081    0.06666667
3    4 ng/mL Drug    Control       OLFM4          0.6421271     0.4567213                  1     0.3962562    0.35238095
4  0.5 ng/mL Drug    Control        LGR5          0.7813548     0.6553367                  1     0.7272935    0.55555556
5  0.5 ng/mL Drug    Control       MKI67          0.8176996     0.3555748                  1     0.3823081    0.73015873
6  0.5 ng/mL Drug    Control       OLFM4          1.2395298     0.6083415                  1     0.3962562    0.55555556

```
</p>
</details>
<br>

Calculate z-score values using the RQVs from the normalized 2^-dCT values.

``` r
ZScore_RQV_Normalize<-quickZScore(data_df = RQV_Normalize,
                                    data_input_num = 12) #Numeric specifying the column with data you want to use to calculate Z-score values
```

<details><summary>Expected output</summary>
<p>

When running this step after quickNormalize(), quickZScore() will automatically identify which gene you normalized to and remove it from the output.

``` r
head(ZScore_RQV_Normalize)
   Sample_Name    Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized        RQV    Z.Score
2          S10      Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309  0.6452716 -0.5918466
7          S11      Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836  0.4331445 -0.6639490
12         S12      Control        LGR5 33.43422 32.91566 33.26767   33.20585         0.26474682  9.951401 0.0010100198    0.0016683323  0.8628574 -0.5178888
17         S13 4 ng/mL Drug        LGR5 29.03057 29.07294 29.10276   29.06876         0.03627536  6.349655 0.0122620626    0.0205502522 10.6285396  2.8014840
21         S14 4 ng/mL Drug        LGR5 30.24761 30.57711 30.56170   30.46214         0.18594941  7.995673 0.0039179838    0.0034887935  1.8043954 -0.1978583
26         S15 4 ng/mL Drug        LGR5 32.18607 31.83747 32.03858   32.02070         0.17498527  9.590382 0.0012972016    0.0013215090  0.6834812 -0.5788591

```

</p>
</details>
<br>

Generate plot of the gene normalized RQV values to get a quick view of the data. Use "?quickPlot" to see more information on how to customize plots.

``` r
qPCR_plot3<-quickPlot(data_df = subset(RQV_Normalize, RQV_Normalize[,3] == "LGR5"),
                     input_num = 12, #Numeric specifying the column with data you want to use to generate your plot
                     control_char = "Control") #Character specifying the control condition

ggsave("RQV_Normalize.png",
       plot = qPCR_plot3,
       units = 'in',
       width = 8,
       height = 8)

```

<details><summary>Expected output</summary>
<p>
Need to subset to specify one gene to plot. Red dotted line marks the average of the specified control condition.
<br>
<img src="https://github.com/jwvillain/quickPCR/blob/main/Figures/RQV_Normalize.png" width="400" height="400">

</p>
</details>
<br>

## Example 4: Normalize to a specified gene and then calculate AUs

For certain datasets, you may need to normalize your 2^-dCT values to another gene. If you are looking at changes in a cell population within a heterogenous tissue (like primary intestinal tissue), you may need to normalize to a another gene from your qPCR data. For example, if looking at changes in LGR5 stem cells from the epithelium, you may want to normalize 2^-dCT values to an epithelial marker (CDH1).
<br>

Normalize your 2^-dCT values to another gene from your qPCR data.

``` r
normalized<-quickNormalize(data_df = processedData,
                           normalizer2_char = "CDH1", #character value specifying which gene you want to normalize your 2^-dCT values to
                           twoToNeg_dCT_num = 10) #Numeric specifying the column with the data you want to process
```

<details><summary>Expected output</summary>
<p>

``` r
head(normalized)
  Sample_Name Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized
2         S10   Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309
3         S10   Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589    0.2850707837
4         S10   Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812    0.1676849304
5         S10   Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231    1.0000000000
1         S10   Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA              NA
7         S11   Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836

```

</p>
</details>
<br>

Calculate AUs by multiplying normalized 2^-dCT by 1000.

``` r
AU_Normalize<-quickAU(data_df = normalized,
                      AU_input_num = 11) #Numeric specifying the column with the 2^-dCT values
```

<details><summary>Expected output</summary>
<p>

For the gene you normalized to in the previous step (in this case CDH1), the normalized AU values will be 1000 across all your samples.
  
``` r
head(AU_Normalize)
  Sample_Name Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized           AU
2         S10   Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309    1.2476309
3         S10   Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589    0.2850707837  285.0707837
4         S10   Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812    0.1676849304  167.6849304
5         S10   Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231    1.0000000000 1000.0000000
1         S10   Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA              NA           NA
7         S11   Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836    0.8374836
```

</p>
</details>
<br>

Calculate p-values using the gene normalized AUs.

``` r
signif<-quickSignif(data_df = AU_Normalize,
            reference_condition_char = "Control", #Character specifying the condition you want to compare to
            test_char = "wilcox", #Character specifying statistical test (wilcox or ttest)
            data_input_num = 12 #Numeric specifying the column with data you want to use to calculate significance
```

<details><summary>Expected output</summary>
<p>

When running this step after quickNormalize(), quickSignif() will automatically identify which gene you normalized to and remove it from the output. 

``` r
head(signif)
       Condition1 Condition2 Gene_Target Average_Condition1 SD_Condition1 Average_Condition2 SD_Condition2 wilcox_pvalue
1    4 ng/mL Drug    Control        LGR5          10.364891      8.438419           1.933497       1.40622    0.06666667
2    4 ng/mL Drug    Control       OLFM4         161.467824    114.846096         251.457721      99.64168    0.35238095
3    4 ng/mL Drug    Control       MKI67         194.590895    155.457605         328.650046     125.64557    0.06666667
4  0.5 ng/mL Drug    Control        LGR5           1.510747      1.267092           1.933497       1.40622    0.55555556
5  0.5 ng/mL Drug    Control       MKI67         268.737023    116.859672         328.650046     125.64557    0.73015873
6  0.5 ng/mL Drug    Control       OLFM4         311.689348    152.972156         251.457721      99.64168    0.55555556

```
</p>
</details>
<br>

Calculate z-score values using the AUs from the normalized 2^-dCT values.

``` r
ZScore_RQV_Normalize<-quickZScore(data_df = AU_Normalize,
                                    data_input_num = 12) #Numeric specifying the column with data you want to use to calculate Z-score values
```

<details><summary>Expected output</summary>
<p>

When running this step after quickNormalize(), quickZScore() will automatically identify which gene you normalized to and remove it from the output.

``` r
head(ZScore_AU_Normalize)
   Sample_Name    Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized         AU    Z.Score
2          S10      Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309  1.2476309 -0.5918466
7          S11      Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836  0.8374836 -0.6639490
12         S12      Control        LGR5 33.43422 32.91566 33.26767   33.20585         0.26474682  9.951401 0.0010100198    0.0016683323  1.6683323 -0.5178888
17         S13 4 ng/mL Drug        LGR5 29.03057 29.07294 29.10276   29.06876         0.03627536  6.349655 0.0122620626    0.0205502522 20.5502522  2.8014840
21         S14 4 ng/mL Drug        LGR5 30.24761 30.57711 30.56170   30.46214         0.18594941  7.995673 0.0039179838    0.0034887935  3.4887935 -0.1978583
26         S15 4 ng/mL Drug        LGR5 32.18607 31.83747 32.03858   32.02070         0.17498527  9.590382 0.0012972016    0.0013215090  1.3215090 -0.5788591
```

</p>
</details>
<br>

Generate plot of the gene normalized RQV values to get a quick view of the data. Use "?quickPlot" to see more information on how to customize plots.

``` r
qPCR_plot4<-quickPlot(data_df = subset(AU_Normalize, AU_Normalize[,3] == "LGR5"),
                     input_num = 12, #Numeric specifying the column with data you want to use to generate your plot
                     control_char = "Control") #Character specifying the control condition

ggsave("AU_Normalize.png",
       plot = qPCR_plot4,
       units = 'in',
       width = 8,
       height = 8)

```

<details><summary>Expected output</summary>
<p>
Need to subset to specify one gene to plot. Red dotted line marks the average of the specified control condition.
<br>
<img src="https://github.com/jwvillain/quickPCR/blob/main/Figures/AU_Normalize.png" width="400" height="400">

</p>
</details>
<br>

## Misc
Customize your plots using quickPlot()
``` r
qPCR_plot4_custom<-quickPlot(data_df = subset(AU_Normalize, AU_Normalize[,3] == "LGR5"),
                      input_num = 12,
                      control_char = "Control",
                      conditionOrder_vec = c("Control","0.5 ng/mL Drug","1 ng/mL Drug","2 ng/mL Drug","4 ng/mL Drug"), #character vector specifying the order in which you want your conditions to appear along the x-axis
                      ymin_num = 0, #Numeric value specifying the minimum you want along the y-axis
                      ymax_num = 25, #Numeric value specifying the maximium you want along the y-axis
                      dotSize_num = 4, #Numeric used to determine the dot size
                      xTitle_char = "Treatment", #Character value with label for your x-axis.
                      yTitle_char = "LGR5 (AU)", #Character value with label for your y-axis.
                      axisTextSize_num = 26, #Numeric value to specify the text size of the axes values. 
                      axisTitleSize_num = 30, #Numeric value to specify the text size of the axes titles.
                      legendPosition_char = "right", #Character value specifying the legend position. Default set to "none."
                      legendTextSize_num = 15, #Numeric value to specify the text size of the legend text.
                      legendTitleSize_num = 18, #Numeric value to specify the text size of the legend title.
                      legendTitle_char = "Color Key", #Character value with label for your legend.
                      jitter_num = 0.25, #Numeric value specifying the horizontal jitter or horizontal range datapoints can appear within a given condition.
                      legendSize_num = 1, #Numeric value specifying the size of the icons in the legend.
                      extraMargin_num = 1.5 #Numeric value specifying the extra margin space you want in your plot (in cm).
                      )

ggsave("AU_Normalize_custom.png",
       plot = qPCR_plot4_custom,
       units = 'in',
       width = 8,
       height = 8)

```

<details><summary>Expected output</summary>
<p>

<br>
<img src="https://github.com/jwvillain/quickPCR/blob/main/Figures/AU_Normalize_custom.png" width="400" height="400">

</p>
</details>
<br>

Customize your plot using commands/tools compatible with ggplot2
``` r
library(ggsignif)
qPCR_plot4_custom2<-qPCR_plot4+
  scale_fill_manual(values = c("#8B92D2","#C3EDD9","#E01542","#95F5FF","#EB6B50"))+
  geom_signif(comparisons = list(c('Control','4 ng/mL Drug'),
                                 c('Control','0.5 ng/mL Drug'),
                                 c('Control','1 ng/mL Drug'),
                                 c('Control','2 ng/mL Drug')),
              step_increase = 0.18,
              test = "wilcox.test",
              textsize=6)+
  ylim(0,34)

ggsave("AU_Normalize_custom2.png",
       plot = qPCR_plot4_custom2,
       units = 'in',
       width = 8,
       height = 8)
```

<details><summary>Expected output</summary>
<p>

<br>
<img src="https://github.com/jwvillain/quickPCR/blob/main/Figures/AU_Normalize_custom2.png" width="400" height="400">

</p>
</details>
<br>

Create a heatmap using Z-score values
``` r
library(ComplexHeatmap)
library(circlize)
library(dplyr)

#establish color palette for heatmap
col_fun = colorRamp2(c(-3,-1.5, 0, 1.5,3), c("purple","blue", "white","orange", "red")) numeric values based off of min, max values of Z-scores
col_fun(seq(-3, 3))


#Create dataframe for pheatmap and populate it with Z-scores calculated by quickZScore()
heatmap_df<-data.frame(matrix(
  ncol = length(unique(ZScore_AU_Normalize$Sample_Name)), 
  nrow = length(unique(ZScore_AU_Normalize$Target_Gene))))

rownames(heatmap_df)<-unique(ZScore_AU_Normalize$Target_Gene)
colnames(heatmap_df)<-unique(ZScore_AU_Normalize$Sample_Name)

for(i in 1:ncol(heatmap_df)){
  temp<-subset(ZScore_AU_Normalize,
               ZScore_AU_Normalize$Sample_Name == colnames(heatmap_df)[i])
  
  for(j in 1:nrow(heatmap_df)){
    temp2<-subset(temp,
                 temp$Target_Gene == rownames(heatmap_df)[j])
    heatmap_df[j,i]<-temp2[1,13]
  }
}

#Create and format annotation file
annotation_df<-data.frame(ZScore_AU_Normalize$Sample_Name,ZScore_AU_Normalize$Condition)
annotation_df<-distinct(annotation_df)
rownames(annotation_df)<-annotation_df$ZScore_AU_Normalize.Sample_Name
annotation_df<-subset(annotation_df, select = -ZScore_AU_Normalize.Sample_Name)
colnames(annotation_df)<-"Condition"

#Set colors for annotation
colours <- list('Condition' = c('Control' = 'lightgray', 
                                '4 ng/mL Drug' = 'purple',
                                '0.5 ng/mL Drug' = "lightblue",
                                '1 ng/mL Drug' = "blue",
                                '2 ng/mL Drug' = "darkblue"))

colAnn <- HeatmapAnnotation(df = annotation_df,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

#Generate and export heatmap
png("ZScore_Heatmap.png",width=13,height=3.25,units="in",res=1200)

Heatmap(
  data.matrix(heatmap_df),
  name = "Z-Score",
  col = col_fun, 
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_row_dend = TRUE,
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  width = unit(250, "mm"),
  height = unit(50, "mm"),
  top_annotation=colAnn
  )

while (!is.null(dev.list()))  dev.off()
```

<details><summary>Expected output</summary>
<p>

<br>
<img src="https://github.com/jwvillain/quickPCR/blob/main/Figures/ZScore_Heatmap.png" width="1300" height="325">

</p>
</details>
