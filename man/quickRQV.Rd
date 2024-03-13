\name{quickRQV}
\alias{quickRQV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{quickRQV
}
\description{
Takes in a dataframe produced by quickProcess() or quickNormalize() and adds a column with relative quantitative values (RQVs).
}
\usage{
quickRQV(data_df,
         control_char,
         RQV_input_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_df}{The dataframe output by the quickProcess() or quickNormalize() function.
}
  \item{control_char}{Character value specifying the control condition to be used for calculating RQVs.
}
  \item{RQV_input_num}{Numeric value specifying the column with the data to be used for calculating RQVs. If running immediately after quickProcess(), this value should be 10. If running immediately after quickNormalize(), this value should be 11.
}
}
\details{
Function calculates RQV by finding the average of the specified values (RQV_input_num) from the control condition and dividing all individual values by this average value. Average of the control RQVs should be equal to 1.
}
\value{
Outputs the data from data_df with an additional column with your RQVs. If you ran quickNormalize() first, all the RQVs for the gene specified by normalizer2_char will be 1.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jonathan W. Villanueva
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##Example with output from quickProcess()
head(processedData)
    Sample_Name      Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation        dCT twoToNeg_dCT
2           S32 0.5 ng/mL TGFB        LGR5 33.38796 33.79592 33.92883   33.70424         0.28185169 11.0094522 0.0004850926
5           S32 0.5 ng/mL TGFB       MKI67 26.94072 26.56007 26.50558   26.66879         0.23706588  3.9740054 0.0636363375
6           S32 0.5 ng/mL TGFB       OLFM4 26.55600 26.31314 26.08473   26.31796         0.23567181  3.6231734 0.0811551581
111         S32 0.5 ng/mL TGFB        CDH1 23.61757 23.58997       NA   23.60377         0.01951837  0.9089851 0.5325595865
1           S32 0.5 ng/mL TGFB        RPS9 22.67593 22.72955 22.67887   22.69478         0.03014343         NA           NA
4           S40   1 ng/mL TGFB        LGR5 32.32980 32.31744 32.46043   32.36922         0.07922760 10.1321055 0.0008911121

# The function is currently defined as
RQV<-quickRQV(data_df = processedData, control_char = "Citric Acid Control", RQV_input_num = 10)

##Example with output from quickNormalize()
head(normalized)
    Sample_Name           Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized
66          S10 Citric Acid Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309
85          S10 Citric Acid Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589    0.2850707837
86          S10 Citric Acid Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812    0.1676849304
103         S10 Citric Acid Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231    1.0000000000
65          S10 Citric Acid Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA              NA
70          S11 Citric Acid Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836

# The function is currently defined as
RQV2<-quickRQV(data_df = normalized, control_char = "Citric Acid Control", RQV_input_num = 11)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
