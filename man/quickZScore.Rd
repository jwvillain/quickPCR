\name{quickZScore}
\alias{quickZScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{quickZScore
}
\description{
Takes in the dataframe produced by quickRQV() or quickAU(), and calculates Z-Score values.
}
\usage{
quickZScore(data_df,
            data_input_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_df}{
The dataframe output by quickRQV() or quickAU().
}
  \item{data_input_num}{
Numeric value specifying the column with the data to be used for calculating p-values.
}
}
\details{
Function calculates the Z-score for each gene by first finding the mean and standard deviation across all conditions. For each individual sample, the Z-score is calculated for each gene by (x-mean)/sd.
}
\value{
Outputs the data from data_df with an additional column with your Z-score values.
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
##Example
head(RQV_Normalize)
  Sample_Name           Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized       RQV
2         S10 Citric Acid Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309 0.6452716
3         S10 Citric Acid Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589    0.2850707837 0.8673992
4         S10 Citric Acid Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812    0.1676849304 0.6668514
5         S10 Citric Acid Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231    1.0000000000 1.0000000
1         S10 Citric Acid Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA              NA        NA
7         S11 Citric Acid Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836 0.4331445

# The function is currently defined as
ZScore_RQV_Normalize<-quickZScore(data_df = RQV_Normalize,
                                   data_input_num = 12)
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
