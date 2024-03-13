\name{quickAU}
\alias{quickAU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{quickAU
}
\description{
Takes in a dataframe produced by quickProcess() or quickNormalize() and adds a column with arbitrary units (AUs).
}
\usage{
quickAU(data_df,
        AU_input_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_df}{The dataframe output by the quickProcess() or quickNormalize() function.
}
  \item{AU_input_num}{Numeric value specifying the column with the data to be used for calculating AUs. If running immediately after quickProcess(), this value should be 10. If running immediately after quickNormalize(), this value should be 11.
}
}
\details{
Function calculates AUs by multiplying 2^-dCTs (AU_input_num) by 1000.
}
\value{
Outputs the data from data_df with an additional column with your AUs.
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
##Example if using data directly from quickProcess()
head(processedData)
  Sample_Name           Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT
2         S10 Citric Acid Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741
3         S10 Citric Acid Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589
4         S10 Citric Acid Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812
5         S10 Citric Acid Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231
1         S10 Citric Acid Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA
7         S11 Citric Acid Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657


# The function is currently defined as
AU_noNormalize<-quickAU(data_df = processedData,
                        AU_input_num = 10)

##Example if using data directly from quickNormalize()
head(normalized)
  Sample_Name           Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT CDH1_Normalized
2         S10 Citric Acid Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741    0.0012476309
3         S10 Citric Acid Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589    0.2850707837
4         S10 Citric Acid Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812    0.1676849304
5         S10 Citric Acid Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231    1.0000000000
1         S10 Citric Acid Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA              NA
7         S11 Citric Acid Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657    0.0008374836

# The function is currently defined as
AU_Normalize<-quickAU(data_df = normalized,
                        AU_input_num = 11)
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
