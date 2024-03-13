\name{quickSignif}
\alias{quickSignif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{quickSignif
}
\description{
Takes in the dataframe produced by quickRQV() or quickAU() and performs a series of pair-wise comparisons for every gene between your experimental conditions and control. Outputs a dataframe summarizing the statistics from these comparisons.
}
\usage{
quickSignif(data_df,
            reference_condition_char,
            test_char,
            data_input_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_df}{The dataframe output by quickRQV() or quickAU().
}
  \item{reference_condition_char}{Character value specifying the sample you want to compare your experimental conditions to (i.e. your control condition).
}
  \item{test_char}{Character value specifying the statistical test you want to use in your pair-wise comparison. Supports "ttest" or "wilcox" as input. Both tests are run as non-paired and two-tailed.
}
  \item{data_input_num}{Numeric value specifying the column with the data to be used for calculating p-values.
}
}
\details{
Function assumes the input is using an output file downstream of quickProcess(). This function assumes your Condition  and Gene Target information are in the 2nd and 3rd columns respectively.
}
\value{
Outputs a dataframe with eight columns. The first two columns specify the conditions being compared. The third column specifies the target gene in the comparison. Columns four and five specify the average and standard deviation for your first condition in the comparison. Columns six and seven specify the average and standard deviation for your second condition in the comparison. The final column provides the calculated p-values using the specified statistical test.
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
head(RQV_noNormalize)
  Sample_Name           Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT       RQV
2         S10 Citric Acid Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741 0.3729184
3         S10 Citric Acid Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589 0.6623044
4         S10 Citric Acid Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812 0.4378348
5         S10 Citric Acid Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231 0.6999969
1         S10 Citric Acid Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA        NA
7         S11 Citric Acid Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657 0.2178802

# The function is currently defined as
quickSignif(data_df = RQV_noNormalize,
            reference_condition_char = "Citric Acid Control",
            test_char = "ttest",
            data_input_num = 11)

quickSignif(data_df = RQV_noNormalize,
            reference_condition_char = "Citric Acid Control",
            test_char = "wilcox",
            data_input_num = 11)
  }
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
