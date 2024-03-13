\name{quickNormalize}
\alias{quickNormalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
quickNormalize
}
\description{
Takes in the dataframe output by the quickProcess() function and normalizes your 2^-dCTs to another gene in your gene list. For example, normalizing markers for intestinal epithelium cell populations to the epithelial marker CDH1.
}
\usage{
quickNormalize(data_df,
               normalizer2_char,
               twoToNeg_dCT_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_df}{The dataframe output by the quickProcess() function.
}
  \item{normalizer2_char}{Character value specifying which gene you are wanting to normalize your 2^-dCTs to (for example, CDH1).
}
  \item{twoToNeg_dCT_num}{Numeric value specifying which column has the 2^-dCT information.
}
}
\details{
Function divides the 2^-dCT values of your gene(s) of interest by the corresponding 2^-dCT values from the specified normalizer gene (normalizer2_char). Normalization done sample by sample. From here, the author can calculate arbitrary units or use the quickRQV() function to calculate RQVs.

Function makes the assumption that all genes in data_df (except for specified normalizers) will be normalized to the gene specified by normalizer2_char. If you are only looking to perform this function on a subset of genes, the user can subset the dataframe before running the function to only analyze the normalizer specified by normalizer2_char and the genes of interest.
}
\value{
Outputs the data from data_df with an additional column with your normalized data.
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
head(processedData)
    Sample_Name      Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation        dCT twoToNeg_dCT
2           S32 0.5 ng/mL TGFB        LGR5 33.38796 33.79592 33.92883   33.70424         0.28185169 11.0094522 0.0004850926
5           S32 0.5 ng/mL TGFB       MKI67 26.94072 26.56007 26.50558   26.66879         0.23706588  3.9740054 0.0636363375
6           S32 0.5 ng/mL TGFB       OLFM4 26.55600 26.31314 26.08473   26.31796         0.23567181  3.6231734 0.0811551581
111         S32 0.5 ng/mL TGFB        CDH1 23.61757 23.58997       NA   23.60377         0.01951837  0.9089851 0.5325595865
1           S32 0.5 ng/mL TGFB        RPS9 22.67593 22.72955 22.67887   22.69478         0.03014343         NA           NA
4           S40   1 ng/mL TGFB        LGR5 32.32980 32.31744 32.46043   32.36922         0.07922760 10.1321055 0.0008911121

## The function is currently defined as
normalized<-quickNormalize(data_df = processedData, normalizer2_char = "CDH1",twoToNeg_dCT_num = 10)
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
