\name{quickProcess}
\alias{quickProcess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{quickProcess
}
\description{
Takes in a dataframe of CT values and an annotation file to format qPCR data and calculate 2^-dCT for downstream analyses.
}
\usage{
quickProcess(data_df,
             normalizer_char,
             sampleName_num,
             targetGene_num,
             CT_num,
             conditionKey_df,
             undetermined_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_df}{Dataframe of CTs from qPCR instrument. At mimimum, the dataframe should have a column for sample name, target gene measured, and CTs.
  }
  \item{normalizer_char}{Character value specifying the normalizer gene (18S, RPS9, etc.).
  }
  \item{sampleName_num}{Numeric value specifying which column from data_df has the sample name information.
  }
  \item{targetGene_num}{Numeric value specifying which column from data_df has the target gene information.
  }
  \item{CT_num}{Numeric value specifying which column from data_df has the CT information.
  }
  \item{conditionKey_df}{Dataframe specifying which sample name goes with which condition. Dataframe should have two columns. This function assumes that the sample name is in the first column and the condition information is in the second column.
  }
  \item{undetermined_num}{Numeric value that will be used to replace "undetermined" values. Default set to 40.
  }
}
\details{
Function is compatible with multiplexing qPCR and different formats for data_df (as long as the minimum requirements described above are met). Also compatible with up to 3 technical replicates for plating.
}
\value{
Dataframe with 10 columns representing Sample_Name, Condition, Target_Gene, CT1, CT2 (if applicable), CT3 (if applicable), Average_CT, CT Standard_Deviation, delta CT (dCT; target gene CT - normalizer gene CT), 2^-dCT.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jonathan W. Villanueva
}
\note{
For "undetermined" values, quickProcess() will print out a warning and set the "undetermined" values to 40 CT
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Import data for analysis
head(data)
  Well Sample_Name Target_Name    Task Reporter Quencher       CT
1   A1         S32        RPS9 UNKNOWN      VIC  NFQ-MGB 22.67593
2   A1         S32        LGR5 UNKNOWN      FAM  NFQ-MGB 33.38796
3   A2         S32        RPS9 UNKNOWN      VIC  NFQ-MGB 22.72955
4   A2         S32        LGR5 UNKNOWN      FAM  NFQ-MGB 33.79592
5   A3         S32        RPS9 UNKNOWN      VIC  NFQ-MGB 22.67887
6   A3         S32        LGR5 UNKNOWN      FAM  NFQ-MGB 33.92883

#Import condition key:
head(conditionKey)
  Sample_name      Condition
1         S32 0.5 ng/mL TGFB
2         S33 0.5 ng/mL TGFB
3         S34 0.5 ng/mL TGFB
4         S35 0.5 ng/mL TGFB
5         S36 0.5 ng/mL TGFB
6         S37   1 ng/mL TGFB

## The function is currently defined as
processedData<-quickProcess(data_df = data,
                            normalizer_char = "RPS9",
                            sampleName_num = 2,
                            targetGene_num = 3,
                            CT_num = 7,
                            conditionKey_df = conditionKey)
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
