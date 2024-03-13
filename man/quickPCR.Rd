\name{quickPCR}
\alias{quickPCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{quickPCR
}
\description{
Package with a series of functions designed to automate qPCR analysis.
}
%- maybe also 'usage' for other objects documented here.
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##Basic analysis
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

# Process data to produce 2^-dCT
processedData<-quickProcess(data_df = data,
                            normalizer_char = "RPS9",
                            sampleName_num = 2,
                            targetGene_num = 3,
                            CT_num = 7,
                            conditionKey_df = conditionKey)

# Calculate relative quantitative values (RQVs) using 2^-dCTs
RQV<-quickRQV(data_df = processedData, control_char = "Citric Acid Control", RQV_input_num = 10)

##Analysis with added normalization
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

# Process data to produce 2^-dCT
processedData<-quickProcess(data_df = data,
                            normalizer_char = "RPS9",
                            sampleName_num = 2,
                            targetGene_num = 3,
                            CT_num = 7,
                            conditionKey_df = conditionKey)

#Normalize 2^-dCTs to specified gene
normalized<-quickNormalize(data_df = processedData, normalizer2_char = "CDH1",twoToNeg_dCT_num = 10)

#Calculate RQVs
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
