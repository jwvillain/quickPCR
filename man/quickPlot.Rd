\name{quickPlot}
\alias{quickPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{quickPlot
}
\description{
Takes in the dataframe produced by quickRQV() or quickAU(), and creates a ggplot to provide a quick visualization of the qPCR data.
}
\usage{
quickPlot(data_df,
          input_num,
          control_char,
          conditionOrder_vec,
          ymin_num,
          ymax_num,
          dotSize_num,
          xTitle_char,
          yTitle_char,
          axisTitleSize_num,
          axisTextSize_num,
          legendTextSize_num,
          legendTitleSize_num,
          legendTitle_char,
          legendPosition_char,
          jitter_num,
          legendSize_num,
          extraMargin_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_df}{
The dataframe output by quickRQV() or quickAU(). Dataframe should be subset so it only has one gene.
}
  \item{input_num}{
Numeric value specifying the column with the data to be plotted.
}
  \item{control_char}{
Character value specifying the control condition.
}
  \item{conditionOrder_vec}{
A vector with your conditions in the order you want them to appear along the x-axis of your plot. Default puts your control condition first and then your experimental conditions in the order they appear in your input dataframe.
}
  \item{ymin_num}{
Numeric value specifying the minimum you want along the y-axis. Default is 0.98 * lowest qPCR value.
}
  \item{ymax_num}{
Numeric value specifying the maximium you want along the y-axis. Default is 1.02 * lowest qPCR value.
}
  \item{dotSize_num}{
Numeric used to determine the dot size.
}
  \item{xTitle_char}{
Character value with label for your x-axis. Default set to "Condition"
}
  \item{yTitle_char}{
Character value with label for your y-axis. Default set to the column name for data_df[,input_num]
}
  \item{axisTitleSize_num}{
Numeric value to specify the text size of the axes titles. Default set to 35.
}
  \item{axisTextSize_num}{
Numeric value to specify the text size of the axes values. Default set to 30.
}
  \item{legendTextSize_num}{
Numeric value to specify the text size of the legend text. Default set to 26.
}
  \item{legendTitleSize_num}{
Numeric value to specify the text size of the legend title. Default set to 29.
}
  \item{legendTitle_char}{
Character value with label for your legend. Default set to "Legend"
}
  \item{legendPosition_char}{
Character value specifying the legend position. Default set to "none." Refer to "legend.position"" from ggplot for more options.
}
  \item{jitter_num}{
Numeric value specifying the horizontal jitter or horizontal range datapoints can appear within a given condition. Default set to 0.2.
}
  \item{legendSize_num}{
Numeric value specifying the size of the icons in the legend. Default set to 3.
}
  \item{extraMargin_num}{
Numeric value specifying the extra margin space you want in your plot (in cm). Default set to 1.
}
}
\details{
Function requires at least these three sources of input: data_df, input_num, and control_char. Function produces a ggplot that provides a quick visual view of the distribution for the qPCR data.
}
\value{
Outputs a ggplot visualizing the qPCR data.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jonathan W. Villanueva
}
\note{
ggplot2 must be loaded for this function to operate properly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##Example
head(RQV_noNormalize)
  Sample_Name           Condition Target_Gene      CT1      CT2      CT3 Average_CT Standard_Deviation       dCT twoToNeg_dCT       RQV
2         S10 Citric Acid Control        LGR5 34.50635 34.58119 34.78944   34.62566         0.14668779 11.123009 0.0004483741 0.3729184
3         S10 Citric Acid Control       MKI67 26.74216 26.89423 26.73264   26.78967         0.09066912  3.287024 0.1024488589 0.6623044
4         S10 Citric Acid Control       OLFM4 27.40610 27.65227 27.60735   27.55524         0.13109659  4.052591 0.0602626812 0.4378348
5         S10 Citric Acid Control        CDH1 25.07976 24.87838       NA   24.97907         0.14239835  1.476416 0.3593804231 0.6999969
1         S10 Citric Acid Control        RPS9 23.47855 23.50337 23.52603   23.50265         0.02374695        NA           NA        NA
7         S11 Citric Acid Control        LGR5 34.75470 35.09823 34.24502   34.69932         0.42929086 11.898335 0.0002619657 0.2178802

quickPlot(data_df = subset(RQV_noNormalize, RQV_noNormalize[,3] == "LGR5"),
          input_num = 11,
          control_char = "Citric Acid Control")
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
