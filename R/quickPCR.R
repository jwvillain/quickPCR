##Background material on how to build R package: http://web.mit.edu/insong/www/pdf/rpackage_instructions.pdf

###################################################################################################
#
#  quickProcess
#
###################################################################################################
quickProcess<-function(data_df,normalizer_char, sampleName_num, targetGene_num, CT_num, conditionKey_df, undetermined_num){
  #if no value is specified for "undetermined_num", will set to 40
  if(!hasArg(undetermined_num)){undetermined_num = 40}
  #subset to remove the NTC and NRT wells - control wells are offset so CTs are set to NA
  data_df<-subset(data_df, !is.na(data_df[,CT_num]))

  #create a column that holds the sample and condition information
  data_df$sample.condition<-paste0(data_df[,sampleName_num],"-",data_df[,targetGene_num])

  #Make sure CT values are numeric
  data_df[,CT_num]<-as.numeric(data_df[,CT_num])

  #if there are NA values introduced (from undetermined values), replace with undetermined_num
  if(any(is.na(data_df[,CT_num]))){
    #print out warning
    print(paste0("Warning: Some non-numeric values detected in CT column. Converting non-numeric values",undetermined_num, " for processing."))
    print(paste0("See the following rows of the input dataframe: "))
    print(which(is.na(data_df[,CT_num])))

    #Replace NA values with undetermined_num
    NA_rows<-which(is.na(data_df[,CT_num]))
    data_df[NA_rows,CT_num] = undetermined_num
  }

  #create vectors to hold relevant information to be exported
  average_CT_vec<-vector(mode="numeric")
  stdev_CT_vec<-vector(mode="numeric")
  gene_vec<-vector(mode = "character")
  sample_vec<-vector(mode = "character")
  CT1_vec<-vector(mode="numeric")
  CT2_vec<-vector(mode = "numeric")
  CT3_vec<-vector(mode = "numeric")
  condition_vec<-vector(mode = "character")

  for(i in 1:length(unique(data_df$sample.condition))){
    #subset for a sample and a given gene
    temp<-subset(data_df, data_df$sample.condition == unique(data_df$sample.condition)[i])

    #calculate average and standard deviation
    averageCT<-mean(temp[,CT_num])
    SD_CT<-sd(temp[,CT_num])

    #store information in vectors
    average_CT_vec<-c(average_CT_vec,averageCT)
    stdev_CT_vec<-c(stdev_CT_vec,SD_CT)
    gene_vec<-c(gene_vec,temp[1,targetGene_num])
    sample_vec<-c(sample_vec,temp[1,sampleName_num])
    condition_vec<-c(condition_vec,conditionKey_df[which(conditionKey_df[,1]==temp[1,sampleName_num]),2])

    #store individual CTs based off of how many technical replicates the user ran - assumes the same CTs run for all samples
    CT1_vec<-c(CT1_vec, temp[1,CT_num])

    if(length(temp)>1){
      CT2_vec<-c(CT2_vec, temp[2,CT_num])
    } else{
      CT2_vec<-c(CT2_vec, NA)
    }

    if(length(temp)>2){
      CT3_vec<-c(CT3_vec, temp[3,CT_num])
    } else{
      CT3_vec<-c(CT3_vec, NA)
    }
  }

  #combine all the vectors into a single dataframe
  export_df<-data.frame(sample_vec, condition_vec,gene_vec,CT1_vec,CT2_vec,CT3_vec,average_CT_vec,stdev_CT_vec)
  colnames(export_df)<-c("Sample_Name", "Condition","Target_Gene","CT1","CT2","CT3","Average_CT","Standard_Deviation")

  #create a dataframe to store all the dCT and 2^-CT data across samples
  df_final<-data.frame(matrix(ncol = 10, nrow = 0))
  colnames(df_final)<-c(colnames(export_df),"dCT","twoToNeg_dCT")

  #calculate dCT and 2^-CT
  for(i in 1:length(unique(export_df$Sample_Name))){
    #subset for a sample
    sample_df<-subset(export_df, export_df$Sample_Name == unique(export_df$Sample_Name)[i])

    #separate out the normalizer
    normalizer_df<-subset(sample_df, sample_df$Target_Gene == normalizer_char)

    #separate out the genes that are not the normalizer
    geneTargets_df<-subset(sample_df, sample_df$Target_Gene != normalizer_char)

    #calculate dCT
    geneTargets_df$dCT<- geneTargets_df$Average_CT - normalizer_df[1,7]

    #Calculate 2^-dCT
    geneTargets_df$twoToNeg_dCT <- 2^(-geneTargets_df$dCT)

    #Add in NA values for normalizer df
    normalizer_df$dCT<-NA
    normalizer_df$twoToNeg_dCT<-NA

    #re-combine into new dataframe
    df_final<-data.frame(rbind(df_final,geneTargets_df))
    df_final<-data.frame(rbind(df_final,normalizer_df))
  }

  df_final
}

###################################################################################################
#
#  quickRQV
#
###################################################################################################
quickRQV<-function(data_df, control_char, RQV_input_num){
  #create a dataframe just for normalizer information
  normalizer_df<-subset(data_df, is.na(data_df$dCT))

  #Remove the RPS9 rows from data_df
  data_df<-subset(data_df, !is.na(data_df$dCT))

  #subset out your control condition
  control_df<-subset(data_df,data_df$Condition==control_char)

  #Create a dataframe to store the RQV information
  df_final<-data.frame(matrix(ncol = (ncol(data_df)+1), nrow = 0))
  colnames(df_final)<-c(colnames(data_df),"RQV")

  #Create a vector of the list of target genes looked at in your qCPR data
  geneList<-unique(data_df$Target_Gene)

  #Create a vector to store the RQV information
  RQV_vec<-vector(mode = "numeric")

  #Calculate RQVs and store in df_final
  for(i in 1:length(geneList)){
    #subset for a given gene from geneList
    temp<-subset(data_df, data_df$Target_Gene == geneList[i])
    temp_control<-subset(control_df, control_df$Target_Gene == geneList[i])

    #Calculate the average 2^-dCT for the control condition
    average_control<-mean(temp_control[,RQV_input_num])

    #Calculate RQVs for a given gene
    temp$RQV <- (temp[,RQV_input_num])/average_control

    #load RQV info into df_final
    df_final<-data.frame(rbind(df_final, temp))
  }

  #store RPS9 information into df_final
  normalizer_df$RQV<-NA
  df_final<-data.frame(rbind(df_final, normalizer_df))

  #sort the dataframe by sample condition
  df_final<-df_final[order(df_final$Sample_Name),]

  #export df_final
  df_final
}

###################################################################################################
#
#  quickNormalize
#
###################################################################################################
quickNormalize<-function(data_df, normalizer2_char,twoToNeg_dCT_num){
  #create a dataframe just for normalizer information
  normalizer_df<-subset(data_df, is.na(data_df$dCT))

  #subset to get rid of your first normalizer (RPS9, 18S) that is used to standardize cDNA concentration
  data_df<-subset(data_df, !is.na(data_df$dCT))

  #Create a dataframe just for your 2nd normalizer
  Normalizer2_df<-subset(data_df, data_df$Target_Gene == normalizer2_char)

  #Create a dataframe to store information to be consolidated
  df_final<-data.frame(matrix(ncol = (ncol(data_df)+1), nrow = 0))
  colnames(df_final)<-c(colnames(data_df),paste0(normalizer2_char,"_Normalized"))

  #Calculate normalized values and store in df_final
  for(i in 1:length(unique(data_df$Sample_Name))){
    #subset for a given sample
    temp<-subset(data_df, data_df$Sample_Name == unique(data_df$Sample_Name)[i])
    temp_normalizer2<-subset(temp, temp$Target_Gene == normalizer2_char)

    #Normalize all your experimental gene values to your 2nd normalizer
    temp$Normalized<-temp$twoToNeg_dCT/temp_normalizer2[1,twoToNeg_dCT_num]

    #Change colnames
    colnames(temp)<-c(colnames(data_df)[1:ncol(data_df)],paste0(normalizer2_char,"_Normalized"))

    #load RQV info into df_final
    df_final<-data.frame(rbind(df_final, temp))
  }

  #store RPS9 information into df_final
  normalizer_df$Normalized<-NA
  colnames(normalizer_df)<-c(colnames(data_df)[1:ncol(data_df)],paste0(normalizer2_char,"_Normalized"))
  df_final<-data.frame(rbind(df_final, normalizer_df))

  #sort the dataframe by sample condition
  df_final<-df_final[order(df_final$Sample_Name),]

  #export df_final
  df_final

}

###################################################################################################
#
#  quickAU
#
###################################################################################################
quickAU<-function(data_df, AU_input_num){
  #calculate arbitrary units (AU)
  data_df$AU<-data_df[,AU_input_num]*1000

  #sort the dataframe by sample condition
  data_df<-data_df[order(data_df$Sample_Name),]

  #export data_df
  data_df
}

###################################################################################################
#
#  quickSignif
#
###################################################################################################
quickSignif<-function(data_df, reference_condition_char,test_char,data_input_num){
  #subset data_df to remove normalizer (18S, RPS9)
  data_df<-subset(data_df, !is.na(data_df$dCT))

  #subset for a "reference" condition (i.e. the condition everything will be compared to. Can think of as the control)
  reference_df<-subset(data_df, data_df[,2]==reference_condition_char)

  #subset to make a dataframe with all the conditions, EXCEPT for your reference condition
  inquiry_df<-subset(data_df, data_df[,2] != reference_condition_char)

  #Create vectors to store relevant information
  condition1<-vector(mode = "character")
  condition2<-vector(mode = "character")
  pvalue<-vector(mode = "numeric")
  gene<-vector(mode = "character")
  average_condition1<-vector(mode = "numeric")
  sd_condition1<-vector(mode = "numeric")
  average_condition2<-vector(mode = "numeric")
  sd_condition2<-vector(mode = "numeric")

  #Perform wilcox for datasets that do not follow a normal distribution
  if(test_char == "wilcox"){
    #for each sample condition
    for(i in 1:length(unique(inquiry_df[,2]))){
      #subset for one specific sample condition
      temp_condition<-subset(inquiry_df, inquiry_df[,2] == unique(inquiry_df[,2])[i])

      #for each gene within a given sample condition
      for(j in 1:length(unique(temp_condition[,3]))){
        #subset for 1 gene
        temp_gene_inquiry_df<-subset(temp_condition, temp_condition[,3]==unique(temp_condition[,3])[j])
        reference_temp_df<-subset(reference_df, reference_df[,3]==unique(temp_condition[,3])[j])

        #below "if" statement is to filter out genes that were used for additional normalization (for example CDH1)
        if(sd(temp_gene_inquiry_df[,data_input_num]) != 0 & sd(reference_temp_df[,data_input_num]) != 0 & sd(temp_gene_inquiry_df[,data_input_num])!=sd(reference_temp_df[,data_input_num])){
          #Add information to vectors
          condition1<-c(condition1, unique(inquiry_df[,2])[i])
          condition2<-c(condition2, reference_condition_char)
          pvalue<- c(pvalue, wilcox.test(temp_gene_inquiry_df[,data_input_num],reference_temp_df[,data_input_num])$p.value)
          gene<-c(gene, unique(temp_condition[,3])[j])
          average_condition1<-c(average_condition1, mean(temp_gene_inquiry_df[,data_input_num]))
          sd_condition1<- c(sd_condition1, sd(temp_gene_inquiry_df[,data_input_num]))
          average_condition2<-c(average_condition2, mean(reference_temp_df[,data_input_num]))
          sd_condition2<-c(sd_condition2, sd(reference_temp_df[,data_input_num]))
        }
      }
    }

    final_df<-data.frame(condition1,condition2,gene,average_condition1,sd_condition1,average_condition2,sd_condition2,pvalue)
    colnames(final_df)<-c("Condition1","Condition2","Gene_Target","Average_Condition1","SD_Condition1","Average_Condition2","SD_Condition2",paste0(test_char,"_pvalue"))

    #export final_df
    data.frame(final_df)
    return(final_df)
  }

  #perform t-test for datasets that follow a normal distribution
  if(test_char == "ttest"){
    #for each sample condition
    for(i in 1:length(unique(inquiry_df[,2]))){
      #subset for one specific sample condition
      temp_condition<-subset(inquiry_df, inquiry_df[,2] == unique(inquiry_df[,2])[i])

      #for each gene within a given sample condition
      for(j in 1:length(unique(temp_condition[,3]))){
        #subset for 1 gene
        temp_gene_inquiry_df<-subset(temp_condition, temp_condition[,3]==unique(temp_condition[,3])[j])
        reference_temp_df<-subset(reference_df, reference_df[,3]==unique(temp_condition[,3])[j])

        #below "if" statement is to filter out genes that were used for additional normalization (for example CDH1)
        if(sd(temp_gene_inquiry_df[,data_input_num]) != 0 & sd(reference_temp_df[,data_input_num]) != 0 & sd(temp_gene_inquiry_df[,data_input_num])!=sd(reference_temp_df[,data_input_num])){
          #Add information to vectors
          condition1<-c(condition1, unique(inquiry_df[,2])[i])
          condition2<-c(condition2, reference_condition_char)
          pvalue<- c(pvalue, t.test(temp_gene_inquiry_df[,data_input_num],reference_temp_df[,data_input_num])$p.value)
          gene<-c(gene, unique(temp_condition[,3])[j])
          average_condition1<-c(average_condition1, mean(temp_gene_inquiry_df[,data_input_num]))
          sd_condition1<- c(sd_condition1, sd(temp_gene_inquiry_df[,data_input_num]))
          average_condition2<-c(average_condition2, mean(reference_temp_df[,data_input_num]))
          sd_condition2<-c(sd_condition2, sd(reference_temp_df[,data_input_num]))
        }
      }
    }

    final_df<-data.frame(condition1,condition2,gene,average_condition1,sd_condition1,average_condition2,sd_condition2,pvalue)
    colnames(final_df)<-c("Condition1","Condition2","Gene_Target","Average_Condition1","SD_Condition1","Average_Condition2","SD_Condition2",paste0(test_char,"_pvalue"))

    #export final_df
    data.frame(final_df)
    return(final_df)
  }

}

###################################################################################################
#
#  quickZScore
#
###################################################################################################
#Assumes data is processed using an upstream function (column 3 is gene name)
quickZScore<-function(data_df,data_input_num){
  #subset data_df to remove normalizer (18S, RPS9)
  data_df<-subset(data_df, !is.na(data_df$dCT))

  #create dataframe to store all the information
  df_final<-data.frame(matrix(ncol = (ncol(data_df)+1), nrow = 0))
  colnames(df_final)<-c(colnames(data_df),"Z.Score")

  #for each gene...
  for(i in 1:length(unique(data_df[,3]))){
    #subset for a gene of interest
    temp_df<-subset(data_df, data_df[,3]==unique(data_df[,3])[i])

    #Following if-statement is to determine whether a second normalization was used using quickNormalize(). If quickNormalize was used, the standard deviation across all your samples for the gene should be 0. If so, skips the normalizer used for quickNormalize().
    if(sd(temp_df[,data_input_num]) != 0){
      #Calculate average and standard deviation for a given gene across all samples
      sd_num<-sd(temp_df[,data_input_num])
      average_num<-mean(temp_df[,data_input_num])

      #Print average and standard deviation for the user:
      print(paste0(unique(data_df[,3])[i],": Average = ",average_num,", Standard Deviation =",sd_num))

      #calculate Z-score (x-mean)/sd
      temp_df$Z.Score<-(temp_df[,data_input_num]-average_num)/sd_num

      #Add to df_final
      df_final<-data.frame(rbind(df_final, temp_df))
    }
  }
  #export
  return(df_final)
}

###################################################################################################
#
#  quickPlot
#
###################################################################################################
quickPlot<-function(data_df,input_num, control_char, conditionOrder_vec,ymin_num,ymax_num,dotSize_num,xTitle_char,yTitle_char,axisTitleSize_num,axisTextSize_num,legendTextSize_num, legendTitleSize_num,legendTitle_char, legendPosition_char,jitter_num,legendSize_num,extraMargin_num){
  #set defaults if not established
  if(!hasArg(extraMargin_num)){extraMargin_num = 1}
  if(!hasArg(legendPosition_char)){legendPosition_char = "none"} # legend position
  if(!hasArg(ymin_num)){ymin_num = min(data_df[,input_num])*0.98} #ymin
  if(!hasArg(ymax_num)){ymax_num = max(data_df[,input_num])*1.02} #ymax
  if(!hasArg(dotSize_num)){dotSize_num = 3} #factor for dot size
  if(!hasArg(axisTitleSize_num)){axisTitleSize_num = 35} #size of axis title
  if(!hasArg(axisTextSize_num)){axisTextSize_num = 30} #size of axis text
  if(!hasArg(legendTextSize_num)){legendTextSize_num = 26} #legend text size
  if(!hasArg(legendTitleSize_num)){legendTitleSize_num = 29} #legend title text size
  if(!hasArg(legendSize_num)){legendSize_num = 3}
  if(!hasArg(xTitle_char)){xTitle_char = "Condition"}
  if(!hasArg(yTitle_char)){yTitle_char = colnames(data_df)[input_num]}
  if(!hasArg(legendTitle_char)){legendTitle_char = "Legend"}
  if(!hasArg(jitter_num)){jitter_num = 0.2}
  if(!hasArg(conditionOrder_vec)){
    #add control condition to start of vector containing all condition names
    conditionOrder_vec<-c(control_char,unique(data_df[,2]))

    #remove the second occurance of control condition
    conditionOrder_vec<-unique(conditionOrder_vec)
  } #order of conditions on x-axis

  #set the condition order to a factor
  condition_order<-factor(conditionOrder_vec)

  #Calculate average value for control condition
  control_df<-subset(data_df,data_df[,2]==control_char)
  control_average<-mean(control_df[,input_num])

  #generate plot
  ggplot(data_df,aes(x=factor(data_df[,2],level=condition_order),
                     y=data_df[,input_num],
                     fill=data_df[,2])) +
    theme_bw()+
    geom_boxplot(data = data_df,
                 width=0.5,
                 alpha = 0.25,
                 outlier.size = 0,
                 aes(x=factor(data_df[,2],level=condition_order),
                     y=data_df[,input_num],
                     fill=NULL))+
    geom_hline(yintercept = control_average,
               linetype = "dashed",
               color = "red")+
    geom_dotplot(binaxis = 'y',
                 stackdir = 'center',
                 stackratio=1.0,
                 dotsize=(ymax_num-ymin_num)*dotSize_num,
                 method='dotdensity',
                 binwidth=1/50,
                 position = position_jitter(width = jitter_num, height = 0)
    ) +
    theme(plot.title = element_text(hjust = 0.3))+
    xlab(xTitle_char) +
    ylab(yTitle_char) +
    theme(panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(color = "white"),
          panel.grid.minor = element_line(color = "white"),
          panel.border=element_rect(color="black", fill=NA, size=0.75))+
    theme(axis.text.x = element_text(angle = 55, hjust = 1,vjust=1)) +
    theme(axis.text=element_text(size=axisTextSize_num),
          axis.title.x=element_text(size=axisTitleSize_num),
          axis.title.y=element_text(size=axisTitleSize_num),
          legend.title = element_text(size=legendTitleSize_num),
          legend.text = element_text(size=legendTextSize_num)) +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))+
    labs(fill = legendTitle_char)+
    ylim(ymin_num,ymax_num)+
    theme(legend.key.size = unit(legendSize_num,"line"),
          legend.position = legendPosition_char)+
    theme(plot.margin=unit(c(extraMargin_num,extraMargin_num,extraMargin_num,extraMargin_num),"cm"))
}
