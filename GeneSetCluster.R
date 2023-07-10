require(GeneSetCluster)
library(limma)
library(data.table)
library(ggplot2)

con="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Outputs/PLEIO_Output_with_Concordant_Specific_Magma_Magma_output/PLEIO_Output_with_Concordant_Specific_Magma.gsa.out"
dis="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Outputs/PLEIO_Output_with_Disconcordant_Specific_Magma_Magma_output/PLEIO_Output_with_Disconcordant_Specific_Magma.gsa.out"


con_df<-fread(con)
dis_df<-fread(dis)

dis_df<-dis_df[,c("FULL_NAME","CommonGenes","P_fdr_bh_corr","Total_NumberGenes" ,"NGENES","P")]
con_df<-con_df[,c("FULL_NAME","CommonGenes","P_fdr_bh_corr","Total_NumberGenes" ,"NGENES","P")]

dis_df<-dis_df[dis_df$P<0.05,]
con_df<-con_df[con_df$P<0.05,]

dis_df<-dis_df[dis_df$P_fdr_bh_corr<0.05,]
con_df<-con_df[con_df$P_fdr_bh_corr<0.05,]


dis_df2<-dis_df[,c("FULL_NAME","P_fdr_bh_corr","Total_NumberGenes","NGENES")]
colnames(dis_df2)<-c("FULL_NAME","Discordant_P_fdr_bh_corr","Discordant_Total_NumberGenes","Discordant_NGENES")
con_df2<-con_df[,c("FULL_NAME","P_fdr_bh_corr","Total_NumberGenes","NGENES")]
colnames(con_df2)<-c("FULL_NAME","Concordant_P_fdr_bh_corr","Concordant_Total_NumberGenes","Concordant_NGENES")

con_dis <- merge(dis_df2, con_df2, by = "FULL_NAME", all = TRUE)


#Creating Combine
IPA.KOvsWT.PathwayObject <- ObjectCreator(Pathways = c(dis_df$`FULL_NAME`,
                                                       con_df$`FULL_NAME`), 
                                          Molecules = c(dis_df$CommonGenes,
                                                       con_df$CommonGenes),
                                          Groups = c(rep("Discordant", times = nrow(dis_df)),
                                                     rep("Concordannt", times = nrow(con_df))),
                                          Source = "IPA",
                                          Type = "Canonical_Pathways",#Optional
                                          structure = "SYMBOL",
                                          organism ="org.Hs.eg.db",
                                          sep = ",")


ShowExperimentdata(Object =IPA.KOvsWT.PathwayObject )
ShowMeta(Object =IPA.KOvsWT.PathwayObject)

#Combine gene sets
man.Great.Object2 <- CombineGeneSets(Object = IPA.KOvsWT.PathwayObject)



methods="kmeans, kmeans_group, Hierarchical and Hierarchical_group"
##clust:  "silhouette","elbow", "gap"
#Optimal number of clusters

pdf("output_figure_gap_kmeans.pdf")
OptimalGeneSets(object = man.Great.Object2,method = "gap",max_cluster = 24,
                cluster_method = "kmeans",main = "Kmeans for 24 clusters")

dev.off()


#Cluster Gene sets
man.Great.Object3 <- ClusterGeneSets(Object = man.Great.Object2, clusters = 11, method = "kmeans")



# Create the plot using PlotGeneSets()
plot <- PlotGeneSets(Object = man.Great.Object3,fontsize = 3,legend = TRUE,
                     annotation.mol = FALSE,RR.max = 60,
                     main = "Great_Background clustered with Kmeans \n Disease Ontology and GO Biological Process")

# Save the plot to a file
ggsave("RRmax60_output_Kmeans_group_plot_fdr_11group.png", plot, width = 8, height = 6, dpi = 300)


result_df<-as.data.frame(man.Great.Object3@Data)
result_df<-result_df[,c("Pathways","Groups","Pval","Ratio","MoleculesCount","GeneSet","GeneSets","RR_name","cluster" ,"mean.RR.cl" ,"sum.RR.cl")]


merged_df <- merge(con_dis, result_df, by.x = "FULL_NAME", by.y = "Pathways", all = TRUE)


write.csv(merged_df,"Pathway_geneset_kmeans.csv")

patcorr<-as.data.frame(man.Great.Object3@Data.RR)
write.csv(patcorr,"Pathway_geneset_correlation_kmeans.csv")


man.Great.Object4 <- (Object = man.Great.Object3, breakup.cluster = 4, sub.cluster=2)   
