library(Seurat)
library(GO.db)

### Freshly isolated RPE

load("~/list.RData")
chip86<-obj.list[['chip125286']]
chip86<-subset(chip86, subset = time=="2W")
obj.list[['chip125286']]<-chip86
obj.list<-obj.list[c(1:4,8:10)]

obj.list0<-obj.list[c(5:7)]
obj.list<-obj.list0

DefaultAssay(obj.list[[1]])

for (i in 1:length(obj.list)) {
  DefaultAssay(obj.list[[i]])<-"RNA"
  #obj.list[[i]]@assays$RNA<-as(object = obj.list[[i]]@assays$RNA, Class = "Assay5")
  # obj.list[[i]]<-subset(obj.list[[i]], subset=nFeature_RNA>500)
  # obj.list[[i]]<-subset(obj.list[[i]], subset=nCount_RNA<1500000)
  obj.list[[i]][["percent.mt"]] <- PercentageFeatureSet(obj.list[[i]], pattern = "^MT-")
}

DefaultAssay(obj.list[[1]])


obj<-merge(obj.list[[1]],obj.list[2:3])

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA","percent.mt","donor"))
obj <- RunPCA(obj)

obj <- FindNeighbors(obj, dims = 2:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")

obj <- RunUMAP(obj, dims = 2:30, reduction = "pca", reduction.name = "umap.unintegrated")


obj <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 4:10)
obj <- FindClusters(obj, resolution = 0.5)

obj <- RunUMAP(obj, dims = 4:15, reduction = "integrated.cca")

DimPlot(obj, reduction = "umap", group.by = "donor", label = T)
DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = T)

VlnPlot(obj, c("nCount_RNA", "nFeature_RNA"))


#nk.markers <- FindConservedMarkers(obj)
#head(nk.markers)

markers<-FindAllMarkers(obj, only.pos = T,assay = "RNA", min.pct = "0.1")
markers<-markers[markers$p_val_adj<0.05,]
table(markers$cluster)










### CCA

obj.list <- lapply(X = obj.list, FUN = SCTransform, assay="RNA")

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
obj.list <- lapply(X = obj.list, FUN = RunPCA, features = features)

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "LogNormalize",
                                  anchor.features = features, dims = 2:30, reduction = "cca", k.anchor = 40)
obj.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 2:30)


obj.combined <- RunPCA(obj.combined, verbose = FALSE)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 2:30)

ElbowPlot(obj.combined)

#PCAPlot(obj.combined, c(4:5))

obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 2:30)
obj.combined <- FindClusters(obj.combined, resolution = 0.5)

#DefaultAssay(obj.combined)<-"RNA"
#obj.combined<-JoinLayers(obj.combined)
obj.combined<-NormalizeData(obj.combined, assay = "RNA", normalization.method = "LogNormalize")

DimPlot(obj.combined, label = T, label.size = 5)
DimPlot(obj.combined, label = T, label.size = 5, group.by = "donor")

VlnPlot(obj.combined, c("nCount_RNA", "nFeature_RNA"))

obj.combined<-NormalizeData(obj.combined, assay = "RNA")

markers<-FindAllMarkers(obj.combined, only.pos = T,assay = "RNA", min.pct = "0.1")
markers<-markers[markers$p_val_adj<0.05,]
table(markers$cluster)






obj.combined<-obj
## WNN
DefaultAssay(obj.combined)<-"nADT"
VariableFeatures(obj.combined) <- rownames(obj.combined[["nADT"]])
obj.combined <- SCTransform(obj.combined,assay="nADT")
# obj.combined@assays$SCT<-CreateSCTAssayObject(counts = obj.combined@assays$ADT@counts,
#                                               scale.data =as.matrix(obj.combined@assays$nADT@counts))
# obj.combined@assays$SCT@var.features<-rownames(obj.combined[["nADT"]])
# obj.combined@assays$SCT@key<-"nSCT"

obj.combined<-RunPCA(obj.combined,reduction.name = 'apca')

ElbowPlot(obj.combined,reduction="apca")

#dims.list = list(3:30, c(3:10)), modality.weight.name = "RNA.weight",k.nn=30

obj.combined<- FindMultiModalNeighbors(
  obj.combined, reduction.list = list("pca", "apca"), 
  dims.list = list(2:30, c(2:15)), modality.weight.name = "RNA.weight",k.nn=30
)

obj.combined <- RunUMAP(obj.combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

obj.combined <- FindClusters(obj.combined, graph.name = "wsnn", algorithm = 3, resolution = 0.19, verbose = FALSE)

DimPlot(obj.combined,label=T,reduction="wnn.umap", pt.size = 0.3, label.size = 5)




VlnPlot(obj.combined, c("ITGB3","CSF1R","TNFRSF14","B3GAT1","CD24"),
        assay = "nADT", pt.size = 0, stack = T)+NoLegend()

VlnPlot(obj.combined, c("TFPI2","GNGT1","IGFBP5","WFDC1","CXCL14","ID3"),
        assay = "RNA", pt.size = 0, stack = T)+NoLegend()

VlnPlot(obj.combined, c("RPE65","BEST1","SERPINF1","RLBP1","OTX2","VEGFA"),
        assay = "RNA", pt.size = 0, stack = T)+NoLegend()

VlnPlot(obj.combined, c("nCount_RNA","nFeature_RNA","percent.mt","RPE65"), stack = T)


DefaultAssay(obj.combined)<-"RNA"
obj.combined<-NormalizeData(obj.combined, assay = "RNA")

markers<-FindAllMarkers(obj.combined, only.pos = T,assay = "RNA", min.pct = "0.1")
markers<-markers[markers$p_val_adj<0.05,]
table(markers$cluster)


adt.markers<-FindAllMarkers(obj.combined,only.pos=T, assay = "nADT")
adt.markers<-adt.markers[adt.markers$p_val_adj<0.05,]
table(adt.markers$cluster)



library(dplyr)
library(ggplot2)

mk<-markers
mk.top <- mk %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
mk.top<-unique(mk.top$gene)
mk.top<-c(mk.top,"TNFRSF14","CD24")
# mk.top<-c(mk.top,"CD24","TNFRSF14")
# 
mk.top<-c("OLR1","ENTPD1","CCR6","CD37","PTPRC","NCR1","CD24","B3GAT1","CD58","TNFRSF14","ITGA1",
          "CD55","CCR5","NRP1","CR1","SELL","TFRC","CD9","KIR3DL1","SIGLEC7")
# 
# mk.top<-setdiff(mk.top,c("TTR"))

mk.top<-c("LINC01833","HTRA1","PLIN2","HNRNPH3","SLC39A6","SIX3","RB1CC1")

DotPlot(obj.combined, features = mk.top, assay = "RNA",
        cols = c("lightblue", "deepskyblue4"))+
  ylab("Cluster")+
  theme(axis.text.x = element_text(angle=90,size=10),axis.title.x = element_blank())


markers1 <- markers %>% group_by(cluster) %>% top_n(n =1000, wt = p_val_adj)
table(markers$cluster)
table(markers$cluster)


## GO enrichment
library(screp)
library(GO.db)
library(hypeR)
library(ggplot2)
library(RColorBrewer)

markers1<-markers[which(markers$p_val_adj<0.01),]

eres<-enrich_test(markers_df=markers1)

goterms<-names(Term(GOTERM))
names(goterms)<-Term(GOTERM)
GOBP <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="GO:BP")

go_vis_sc<-GO_visualization(eres$Enriched_df,markers_df=markers1,
                            GOcats=GOBP,goterms=goterms,numcats=15)

go_vis_sc$plot

y<-c("sensory perception of light stimulus","response to extracellular stimulus","homeostatic process",
     "locomotion","response to oxidative stress","pigment metabolic process","cell morphogenesis",
     "ion transport","neurogenesis","tissue development","cell cycle","cell division","cell migration","process utilizing autophagic mechanism",
     "programmed cell death","response to wounding","establishment or maintenance of cell polarity",
     "regulation of cell population proliferation","cell morphogenesis","positive regulation of cell adhesion",
     "vacuolar localization","glial cell migration","cellular response to stress",
     "response to lipid","lipid metabolic process","response to endogenous stimulus",
     "phosphorylation","lipid biosynthetic process","regulation of cell differentiation",
     "positive regulation of locomotion","inflammatory response","rhythmic process",
     "pigment biosynthetic process","tissue migration","positive regulation of protein maturation")
y<-sort(unique(y))

y<-c("homeostatic process",
     "cell cycle","cell division","regulation of cell population proliferation",
     "cell migration","tissue migration","glial cell migration",
     "positive regulation of locomotion","locomotion","positive regulation of cell adhesion",
     "lipid metabolic process","lipid biosynthetic process","response to lipid",
     "inflammatory response","response to endogenous stimulus","response to oxidative stress",
     "cellular response to stress","response to extracellular stimulus","response to wounding",
     "process utilizing autophagic mechanism","pigment biosynthetic process",
     "establishment or maintenance of cell polarity","cell morphogenesis")

y<-c("locomotion","cell migration")

p<-y
a<-go_vis_sc$GO_df
p<-unique(a$GO_Names)
a$fdr<--log10(a$FDR)

p<-a %>% group_by(Clusters) %>% top_n(n = 1, wt = fdr)
p<-p$GO_Names

#p<-c("p<-regulation of organelle organization","")

test<-GO_viz_choose(go_vis_sc,markers_df=markers1,
                    chosen_cats=unique(p),
                    goterms=goterms,
                    species_x="Homo sapiens")
test$plot+ guides(colour = guide_legend(override.aes = list(size=5)))



# upset

obj.combined[["RNA"]]<-as(object = obj.combined[["RNA"]], Class = "Assay")

library(UpSetR)

RPEsignature<-read.csv("~/human_RPE_signature_genes_26517551.csv",as.is=T)
RPEsignature<-RPEsignature$Symbol

probs<-function(object) {
  require(foreach)
  require(Seurat)
  clusters<-object@meta.data$seurat_clusters
  rawdata<-object@assays$RNA@counts
  rownames(rawdata)<-rownames(object)
  cluster.probs<-foreach(i=0:max(as.numeric(levels(clusters))),.combine='cbind') %dopar% {
    clust.mat<-rawdata[,clusters==i]
    apply(clust.mat,1,function(x) length(which(x>0))/length(x))
  }
  
  delta.probs<-foreach(i=1:dim(cluster.probs)[2],.combine='cbind') %do% {
    fin<-foreach(m=1:dim(cluster.probs)[2],.combine='cbind') %dopar% {
      cluster.probs[,i]-cluster.probs[,m]
    }
    fin<-apply(fin,1,sum)
    return(fin)
  }
  c.names<-paste("Cluster",as.numeric(levels(clusters)),sep="_")
  colnames(delta.probs)<-c.names
  colnames(cluster.probs)<-c.names
  fin<-list(cluster.probs,delta.probs)
  return(fin)
}

DefaultAssay(obj.combined)<-"RNA"
test<-probs(obj.combined)

y<-intersect(rownames(obj.combined@assays$RNA),RPEsignature)

z<-setdiff(RPEsignature,y)

nc<-7
length(which((apply(test[[1]][intersect(rownames(test[[1]]),z),0:nc],1,max)>0)))
length(which((apply(test[[1]][intersect(rownames(test[[1]]),RPEsignature),0:nc],1,max)>0)))

id<-which((apply(test[[1]][intersect(rownames(test[[1]]),y),],1,max)>0.2))
m1<-test[[1]][y,]

#m1<-test[[1]][names((apply(test[[1]][intersect(rownames(test[[1]]),RPEsignature),],1,max)>0)),]


m<-foreach(i=1:ncol(m1)) %do% {names(which(m1[,i]>0.3))}

names(m)<-paste("Cluster",levels(obj.combined$seurat_clusters),sep=" ")

d<-foreach(i=1:length(y), .combine = rbind) %do% {
  result <- lapply(m, function(x) x[which(x == y[i])])
  result
}

dd<-fromList(m)

rownames(d)<-y
# 
# dd<-d
# 
# dd[!!sapply(dd, length)]<-1
# dd[!sapply(dd, length)]<-0
# dd<-data.frame(dd)

#rownames(dd)<-y

# m2<-m1
# m2[m2>0.5]<-1
# m2[m2<1]<-0
# m1<-as.data.frame(m1)

upfig<-upset(dd,nsets=20,
             matrix.dot.alpha=1,
             order.by="freq",
             show.numbers=T,matrix.color = "black",line.size = 0, sets.bar.color = "black"
             ,att.color = "black", main.bar.color = "black", shade.alpha = 0.5)

upfig


rpe.exp<-unique(unlist(m))
length(rpe.exp)

### cultured

objc<-subset()

y<-intersect(rownames(obj.combined@assays$RNA),RPEsignature)

z<-setdiff(RPEsignature,y)

nc<-7
length(which((apply(test[[1]][intersect(rownames(test[[1]]),z),0:nc],1,max)>0)))
length(which((apply(test[[1]][intersect(rownames(test[[1]]),RPEsignature),0:nc],1,max)>0)))

id<-which((apply(test[[1]][intersect(rownames(test[[1]]),y),],1,max)>0.3))
m1<-test[[1]][y,]

#m1<-test[[1]][names((apply(test[[1]][intersect(rownames(test[[1]]),RPEsignature),],1,max)>0)),]


m<-foreach(i=1:ncol(m1)) %do% {names(which(m1[,i]>0.3))}

names(m)<-paste("Cluster",levels(obj.combined$seurat_clusters),sep=" ")

d<-foreach(i=1:length(y), .combine = rbind) %do% {
  result <- lapply(m, function(x) x[which(x == y[i])])
  result
}

dd<-fromList(m)

rownames(d)<-y

upfig<-upset(dd,nsets=20,
             matrix.dot.alpha=1,
             order.by="freq",
             show.numbers=T,matrix.color = "black",line.size = 0, sets.bar.color = "black"
             ,att.color = "black", main.bar.color = "black", shade.alpha = 0.5)

upfig


rpe.exp<-unique(unlist(m))







#### heatmap

tr<-markers

mk.top <- tr %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
mk.top <- RPEsignature
mk.top<-intersect(tr$gene,RPEsignature)



DotPlot(obj.combined, features = mk.top, assay = "RNA",
        cols = c("lightblue", "deepskyblue4"))+
  ylab("Cluster")+
  theme(axis.text.x = element_text(angle=90,size=10),axis.title.x = element_blank())


DefaultAssay(obj.combined)<-"RNA"
obj.combined<-NormalizeData(obj.combined, assay = "RNA")
obj.combined<-ScaleData(obj.combined, assay = "RNA")

my_colors <- colorRampPalette(c("blue","white", "red"))(100)
hmap<-DoHeatmap(obj.combined, features = mk.top, assay = "RNA", slot = "data")+ scale_fill_gradientn(colours = my_colors)

my_colors <- colorRampPalette(c("black","lightskyblue", "maroon1"))(100)
pdf("heatmap.pdf", width=5, height=9)
hmap+ scale_fill_gradientn(colours = my_colors)
dev.off()



load("~/p0.Rdata")
g<-c("HTRA1","PLIN2","HNRNPH3","SLC39A6","SIX3","RB1CC1")
DotPlot(obj0, features = g, assay = "RNA",
        cols = c("lightblue", "deepskyblue4"))+
  ylab("Cluster")+
  theme(axis.text.x = element_text(angle=90,size=10),axis.title.x = element_blank())
