# analyze RPE morphology
# developer: Farhad Farjood
# cellshaper function takes the output of imageJ ridge detection on Phalloidin images of RPE monolayers and
# returns the results of morphometric analysis

library(EBImage)
library(RColorBrewer)
library(ggplot2)
library(plyr)
library(doParallel)
library(twosamples)
library(gridExtra)

# function to find neighbor cells
cellneighbor<-function(image.label,id){
  require(doParallel)
  require(foreach)
  
  registerDoParallel(makeCluster(3))
  
  y<-image.label
  # border objects
  edgeid<-unique(c(y[1,], y[,1], y[nrow(y),], y[,ncol(y)]))
  # edgeid<-setdiff(edgeid,0)
  # edges<-which(y %in% edgeid)
  
  z<-foreach (i=as.numeric(id), .packages="raster") %dopar% {
    matrixid<-which(y@.Data==i, arr.ind=TRUE)
    matrixid.n<-matrixid; matrixid.n[,"row"]<-matrixid.n[,"row"]-2
    matrixid.w<-matrixid; matrixid.w[,"col"]<-matrixid.w[,"col"]-2
    matrixid.s<-matrixid; matrixid.s[,"row"]<-matrixid.s[,"row"]+2
    matrixid.e<-matrixid; matrixid.e[,"col"]<-matrixid.e[,"col"]+2
    matrixid.nw<-matrixid; matrixid.nw[,"row"]<-matrixid.nw[,"row"]-2;matrixid.nw[,"col"]<-matrixid.nw[,"col"]-2
    matrixid.sw<-matrixid; matrixid.sw[,"row"]<-matrixid.sw[,"row"]+2;matrixid.sw[,"col"]<-matrixid.sw[,"col"]-2
    matrixid.se<-matrixid; matrixid.se[,"row"]<-matrixid.se[,"row"]+2;matrixid.se[,"col"]<-matrixid.se[,"col"]+2
    matrixid.ne<-matrixid; matrixid.ne[,"row"]<-matrixid.ne[,"row"]-2;matrixid.ne[,"col"]<-matrixid.ne[,"col"]+2
    
    ids<-rbind(matrixid.n,matrixid.w,matrixid.s,matrixid.w,
               matrixid.nw,matrixid.sw,matrixid.se,matrixid.ne)
    
    for (j in 1:nrow(ids)){
      if (ids[j,1]<1){
        ids[j,1]<-1
      } else if (ids[j,1]>nrow(y)){
        ids[j,1]<-nrow(y)
      } else if (ids[j,2]<1){
        ids[j,2]<-1
      } else if (ids[j,2]>ncol(y)){
        ids[j,2]<-ncol(y)
      }
    }
    
    neighbor.id<-unique(y@.Data[ids])
    
    # remove border objects
    neighbor.id<-setdiff(neighbor.id, c(i,0,edgeid))
    neighbor.id
  }
  names(z)<-id
  z
}

# function to perform morphometric analysis of RPE cells
# using fluorescence microscope images of cell membrane using the EBImage package
cellshaper<-function(image,brush_size=2,feature="m.eccentricity",nbin=4,refcut=NULL,colpal="PuBu",ref){
  
  require(RColorBrewer)
  require(EBImage)
  
  # load again and label
  x<-readImage(image)
  #display(x, title='Binary')
  x = opening(x, makeBrush(brush_size, shape='disc'))
  y <- bwlabel(x)
  #display(normalize(y), title='Segmented')
  z<-colorLabels(y, normalize = T)
  #display(z, title='Colored segmentation')
  
  ftrs<-cbind(computeFeatures.moment(y,x),computeFeatures.shape(y,x))
  
  # change labels
  prop<-as.data.frame(ftrs)
  prop$compactness<-(4*pi*prop$s.area)/(prop$s.perimeter^2)
  
  pix<-0.493^2
  prop$area<-prop$s.area*pix
  
  prop$minoraxis<-sqrt((1-prop$m.eccentricity^2)*prop$m.majoraxis^2)
  prop$aspectratio<-prop$m.majoraxis/prop$minoraxis
  prop$elongation<-prop$minoraxis/prop$m.majoraxis
  
  # only use for reference 2w
  if (ref){
    refcut<-quantile(prop[,feature], seq(0,1,1/nbin))
    if (max(refcut)<1){
      refcut[1]<-0
      refcut[nbin+1]<-1
    }
  }
  q<-quantile(prop[,feature], seq(0,1,1/nbin))
  if (max(q)>max(refcut)){refcut[nbin+1]<-max(q)}
  # use 2-week cutoffs for binning
  
  circbin<-prop
  circbin$quantile<-cut(prop[,feature], refcut, labels<-c(1:nbin))
  circbin$id<-rownames(circbin)
  
  # #binning with size
  # dbin<-prop%>% mutate(quantile = ntile(s.area,4))%>%group_by(quantile)
  # dbin$id<-rownames(dbin)
  # 
  # # binning with morphology
  # circbin<-prop%>% mutate(quantile = ntile(circ,4))%>%group_by(quantile)
  # circbin$id<-rownames(circbin)
  
  #create color palette
  pal <- heat.colors(n = nbin+1)
  pal<-brewer.pal(n=nbin+1,colpal)
  palrgb <- col2rgb(pal)/255
  
  z2<-z@.Data
  l<-levels(as.factor(circbin$quantile))
  
  for (i in 1:length(l)) {
    mx<-which(circbin$quantile==as.numeric(l[i]))
    mx<-which(y %in% mx)
    z2[,,1][mx]<-palrgb[1,i]
    z2[,,2][mx]<-palrgb[2,i]
    z2[,,3][mx]<-palrgb[3,i]
  }
  
  # remove border objects
  edgeid<-unique(c(y[1,], y[,1], y[nrow(y),], y[,ncol(y)]))
  edges<-which(y %in% edgeid)
  z2[,,1][edges]<-0
  z2[,,2][edges]<-0
  z2[,,3][edges]<-0
  circbin<-circbin[setdiff(1:nrow(circbin),edgeid),]
  
  z@.Data<-z2
  display(z, title='Colored segmentation')
  
  nbr<-cellneighbor(y,circbin$id)
  circbin$neighbor<-unlist(lapply(nbr,function(x) length(x)))
  l<-list(bins=refcut,image=z,labels=y,features=circbin, neighborhood=nbr)
  return(l)
}

#function to label a cell in an Image class image
labelcells<-function(image,image.label,cellid){
  z2<-image@.Data
  mx<-which(image.label %in% cellid)
  z2[,,1][mx]<-1
  z2[,,2][mx]<-0
  z2[,,3][mx]<-0
  
  image@.Data<-z2
  display(image, title='Colored segmentation')
  return(image)
}

image.intensity<-function(img,mask){
  
  labels<-mask$labels
  intensities <- NULL
  for (i in 1:max(labels)) {
    # Extract pixels belonging to each labeled object
    obj_pixels <- img[labels == i]
    # Calculate intensity (mean, sum, median, etc.) of the object
    intensity <- mean(obj_pixels)  # You can replace 'mean' with other functions like 'sum', 'median', etc.
    # Store the intensity value
    intensities <- c(intensities, intensity)
  }
  
  names(intensities)<-1:length(intensities)
  
  img<-labels
  for (i in 1:max(labels)){
    img[img==i]<-intensities[i]
  }
  
  mask$features$intensity<-intensities[mask$features$id]
  mask$features$norm.intensity<-intensities[mask$features$id]/max(intensities)
  return(list(mask,img))
}

# data analysis example

refcut1 = c(0,150,250,350,1e10)
nbin1<-length(refcut1)-1



#### 57 / 24

R1.57.24<-cellshaper(image = "~/Flatmount/RPE1_ActinGreen_CD24_CD57_3/final.tif",brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")
R1.57<-readImage("~/Flatmount/RPE1_ActinGreen_CD24_CD57_3/cd57.jpg")
R1.24<-readImage("~/Flatmount/RPE1_ActinGreen_CD24_CD57_3/cd24.jpg")

display(R1.57.24[[2]])


pixel.scale<-0.512

af24<-image.intensity(R1.24,R1.57.24)
af57<-image.intensity(R1.57,R1.57.24)

display(af24[[2]])

plot(af24[[1]]$features$intensity,af24[[1]]$features$area, pch = 20, cex=0.5)
plot(density(af24[[1]]$features$intensity))


ids<-af24[[1]]$features$id[which(af24[[1]]$features$intensity>0.05)]
k<-labelcells(R1.57.24[[2]],image.label =R1.57.24[[3]] ,cellid =ids)
display(k)


#ECDF
set.seed(314159)
vec1=af24[[1]]$features$area[which(af24[[1]]$features$intensity>0.05)]
vec2=af24[[1]]$features$area[which(af24[[1]]$features$intensity<0.05)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df24<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df24$cd<-df24$cd*0.512

ec24<-ggplot(df24, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
ec24


p24 <- ggplot(df24, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("CD24")
p24

h241 <- ggplot(df24, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
h241


#57

plot(density(af57[[1]]$features$intensity))

ids<-af57[[1]]$features$id[which(af57[[1]]$features$intensity>0.08)]
k<-labelcells(R1.57.24[[2]],image.label =R1.57.24[[3]] ,cellid =ids)
display(k)


set.seed(314159)
vec1=af57[[1]]$features$area[which(af57[[1]]$features$intensity>0.08)]
vec2=af57[[1]]$features$area[which(af57[[1]]$features$intensity<0.08)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df57<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df57$cd<-df57$cd*0.512

ec57<-ggplot(df57, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()



p57 <- ggplot(df57, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("CD57")
p57


h571 <- ggplot(df57, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
h571



##------------- P1 CD57 ---------------

R2.57m<-cellshaper(image = "~/Flatmount/CD57_270_macula/mask.tif",brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")
R2.57<-readImage("~/Flatmount/CD57_270_macula/CD57_270_macula-Image Calculator-27-Add Channels-28_c1.jpg")

display(R2.57m[[2]])

af57p1<-image.intensity(R2.57,R2.57m)


plot(af57p1[[1]]$features$intensity,af57p1[[1]]$features$area, pch = 20, cex=0.5)
plot(density(af57p1[[1]]$features$intensity))


ids<-af57p1[[1]]$features$id[which(af57p1[[1]]$features$intensity>0.05)]
k<-labelcells(R2.57m[[2]],image.label =R2.57m[[3]] ,cellid =ids)
display(k)

set.seed(314159)
vec1=af57p1[[1]]$features$area[which(af57p1[[1]]$features$intensity>0.05)]
vec2=af57p1[[1]]$features$area[which(af57p1[[1]]$features$intensity<0.05)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df57p1<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df57p1$cd<-df57p1$cd*0.512

#ECDF
ec57p1<-ggplot(df57p1, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
ec57p1



p57p2 <- ggplot(df57p1, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+  ylab("CD57p1")
p57p2



h572 <- ggplot(df57p1, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
h572

##### end P1 CD57

### start CD270

R1.270<-cellshaper(image = "~/Flatmount/RPE1_CD57_270_10/ridge-pruned.tif",brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")
R1.27<-readImage("~/Flatmount/RPE1_CD57_270_10/RPE1_CD57_270_10_c3.jpg")

a<-readImage("~/Flatmount/RPE1_CD57_270_10/ridge-pruned.tif")
  
af270<-image.intensity(R1.27,R1.270)


plot(af270[[1]]$features$intensity,af270[[1]]$features$area, pch = 20, cex=0.5)
plot(density(af270[[1]]$features$intensity))


ids<-af270[[1]]$features$id[which(af270[[1]]$features$intensity>0.06)]
k<-labelcells(R1.270[[2]],image.label =R1.270[[3]] ,cellid =ids)
display(k)

set.seed(314159)
vec1=af270[[1]]$features$area[which(af270[[1]]$features$intensity>0.06)]
vec2=af270[[1]]$features$area[which(af270[[1]]$features$intensity<0.06)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df270<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df270$cd<-df270$cd*0.512

#ECDF
ec270<-ggplot(df270, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
ec270



p270 <- ggplot(df270, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+  ylab("CD270")
p270


h270 <- ggplot(df270, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
h270

#------ end ---------

#------ CD24 P3 -----

RPE361.CD24<-cellshaper(image = "~/calculated/CD24_7-Image Calculator-24-Add Channels-26/ridge-pruned.tif",
                        brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")

RRPE361.24<-readImage("~/calculated/CD24_7-Image Calculator-24-Add Channels-26/CD24_7-Image Calculator-24-Add Channels-26_c1.jpg")

display(RPE361.CD24[[2]])
display(RRPE361.24)

pixel.scale<-0.512



a36124<-image.intensity(RRPE361.24,RPE361.CD24)


plot(a36124[[1]]$features$intensity,a36124[[1]]$features$area, pch = 20, cex=0.5)
plot(density(a36124[[1]]$features$intensity))

ids<-a36124[[1]]$features$id[which(a36124[[1]]$features$intensity>0.02)]
k<-labelcells(RPE361.CD24[[2]],image.label =RPE361.CD24[[3]] ,cellid =ids)
display(k)

k<-labelcells(RPE361.CD24[[2]],image.label =RPE361.CD24[[3]] ,cellid =338)
display(k)


#ECDF
set.seed(314159)
vec1=a36124[[1]]$features$area[which(a36124[[1]]$features$intensity>0.015)]
vec2=a36124[[1]]$features$area[which(a36124[[1]]$features$intensity<0.015)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df242<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df242$cd<-df242$cd*0.512

ec24.2<-ggplot(df242, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
ec24.2


ap24 <- ggplot(df242, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+  ylab("24.2")
ap24


h242 <- ggplot(df242, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
h242





# CD61

RPE.CD61<-cellshaper(image = "~/calculated/CD61_human_flatmount-Image Calculator-12-Add Channels-21/ridge-pruned.tif",
                        brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")

RRPE.61<-readImage("~/calculated/CD61_human_flatmount-Image Calculator-12-Add Channels-21/CD61.jpg")
display(RRPE.61)
display(RPE.CD61$image)

pixel.scale<-1.024



a61<-image.intensity(RRPE.61,RPE.CD61)


plot(a61[[1]]$features$intensity,a61[[1]]$features$area, pch = 20, cex=0.5)
plot(density(a61[[1]]$features$intensity))


ids<-a61[[1]]$features$id[which(a61[[1]]$features$intensity>0.04)]
k<-labelcells(RPE.CD61[[2]],image.label =RPE.CD61[[3]] ,cellid =ids)
display(k)

#ECDF
set.seed(314159)
vec1=a61[[1]]$features$area[which(a61[[1]]$features$intensity>0.04)]
vec2=a61[[1]]$features$area[which(a61[[1]]$features$intensity<0.04)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df61<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df61$cd<-df61$cd*0.512

ec61<-ggplot(df61, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
ec61


p61 <- ggplot(df61, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("CD61")
p61


h61 <- ggplot(df61, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
h61



# flatmount combined

p57p2+ylim(c(0,400))+p57+ylim(c(0,400))+ap24+ylim(c(0,400))+p24+ylim(c(0,400))+p270+ylim(c(0,400))+p61+ylim(c(0,400))

h571+xlim(c(0,400))+h572+xlim(c(0,400))+h241+xlim(c(0,400))+h242+xlim(c(0,400))+h270+xlim(c(0,400))+h61+xlim(c(0,400))



df270$ab<-df270$meta; df270$ab[which(df270$ab=="pos")]<-"CD270"
df57$ab<-df57$meta; df57$ab[which(df57$ab=="pos")]<-"CD57"
df57p1$ab<-df57p1$meta; df57p1$ab[which(df57p1$ab=="pos")]<-"CD57"
df24$ab<-df24$meta; df24$ab[which(df24$ab=="pos")]<-"CD24"
df242$ab<-df242$meta; df242$ab[which(df242$ab=="pos")]<-"CD24"
df61$ab<-df61$meta; df61$ab[which(df61$ab=="pos")]<-"CD61"


df<-rbind(df270,df57,df24,df242,df61)
df<-df[which(df$ab!="neg"),]

p <- ggplot(df, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
p


df<-df57p1
df$ab[which(df$ab=="neg")]<-"neg1"
#df$ab[which(df$ab=="CD57")]<-"CD57-1"

df<-rbind(df,df57)

p1 <- ggplot(df, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,400))+
  theme_classic()
p1

# CD24


df<-df24
df$ab[which(df$ab=="neg")]<-"neg1"
#df$ab[which(df$ab=="CD24")]<-"CD24-1"

df<-rbind(df,df242)

p <- ggplot(df, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,400))+
  theme_classic()
p

p1+p

aggregate(cd ~ ab, data = df57, FUN = mean)
aggregate(cd ~ ab, data = df57p1, FUN = mean)

aggregate(cd ~ ab, data = df242, FUN = mean)
aggregate(cd ~ ab, data = df242, FUN = sd)

aggregate(cd ~ ab, data = df24, FUN = mean)
aggregate(cd ~ ab, data = df24, FUN = sd)

refcut1 = c(0,150,250,350,1e10)
nbin1<-length(refcut1)-1


#### cultured

RPE353.57<-cellshaper(image = "~/Cultured/RPE353_57_270-3/Actin_ridge-pruned.tif",
                      brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")
display(RPE353.57[[2]])

img<-R3.270<-readImage("~/Cultured/RPE353_57_270-3/RPE353_57_270-3_c3.jpg")
display(img)


c270<-image.intensity(R3.270,RPE353.57)
temp<-c270

plot(temp[[1]]$features$intensity,temp[[1]]$features$area, pch = 20, cex=0.5)
plot(density(temp[[1]]$features$intensity))


ids<-temp[[1]]$features$id[which(temp[[1]]$features$intensity>0.05)]
k<-labelcells(RPE353.57[[2]],image.label =RPE353.57[[3]] ,cellid =ids)
display(k)

#ECDF
set.seed(314159)
vec1=c270[[1]]$features$area[which(c270[[1]]$features$intensity>0.05)]
vec2=c270[[1]]$features$area[which(c270[[1]]$features$intensity<0.05)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df3270<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df3270$cd<-df3270$cd*0.512

cs270<-ggplot(df3270, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
cs270


pc270 <- ggplot(df3270, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("c270")
pc270

h53.270 <- ggplot(df3270, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(0,500))+
  theme_classic()
h53.270

#cultured 57


img2<-R3.57<-readImage("~/Cultured/RPE353_57_270-3/RPE353_57_270-3-Shading.jpg")
display(img2)



c57<-image.intensity(R3.57,RPE353.57)
temp<-c57

plot(temp[[1]]$features$intensity,temp[[1]]$features$area, pch = 20, cex=0.5)
plot(density(temp[[1]]$features$intensity))


ids<-temp[[1]]$features$id[which(temp[[1]]$features$intensity>0.12)]
k<-labelcells(RPE353.57[[2]],image.label =RPE353.57[[3]] ,cellid =ids)
display(k)

#ECDF
set.seed(314159)
vec1=temp[[1]]$features$area[which(temp[[1]]$features$intensity>0.12)]
vec2=temp[[1]]$features$area[which(temp[[1]]$features$intensity<0.12)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df357<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df357$cd<-df357$cd*0.512

cs57<-ggplot(df357, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
cs57


pc57 <- ggplot(df357, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("c57")
pc57

p357 <- ggplot(df357, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(0,500))+
  theme_classic()
p357


# RPE354 CD24

RPE354.24<-cellshaper(image = "~/Cultured/RPE354_24/mask.tif",
                      brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")
display(RPE354.24[[2]])

img<-R4.24<-readImage("~/Cultured/RPE354_24/RPE354_24_c3.jpg")
display(img)


c424<-image.intensity(R4.24,RPE354.24)

temp<-c424
mask<-RPE354.24

plot(temp[[1]]$features$intensity,temp[[1]]$features$area, pch = 20, cex=0.5)
plot(density(temp[[1]]$features$intensity))


ids<-temp[[1]]$features$id[which(temp[[1]]$features$intensity>0.04)]
k<-labelcells(mask[[2]],image.label =mask[[3]] ,cellid =ids)
display(k)

#ECDF
set.seed(314159)
vec1=temp[[1]]$features$area[which(temp[[1]]$features$intensity>0.04)]
vec2=temp[[1]]$features$area[which(temp[[1]]$features$intensity<0.04)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df424<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df424$cd<-df424$cd*0.512

cs424<-ggplot(df424, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
cs424


pc424 <- ggplot(df424, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("c24")
pc424

p <- ggplot(df424, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(0,500))+
  theme_classic()
p


# RPE355 CD24

RPE355.24<-cellshaper(image = "~/Cultured/RPE355_24-3/label-pruned.tif",
                      brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")
display(RPE355.24[[2]])

img<-R5.24<-readImage("~/Cultured/RPE355_24-3/RPE355_24-3_c3.jpg")
display(img)


c24<-image.intensity(R5.24,RPE355.24)

temp<-c24
mask<-RPE355.24

plot(temp[[1]]$features$intensity,temp[[1]]$features$area, pch = 20, cex=0.5)
plot(density(temp[[1]]$features$intensity))


ids<-temp[[1]]$features$id[which(temp[[1]]$features$intensity>0.04)]
k<-labelcells(mask[[2]],image.label =mask[[3]] ,cellid =ids)
display(k)

#ECDF
set.seed(314159)
vec1=temp[[1]]$features$area[which(temp[[1]]$features$intensity>0.05)]
vec2=temp[[1]]$features$area[which(temp[[1]]$features$intensity<0.05)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df524<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df524$cd<-df524$cd*0.512

cs24<-ggplot(df524, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
cs24


pc24 <- ggplot(df524, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("c24")
pc24

p <- ggplot(df524, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(0,500))+
  theme_classic()
p



# RPE355 CD270/CD57

RPE355.57<-cellshaper(image = "~/Cultured/RPE355_57_270_20x_21/mask.tif",
                      brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")
display(RPE355.57[[2]])

img<-R5.270<-readImage("~/Cultured/RPE355_57_270_20x_21/RPE355_57_270_20x_2_c3.jpg")
display(img)

c5270<-image.intensity(R5.270,RPE355.57)

temp<-c5270
mask<-RPE355.57

plot(temp[[1]]$features$intensity,temp[[1]]$features$area, pch = 20, cex=0.5)
plot(density(temp[[1]]$features$intensity))


ids<-temp[[1]]$features$id[which(temp[[1]]$features$intensity>0.03)]
k<-labelcells(mask[[2]],image.label =mask[[3]] ,cellid =ids)
display(k)

#ECDF
set.seed(314159)
vec1=c5270[[1]]$features$area[which(c5270[[1]]$features$intensity>0.03)]
vec2=c5270[[1]]$features$area[which(c5270[[1]]$features$intensity<0.03)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df5270<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df5270$cd<-df5270$cd*0.512

cs5270<-ggplot(df5270, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
cs5270


pc5270 <- ggplot(df5270, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("c5270")
pc5270

p5270 <- ggplot(df5270, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(0,500))+
  theme_classic()
p5270

#cultured 355 CD57

img2<-R5.57<-readImage("~/Cultured/RPE355_57_270_20x_21/RPE355_57_270_20x_2_c4.jpg")
display(img2)

c557<-image.intensity(R5.57,RPE355.57)
temp<-c557
mask<-RPE355.57

plot(temp[[1]]$features$intensity,temp[[1]]$features$area, pch = 20, cex=0.5)
plot(density(temp[[1]]$features$intensity))


ids<-temp[[1]]$features$id[which(temp[[1]]$features$intensity>0.06)]
k<-labelcells(mask[[2]],image.label =mask[[3]] ,cellid =ids)
display(k)

#ECDF
set.seed(314159)
vec1=temp[[1]]$features$area[which(temp[[1]]$features$intensity>0.06)]
vec2=temp[[1]]$features$area[which(temp[[1]]$features$intensity<0.06)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df557<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df557$cd<-df557$cd*0.512

cs557<-ggplot(df557, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
cs557


pc557 <- ggplot(df557, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("c57")
pc557

p557 <- ggplot(df557, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(0,500))+
  theme_classic()
p557



## CD61 cultured RPE358

RPE358.61<-cellshaper(image = "~/Cultured/raw_figures/ITGB3_mask.tif",
                      brush_size = 2, nbin = nbin1,ref = F, refcut = refcut1, colpal = "PuBu", feature = "area")
display(RPE358.61[[2]])

img<-R8.61<-readImage("~/Cultured/raw_figures/1_1-Create Image Subset-13/c2.jpg")
display(img)


c61<-image.intensity(R8.61,RPE358.61)

temp<-c61
mask<-RPE358.61

plot(temp[[1]]$features$intensity,temp[[1]]$features$area, pch = 20, cex=0.5)
plot(density(temp[[1]]$features$intensity))


ids<-temp[[1]]$features$id[which(temp[[1]]$features$intensity>0.07)]
k<-labelcells(mask[[2]],image.label =mask[[3]] ,cellid =ids)
display(k)

#ECDF
set.seed(314159)
vec1=temp[[1]]$features$area[which(temp[[1]]$features$intensity>0.07)]
vec2=temp[[1]]$features$area[which(temp[[1]]$features$intensity<0.07)]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df861<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("neg", length(vec2))))
df861$cd<-df861$cd*0.512

cs61<-ggplot(df861, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
cs61


pc61 <- ggplot(df861, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,750))+ylab("c61")
pc61



p <- ggplot(df861, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(0,500))+
  theme_classic()
p



# combined
pc270+pc57
pc24+pc61



df3270$ab<-df3270$meta; df3270$ab[which(df3270$ab=="pos")]<-"CD270"
df5270$ab<-df5270$meta; df5270$ab[which(df5270$ab=="pos")]<-"CD270"

df357$ab<-df357$meta; df357$ab[which(df357$ab=="pos")]<-"CD57"
df557$ab<-df557$meta; df557$ab[which(df557$ab=="pos")]<-"CD57"

df424$ab<-df424$meta; df424$ab[which(df424$ab=="pos")]<-"CD24"
df524$ab<-df524$meta; df524$ab[which(df524$ab=="pos")]<-"CD24"
df861$ab<-df861$meta; df861$ab[which(df861$ab=="pos")]<-"CD61"




df<-rbind(df3270,df5270,df357,df557,df424,df524,df861)
df<-df[which(df$ab!="neg"),]

p <- ggplot(df, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
p




data61<-df861
data270<-rbind(df3270,df5270)
data57<-rbind(df357,df557)
data24<-rbind(df424,df524)

cul61 <- ggplot(data61, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
cul270 <- ggplot(data270, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
cul24 <- ggplot(data24, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
cul57 <- ggplot(data57, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()

cul57+cul24+cul270+cul61

grid.arrange(cul57,cul24,cul270,cul61, ncol=4)

table(data61$ab)/nrow(data61)
table(data270$ab)/nrow(data270)
table(data57$ab)/nrow(data57)
table(data24$ab)/nrow(data24)



temp<-data270
l<-levels(as.factor(temp$ab))

set.seed(314159)
vec1=temp[which(temp$ab==l[1]),"cd"]
vec2=temp[which(temp$ab=="neg"),"cd"]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)


c<-ggplot(temp, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
c

p <- ggplot(temp, aes(x=log(cd), fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(-2,10))+
  theme_classic()
p

ggplot(temp, aes(x=cd, fill=meta)) +
  geom_density(alpha=0.4)+xlim(c(0,500))+
  theme_classic()


df<-rbind(data270,data57,data24,data61)
df<-df[which(df$ab!="neg"),]

p <- ggplot(df, aes(x=cd, fill=ab)) +
  geom_density(alpha=0.5)+xlim(c(0,500))+
  theme_classic()
p


#### compare cultured vs native


df270
df57
df57p1
df24
df242
df61


data61
data270
data57
data24


# CD61 comparison

nc61<-df61[which(df61$meta=="pos"),]
nc61$meta<-"native"
nc61<-rbind(nc61,data61[which(data61$meta=="pos"),])


#ECDF
set.seed(314159)
vec1=nc61$cd[which(nc61$meta=="native")]
vec2=nc61$cd[which(nc61$meta=="pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)


n61<-ggplot(nc61, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
n61


box61 <- ggplot(nc61, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,400))+ylab("CD61")
box61


# CD24 comparison

nc24<-rbind(df24,df242)
nc24<-nc24[which(nc24$meta=="pos"),]
nc24$meta<-"native"
nc24<-rbind(nc24,data24[sample(which(data24$meta=="pos"),size = nrow(nc24), replace = F),])



#ECDF
set.seed(314159)
vec1=nc24$cd[which(nc24$meta=="native")]
vec2=nc24$cd[which(nc24$meta=="pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)


n24<-ggplot(nc24, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
n24


box24 <- ggplot(nc24, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,400))+ylab("CD24")
box24


# CD270 comparison

nc270<-df270[which(df270$meta=="pos"),]
nc270$meta<-"native"
nc270<-rbind(nc270,data270[sample(which(data270$meta=="pos"), size = nrow(nc270), replace = F),])

table(nc270$meta)

#ECDF
set.seed(314159)
vec1=nc270$cd[which(nc270$meta=="native")]
vec2=nc270$cd[which(nc270$meta=="pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)


n270<-ggplot(nc270, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
n270


box270 <- ggplot(nc270, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,400))+ylab("CD270")
box270

# CD57 comparison

nc57<-df57[which(df57$meta=="pos"),]
nc57$meta<-"native"
nc57<-rbind(nc57,data57[sample(which(data57$meta=="pos"), size = nrow(nc57), replace = F),])

table(nc57$meta)

#ECDF
set.seed(314159)
vec1=nc57$cd[which(nc57$meta=="native")]
vec2=nc57$cd[which(nc57$meta=="pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)


n57<-ggplot(nc57, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
n57


box57 <- ggplot(nc57, aes(x=meta, y=cd)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")+ ylim(c(0,400))+ylab("CD57")
box57


grid.arrange(box57,box270,box24,box61, ncol=4)







###
w2c<-100*table(ddf[ddf$time=="2w",])/nrow(ddf[ddf$time=="2w",])
w4c<-100*table(ddf[ddf$time=="4w",])/nrow(ddf[ddf$time=="4w",])
w8c<-100*table(ddf[ddf$time=="8w",])/nrow(ddf[ddf$time=="8w",])

ddfc<-reshape2::melt(rbind(w2c,w4c,w8c))
colnames(ddfc)<-c("time","quantile","percentage")

ggplot(ddfc, aes(fill=as.factor(quantile), y=percentage, x=time)) +
  geom_bar(position="stack", stat="identity")+scale_fill_brewer(palette="PuBu")+
  theme_classic()


specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)


# create a ggplot graph

c28<-cutrf3182[[4]]
c28$time<-"2w"
c48<-cutrf3184[[4]]
c48$time<-"4w"
c88<-cutrf3188[[4]]
c88$time<-"8w"

c29<-cutrf3192[[4]]
c29$time<-"2w"
c49<-cutrf3194[[4]]
c49$time<-"4w"
c89<-cutrf3198[[4]]
c89$time<-"8w"


c22<-cutrf3222[[4]]
c22$time<-"2w"
c42<-cutrf3224[[4]]
c42$time<-"4w"
c82<-cutrf3228[[4]]
c82$time<-"8w"


# ggplot figures

df<-rbind(c28,c48,c88,c29,c49,c89,c22,c42,c82)

p2 <- ggplot(df, aes(x=quantile, y=neighbor, fill=quantile)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")
p2


p <- ggplot(df, aes(x=quantile, y=neighbor, fill=time))+
  stat_boxplot(geom = "errorbar",width = 0.15)+
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, outlier.alpha = 0.1)+
  theme_classic()+
  scale_fill_brewer(palette="PuBu")
p


p <- ggplot(df, aes(quantile, fill=quantile)) +
  geom_density(alpha=0.4, xlim = c(-0.8, 4.8))+
  theme_classic()
p


mu <- ddply(df, "time", summarise, grp.mean=median(compactness))
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")
