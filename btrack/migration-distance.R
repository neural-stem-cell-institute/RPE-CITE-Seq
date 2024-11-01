# calculate migration distance from btrack output

library(foreach)
#library(rhdf5)
library(ggplot2)
library(ggridges)
library(RColorBrewer)
library(reshape2)
library(gridExtra)


path_dist <- function(x, y) {
  d<-c()
  for (i in 2:length(x)){
    d[i-1]<-sqrt((x[i]-x[i-1])^2+(y[i]-y[i-1])^2)
  }
  d<-sum(d)
  dt<-sqrt((x[length(x)]-x[1])^2+(y[length(x)]-y[1])^2)
  
  return(c(d,dt))
}

ref.dist<-function(x, y,ref){
  d<-c()
  for (i in 1:length(x)){
    d[i]<-sqrt((x[i]-ref[1])^2+(y[i]-ref[2])^2)
  }
  return(d)
}


neg<-"~/Desktop/R/napari/115-tracking/tracks_318-115neg.csv"
pos<-"~/Desktop/R/napari/115-tracking/tracks_318-115pos.csv"

btrackv1<- function(neg,pos, line, dist.factor=NULL) {
  
  track<-read.csv(neg, header=F)
  track_pos<-read.csv(pos, header=F)
  
  colnames(track)<-c("all","ID", "t", "x", "y", "root", "parent", "children", "generation", "fate")
  colnames(track_pos)<-c("all","ID", "t", "x", "y", "root", "parent", "children", "generation", "fate")
  
  
  track<-foreach(i=1:nrow(track)) %do% {
    
    t<-as.numeric(unlist(strsplit(track$t[i], split=",|\\[|\\]| ")))
    t<-t[!is.na(t)]
    
    x<-as.numeric(unlist(strsplit(track$x[i], split=",|\\[|\\]| ")))
    x<-x[!is.na(x)]
    
    y<-as.numeric(unlist(strsplit(track$y[i], split=",|\\[|\\]| ")))
    y<-y[!is.na(y)]
    
    z<-rep(0, length(x))
    
    label<-rep(i, length(x))
    
    p<-rep(track$parent[i], length(x))
    
    r<-rep(track$root[i], length(x))
    
    g<-rep(track$generation[i], length(x))
    
    tr<-unlist(strsplit(track$all[i], split="[(]"))
    a<-as.numeric(unlist(strsplit(tr[14], split=",|\\[|\\]| ")))
    a<-a[!is.na(a)]
    
    f<-rep(track$fate[i], length(x))
    
    d<-data.frame(t=t,x=x,y=y,z=z,
                  label=label,
                  area=a,
                  parent=p,
                  root=r,
                  generation=g,
                  fate=f)
  }
  track_pos<-foreach(i=1:nrow(track_pos)) %do% {
    
    t<-as.numeric(unlist(strsplit(track_pos$t[i], split=",|\\[|\\]| ")))
    t<-t[!is.na(t)]
    
    x<-as.numeric(unlist(strsplit(track_pos$x[i], split=",|\\[|\\]| ")))
    x<-x[!is.na(x)]
    
    y<-as.numeric(unlist(strsplit(track_pos$y[i], split=",|\\[|\\]| ")))
    y<-y[!is.na(y)]
    
    z<-rep(0, length(x))
    
    label<-rep(i, length(x))
    
    p<-rep(track_pos$parent[i], length(x))
    
    r<-rep(track_pos$root[i], length(x))
    
    g<-rep(track_pos$generation[i], length(x))
    
    tr<-tr<-unlist(strsplit(track_pos$all[i], split="[(]"))
    a<-as.numeric(unlist(strsplit(tr[14], split=",|\\[|\\]| ")))
    a<-a[!is.na(a)]
    
    f<-rep(track_pos$fate[i], length(x))
    
    d<-data.frame(t=t,x=x,y=y,z=z,
                  label=label,
                  area=a,
                  parent=p,
                  root=r,
                  generation=g,
                  fate=f)
  }
  
  a<-unlist(lapply(1:length(track), function(x) track[[x]]$area[1]))
  l<-unlist(lapply(1:length(track), function(x) nrow(track[[x]])))
  track<-track[a>10&l>1]
  
  a<-unlist(lapply(1:length(track_pos), function(x) track_pos[[x]]$area[1]))
  l<-unlist(lapply(1:length(track_pos), function(x) nrow(track_pos[[x]])))
  track_pos<-track_pos[a>10&l>1]
  
  track.n<-track
  
  for (i in 1:length(track.n)){
    track.n[[i]]$label<-rep(track.n[[i]]$label[1],nrow(track.n[[i]]))
  }
  
  
  np.track<-foreach(i=1:length(track.n),.combine=rbind) %do% {
    track.n[[i]]
  }
  
  np.track<-np.track[,c("label","z","y","x")]
  
  
  #calculate distance
  track.dist<-foreach(i=1:length(track),.combine=rbind) %do% {
    path_dist(track[[i]]$x,track[[i]]$y)}
  
  track.dist<-data.frame(track.dist)
  colnames(track.dist)<-c("path.distance","overall.distance")
  
  track.dist<-track.dist*dist.factor
  
  track.dist.pos<-foreach(i=1:length(track_pos),.combine=rbind) %do% {
    path_dist(track_pos[[i]]$x,track_pos[[i]]$y)}
  
  track.dist.pos<-data.frame(track.dist.pos)
  colnames(track.dist.pos)<-c("path.distance","overall.distance")
  track.dist.pos<-track.dist.pos*dist.factor
  
  
  dim(track.dist)
  dim(track.dist.pos)
  
  n<-min(nrow(track.dist),nrow(track.dist.pos))
  
  track.dist<-track.dist[sample(rownames(track.dist), n),]
  track.dist.pos<-track.dist.pos[sample(rownames(track.dist.pos), n),]
  
  plot(density(track.dist$path.distance), col="blue")
  lines(density(track.dist.pos$path.distance), col="red")
  
  plot(density(track.dist$overall.distance), col="blue", xlim=c(0,300))
  lines(density(track.dist.pos$overall.distance), col="red")
  
  d<-data.frame(neg=track.dist$path.distance, pos=track.dist.pos$path.distance, line=line)
  return(d)
  
}

btrack<- function(neg,pos, line, dist.factor=NULL) {
  
  track<-read.csv(neg, header=F)
  track_pos<-read.csv(pos, header=F)
  
  colnames(track)<-c("all","ID", "t", "x", "y", "root", "parent", "children", "generation", "fate")
  colnames(track_pos)<-c("all","ID", "t", "x", "y", "root", "parent", "children", "generation", "fate")
  
  
  track<-foreach(i=1:nrow(track)) %do% {
    
    t<-as.numeric(unlist(strsplit(track$t[i], split=",|\\[|\\]| ")))
    t<-t[!is.na(t)]
    
    x<-as.numeric(unlist(strsplit(track$x[i], split=",|\\[|\\]| ")))
    x<-x[!is.na(x)]
    
    y<-as.numeric(unlist(strsplit(track$y[i], split=",|\\[|\\]| ")))
    y<-y[!is.na(y)]
    
    z<-rep(0, length(x))
    
    label<-rep(i, length(x))
    
    p<-rep(track$parent[i], length(x))
    
    r<-rep(track$root[i], length(x))
    
    g<-rep(track$generation[i], length(x))
    
    tr<-unlist(strsplit(track$all[i], split=":"))
    a<-as.numeric(unlist(strsplit(tr[12], split=",|\\[|\\]| ")))
    a<-a[!is.na(a)]
    
    f<-rep(track$fate[i], length(x))
    
    d<-data.frame(t=t,x=x,y=y,z=z,
                  label=label,
                  area=a,
                  parent=p,
                  root=r,
                  generation=g,
                  fate=f)
  }
  track_pos<-foreach(i=1:nrow(track_pos)) %do% {
    
    t<-as.numeric(unlist(strsplit(track_pos$t[i], split=",|\\[|\\]| ")))
    t<-t[!is.na(t)]
    
    x<-as.numeric(unlist(strsplit(track_pos$x[i], split=",|\\[|\\]| ")))
    x<-x[!is.na(x)]
    
    y<-as.numeric(unlist(strsplit(track_pos$y[i], split=",|\\[|\\]| ")))
    y<-y[!is.na(y)]
    
    z<-rep(0, length(x))
    
    label<-rep(i, length(x))
    
    p<-rep(track_pos$parent[i], length(x))
    
    r<-rep(track_pos$root[i], length(x))
    
    g<-rep(track_pos$generation[i], length(x))
    
    tr<-tr<-unlist(strsplit(track_pos$all[i], split=":"))
    a<-as.numeric(unlist(strsplit(tr[12], split=",|\\[|\\]| ")))
    a<-a[!is.na(a)]
    
    f<-rep(track_pos$fate[i], length(x))
    
    d<-data.frame(t=t,x=x,y=y,z=z,
                  label=label,
                  area=a,
                  parent=p,
                  root=r,
                  generation=g,
                  fate=f)
  }
  
  a<-unlist(lapply(1:length(track), function(x) track[[x]]$area[1]))
  l<-unlist(lapply(1:length(track), function(x) nrow(track[[x]])))
  track<-track[a>10&l>1]
  
  a<-unlist(lapply(1:length(track_pos), function(x) track_pos[[x]]$area[1]))
  l<-unlist(lapply(1:length(track_pos), function(x) nrow(track_pos[[x]])))
  track_pos<-track_pos[a>10&l>1]
  
  track.n<-track
  
  for (i in 1:length(track.n)){
    track.n[[i]]$label<-rep(track.n[[i]]$label[1],nrow(track.n[[i]]))
  }
  
  
  np.track<-foreach(i=1:length(track.n),.combine=rbind) %do% {
    track.n[[i]]
  }
  
  np.track<-np.track[,c("label","z","y","x")]
  
  
  #calculate distance
  track.dist<-foreach(i=1:length(track),.combine=rbind) %do% {
    path_dist(track[[i]]$x,track[[i]]$y)}
  
  track.dist<-data.frame(track.dist)
  colnames(track.dist)<-c("path.distance","overall.distance")
  
  track.dist<-track.dist*dist.factor
  
  track.dist.pos<-foreach(i=1:length(track_pos),.combine=rbind) %do% {
    path_dist(track_pos[[i]]$x,track_pos[[i]]$y)}
  
  track.dist.pos<-data.frame(track.dist.pos)
  colnames(track.dist.pos)<-c("path.distance","overall.distance")
  track.dist.pos<-track.dist.pos*dist.factor
  
  
  dim(track.dist)
  dim(track.dist.pos)
  
  n<-min(nrow(track.dist),nrow(track.dist.pos))
  
  track.dist<-track.dist[sample(rownames(track.dist), n),]
  track.dist.pos<-track.dist.pos[sample(rownames(track.dist.pos), n),]
  
  plot(density(track.dist$path.distance), col="blue")
  lines(density(track.dist.pos$path.distance), col="red")
  
  plot(density(track.dist$overall.distance), col="blue", xlim=c(0,300))
  lines(density(track.dist.pos$overall.distance), col="red")

  d<-data.frame(neg=track.dist$path.distance, pos=track.dist.pos$path.distance, line=line)
  return(d)
  
}


track<-read.csv("~/Desktop/R/napari/tracks_356-61neg.csv", header=F)
track_pos<-read.csv("~/Desktop/R/napari/tracks_356-61pos.csv", header=F)

colnames(track)<-c("all","ID", "t", "x", "y", "root", "parent", "children", "generation", "fate")
colnames(track_pos)<-c("all","ID", "t", "x", "y", "root", "parent", "children", "generation", "fate")


track<-foreach(i=1:nrow(track)) %do% {

  t<-as.numeric(unlist(strsplit(track$t[i], split=",|\\[|\\]| ")))
  t<-t[!is.na(t)]
  
  x<-as.numeric(unlist(strsplit(track$x[i], split=",|\\[|\\]| ")))
  x<-x[!is.na(x)]
  
  y<-as.numeric(unlist(strsplit(track$y[i], split=",|\\[|\\]| ")))
  y<-y[!is.na(y)]
  
  z<-rep(0, length(x))
  
  label<-rep(i, length(x))
  
  p<-rep(track$parent[i], length(x))
  
  r<-rep(track$root[i], length(x))
  
  g<-rep(track$generation[i], length(x))
  
  tr<-tr<-unlist(strsplit(track$all[i], split="[(]"))
  a<-as.numeric(unlist(strsplit(tr[14], split=",|\\[|\\]| ")))
  a<-a[!is.na(a)]
  
  f<-rep(track$fate[i], length(x))
  
  d<-data.frame(t=t,x=x,y=y,z=z,
                label=label,
                area=a,
                parent=p,
                root=r,
                generation=g,
                fate=f)
}
track_pos<-foreach(i=1:nrow(track_pos)) %do% {
  
  t<-as.numeric(unlist(strsplit(track_pos$t[i], split=",|\\[|\\]| ")))
  t<-t[!is.na(t)]
  
  x<-as.numeric(unlist(strsplit(track_pos$x[i], split=",|\\[|\\]| ")))
  x<-x[!is.na(x)]
  
  y<-as.numeric(unlist(strsplit(track_pos$y[i], split=",|\\[|\\]| ")))
  y<-y[!is.na(y)]
  
  z<-rep(0, length(x))
  
  label<-rep(i, length(x))
  
  p<-rep(track_pos$parent[i], length(x))
  
  r<-rep(track_pos$root[i], length(x))
  
  g<-rep(track_pos$generation[i], length(x))
  
  tr<-tr<-unlist(strsplit(track_pos$all[i], split="[(]"))
  a<-as.numeric(unlist(strsplit(tr[14], split=",|\\[|\\]| ")))
  a<-a[!is.na(a)]
  
  f<-rep(track_pos$fate[i], length(x))
  
  d<-data.frame(t=t,x=x,y=y,z=z,
                label=label,
                area=a,
                parent=p,
                root=r,
                generation=g,
                fate=f)
}




a<-unlist(lapply(1:length(track), function(x) track[[x]]$area[1]))
l<-unlist(lapply(1:length(track), function(x) nrow(track[[x]])))
track<-track[a>10&l>1]

a<-unlist(lapply(1:length(track_pos), function(x) track_pos[[x]]$area[1]))
l<-unlist(lapply(1:length(track_pos), function(x) nrow(track_pos[[x]])))
track_pos<-track_pos[a>10&l>1]



j<-0
idds<-1
while (length(idds)>0) {
  
  a<-unlist(lapply(1:length(track), function(x) track[[x]]$area[nrow(track[[x]])]))
  track<-track[a>10]
  a<-unlist(lapply(1:length(track), function(x) track[[x]]$area[nrow(track[[x]])]))
  d<-unlist(lapply(1:length(track), function(x) nrow(track[[x]])))
  t1<-unlist(lapply(1:length(track), function(x) track[[x]]$t[1]))
  tend<-unlist(lapply(1:length(track), function(x) track[[x]]$t[nrow(track[[x]])]))
  x1<-unlist(lapply(1:length(track), function(x) track[[x]]$x[1]))
  xend<-unlist(lapply(1:length(track), function(x) track[[x]]$x[nrow(track[[x]])]))
  y1<-unlist(lapply(1:length(track), function(x) track[[x]]$y[1]))
  yend<-unlist(lapply(1:length(track), function(x) track[[x]]$y[nrow(track[[x]])]))
  l<-unlist(lapply(1:length(track), function(x) track[[x]]$label[1]))
  f1<-unlist(lapply(1:length(track), function(x) track[[x]]$fate[1]))
  
  
  
  di<-which(tend!=max(tend)&f1!="Fates.DIVIDE"&f1!="Fates.TERMINATE_BORDER")
  ids<-c()
  for (i in di){
      candidates<-which(t1==tend[i]+1)
      candidates<-setdiff(candidates,ids)
        if (length(candidates)>0){
            d<-ref.dist(x1[candidates],y1[candidates],c(xend[i],yend[i]))
            area<-a[candidates]
            ai<-which(area>a[i]*0.8&area<a[i]*1.2)
            d<-d[ai]
            if (length(d)>0){
              if (min(d)<100){
                id<-candidates[which(d==min(d))]
                track[[i]]<-rbind(track[[i]],track[[id]])
                ids<-c(ids,id)
              }
            }
        }
  }
  idds<-ids
  track<-track[setdiff(1:length(track),ids)]
  j=j+1
  print(j)
}

track.n<-track

for (i in 1:length(track.n)){
  track.n[[i]]$label<-rep(track.n[[i]]$label[1],nrow(track.n[[i]]))
}



np.track<-foreach(i=1:length(track.n),.combine=rbind) %do% {
  track.n[[i]]
}

np.track<-np.track[,c("label","z","y","x")]

write.table(np.track, file="~/Desktop/R/napari/nptrack.txt", row.names=F, sep=",")

#track[[i]]$t[1]<10 & 

# id=c()
# for (i in 1:length(track)){
#   if (track[[i]]$t[1]<10 & nrow(track[[i]])>10){
#     id=c(id,i)
#   }
# }
# track<-track[id]
# 
# id=c()
# for (i in 1:length(track_pos)){
#   if (track_pos[[i]]$t[1]<10 & nrow(track_pos[[i]])>10){
#     id=c(id,i)
#   }
# }
# track_pos<-track_pos[id]
# 


#calculate distance
track.dist<-foreach(i=1:length(track),.combine=rbind) %do% {
  path_dist(track[[i]]$x,track[[i]]$y)}

track.dist<-data.frame(track.dist)
colnames(track.dist)<-c("path.distance","overall.distance")

#dist.factor<-0.908
#track.dist<-track.dist*dist.factor

track.dist.pos<-foreach(i=1:length(track_pos),.combine=rbind) %do% {
  path_dist(track_pos[[i]]$x,track_pos[[i]]$y)}

track.dist.pos<-data.frame(track.dist.pos)
colnames(track.dist.pos)<-c("path.distance","overall.distance")



dim(track.dist)
dim(track.dist.pos)

n<-min(nrow(track.dist),nrow(track.dist.pos))

track.dist<-track.dist[sample(rownames(track.dist), n),]
track.dist.pos<-track.dist.pos[sample(rownames(track.dist.pos), n),]

plot(density(track.dist$path.distance), col="blue")
lines(density(track.dist.pos$path.distance), col="red")

plot(density(track.dist$overall.distance), col="blue", xlim=c(0,300))
lines(density(track.dist.pos$overall.distance), col="red")


RPE318<-data.frame(CD61neg=track.dist$path.distance, CD61pos=track.dist.pos$path.distance, line="RPE318")
RPE319<-data.frame(CD61neg=track.dist$path.distance, CD61pos=track.dist.pos$path.distance, line="RPE319")
RPE322<-data.frame(CD61neg=track.dist$path.distance, CD61pos=track.dist.pos$path.distance, line="RPE322")


dd<-rbind(melt(RPE318),melt(RPE319),melt(RPE322))
d<-rbind(melt(R356),melt(R358))

p2<-ggplot(dd, aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,9))
p2+ theme_minimal()


#ECDF
library(twosamples)
set.seed(31415)
vec1=dd$value[which(dd$variable=="CD61neg")]
vec2=dd$value[which(dd$variable=="CD61pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()

ec

# P0

p0<-ggplot(d, aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,9))
p0+ theme_minimal()


#ECDF
library(twosamples)
set.seed(31415)
vec1=d$value[which(d$variable=="CD61neg")]
vec2=d$value[which(d$variable=="CD61pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()











dp <- ggplot(melt(RPE318), aes(x=variable, y=value, fill=line)) + 
  geom_violin(trim=FALSE)+
  labs(title="Plot of length  by dose",x="Dose (mg)", y = "Length")
dp + theme_classic()
dp+scale_fill_brewer(palette="Blues") + theme_classic()


p<-ggplot(d[d$line=="RPE356",], aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,9))
p+ theme_minimal()

p+ theme_classic()

p + scale_color_brewer(palette="Dark2") + theme_minimal()
p + scale_color_brewer(palette="Accent") + theme_classic()


#ECDF
library(twosamples)
set.seed(31415)
vec1=RPE319$CD61neg
vec2=RPE319$CD61pos
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()

ec





cul<-data.frame(neg=c(mean(RPE318$CD61neg),mean(RPE319$CD61neg),mean(RPE322$CD61neg)),
                pos=c(mean(RPE318$CD61pos),mean(RPE319$CD61pos),mean(RPE322$CD61pos)))

cul<-melt(cul)


R356<-data.frame(CD61neg=track.dist$path.distance, CD61pos=track.dist.pos$path.distance, line="RPE356")
R358<-data.frame(CD61neg=track.dist$path.distance, CD61pos=track.dist.pos$path.distance, line="RPE358")

pr<-data.frame(neg=c(mean(R356$CD61neg),mean(R358$CD61neg)),
               pos=c(mean(R356$CD61pos),mean(R358$CD61pos)))
pr<-melt(pr)

d<-cul %>% group_by(variable) %>% summarise( 
  n=n(),
  mean=mean(value),
  sd=sd(value))
  
  
culp<-ggplot(d)+
  geom_bar(aes(x=variable, y=mean), stat="identity", fill="black", alpha=1)+
  geom_errorbar( aes(x=variable, ymin=mean, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1.3)+
  ylim(c(0,350))+
  theme_classic()



dd<-rbind(melt(RPE318),melt(RPE319),melt(RPE322))
d<-rbind(melt(R356),melt(R358))

#ECDF cultured

library(twosamples)
set.seed(31415)
vec1=RPE319$CD61neg
vec2=RPE319$CD61pos
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()

ec








d<-pr %>% group_by(variable) %>% summarise( 
  n=n(),
  mean=mean(value),
  sd=sd(value))


prp<-ggplot(d)+
  geom_bar(aes(x=variable, y=mean), stat="identity", fill="black", alpha=1)+
  geom_errorbar( aes(x=variable, ymin=mean, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1.3)+
  ylim(c(0,350))+
  theme_classic()




mean(track.dist$path.distance)
sd(track.dist$path.distance)
mean(track.dist.pos$path.distance)
sd(track.dist.pos$path.distance)


a<-data.frame(cd61=c("pos","pos","neg","neg"), cat=c("<150",">150","<150",">150"),
              dist=c(mean(track.dist.pos[which(track.dist.pos$path.distance<150),"path.distance"]),mean(track.dist.pos[which(track.dist.pos$path.distance>150),"path.distance"]),
                                                                             mean(track.dist[which(track.dist$path.distance<150),"path.distance"]),mean(track.dist[which(track.dist$path.distance>150),"path.distance"])),
              freq=c(length(which(track.dist.pos$path.distance<150)),length(which(track.dist.pos$path.distance>150)),
                     length(which(track.dist$path.distance<150)),length(which(track.dist$path.distance>150))),
              sd=c(sd(track.dist.pos[which(track.dist.pos$path.distance<150),"path.distance"]),sd(track.dist.pos[which(track.dist.pos$path.distance>150),"path.distance"]),
                     sd(track.dist[which(track.dist$path.distance<150),"path.distance"]),sd(track.dist[which(track.dist$path.distance>150),"path.distance"])))


mean(track.dist$path.distance)
sd(track.dist$path.distance)
mean(track.dist.pos$path.distance)
sd(track.dist.pos$path.distance)



  
a<-a18  
a$freq<-(a18$freq+a22$freq+a19$freq)/3
a$dist<-(a18$dist+a19$dist+a22$dist)/3
a$sd<-(a18$sd+a19$sd+a22$sd)/3
a$sdf<-c(sd(c(a18$freq[1],a19$freq[1],a22$freq[1])),
         sd(c(a18$freq[2],a19$freq[2],a22$freq[2])),
         sd(c(a18$freq[3],a19$freq[3],a22$freq[3])),
         sd(c(a18$freq[4],a19$freq[4],a22$freq[4])))

a<-a[c(3,4,1,2),]
a$cd61<-factor(a$cd61, levels = c("neg","pos"))
a$cat<-factor(a$cat, levels = c(">150","<150"))
a$fpercent<-c(100*(a$freq[1:2]/sum(a$freq[1:2])),100*(a$freq[3:4]/sum(a$freq[3:4])))

ggplot(a, aes(x=cat, y=freq, fill=cd61)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  theme_classic()

brewer_palette = brewer.pal(2, "Set1")
ggplot(a, aes(x = cd61, y = fpercent, fill = cat)) +
  geom_bar(stat = "identity") +
  labs(y = "frequency") +
  scale_fill_manual(values = brewer_palette)+
  theme_classic()


ggplot(a, aes(x=cat, y=dist, fill=cd61)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=dist, ymax=dist+sd), width=.2,
                position=position_dodge(.9))+
  theme_classic()

  ggplot(a, aes(x=cd61, y=freq, fill=cat)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=dist, ymax=dist+sd), width=.2,
                position=position_dodge(.9))+
  theme_classic()



ggplot(a, aes(x=line, y=integ, fill=cd)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=integ, ymax=integ+sd), width=.2,
                position=position_dodge(.9))+
  theme_classic()


ggplot(a, aes(x=cd61, y=dist, fill=cd61)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=dist, ymax=dist+sd), width=.2,
                position=position_dodge(.9))+
  theme_classic()



###-------------------- CD61 ---------------------------
CD61.318<-btrackv1(neg="~/Desktop/R/napari/tracks_318-61neg.csv",
                  pos="~/Desktop/R/napari/tracks_318-61pos.csv",
                  line = "RPE318",
                  dist.factor =0.908)

CD61.319<-btrackv1(neg="~/Desktop/R/napari/319_neg.csv",
                  pos="~/Desktop/R/napari/319_pos.csv",
                  line = "RPE319",
                  dist.factor =0.908)

CD61.322<-btrackv1(neg="~/Desktop/R/napari/tracks_322-61neg.csv",
                  pos="~/Desktop/R/napari/tracks_322-61pos.csv",
                  line = "RPE322",
                  dist.factor =0.908)

d261<-rbind(melt(CD61.318),melt(CD61.319),melt(CD61.322))

p261<-ggplot(d261, aes(x=(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-20,1000))
p261+ theme_minimal()


p261<-ggplot(d261, aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,10))
p261+ theme_minimal()


#ECDF
library(twosamples)
set.seed(31415)
vec1=d261$value[which(d261$variable=="neg")]
vec2=d261$value[which(d261$variable=="pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()

ec

#P0

CD61.356<-btrack(neg="~/Desktop/R/napari/115-tracking/tracks_356-61neg.csv",
                  pos="~/Desktop/R/napari/115-tracking/tracks_356-61pos.csv",
                  line = "RPE356",
                  dist.factor =1.36)

CD61.358<-btrack(neg="~/Desktop/R/napari/115-tracking/tracks_358-61neg.csv",
                  pos="~/Desktop/R/napari/115-tracking/tracks_358-61pos.csv",
                  line = "RPE358",
                  dist.factor =0.908)


d061<-rbind(melt(CD61.356),melt(CD61.358))

p061<-ggplot(d061, aes(x=(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-50,1000))
p061+ theme_minimal()


p061<-ggplot(d061, aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,10))
p061+ theme_minimal()


#ECDF
library(twosamples)
set.seed(31415)
vec1=d061$value[which(d061$variable=="neg")]
vec2=d061$value[which(d061$variable=="pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()
ec










###-------------------- CD115 --------------------------

CD115.318<-btrack(neg="~/Desktop/R/napari/115-tracking/tracks_318-115neg.csv",
                  pos="~/Desktop/R/napari/115-tracking/tracks_318-115pos.csv",
                  line = "RPE318",
                  dist.factor =0.908)

CD115.319<-btrack(neg="~/Desktop/R/napari/115-tracking/tracks_319-115neg.csv",
                  pos="~/Desktop/R/napari/115-tracking/tracks_319-115pos.csv",
                  line = "RPE319",
                  dist.factor =0.908)

CD115.322<-btrack(neg="~/Desktop/R/napari/115-tracking/tracks_322-115neg.csv",
                  pos="~/Desktop/R/napari/115-tracking/tracks_322-115pos.csv",
                  line = "RPE322",
                  dist.factor =0.908)


d115<-rbind(melt(CD115.318),melt(CD115.319),melt(CD115.322))

p115<-ggplot(d115, aes(x=(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-20,1000))
p115+ theme_minimal()

p115<-ggplot(d115, aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,10))
p115+ theme_minimal()



#ECDF
library(twosamples)
set.seed(31415)
vec1=d115$value[which(d115$variable=="neg")]
vec2=d115$value[which(d115$variable=="pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()

ec




# P0



CD115.356<-btrack(neg="~/Desktop/R/napari/115-tracking/tracks_356-115neg.csv",
                  pos="~/Desktop/R/napari/115-tracking/tracks_356-115pos.csv",
                  line = "RPE356",
                  dist.factor =1.36)

CD115.358<-btrack(neg="~/Desktop/R/napari/115-tracking/tracks_358-115neg.csv",
                  pos="~/Desktop/R/napari/115-tracking/tracks_358-115pos.csv",
                  line = "RPE358",
                  dist.factor =0.908)


d015<-rbind(melt(CD115.356),melt(CD115.358))

p015<-ggplot(d015, aes(x=(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-50,1000))
p015+ theme_minimal()

p015<-ggplot(d015, aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,10))
p015+ theme_minimal()


#ECDF
library(twosamples)
set.seed(31415)
vec1=d015$value[which(d015$variable=="neg")]
vec2=d015$value[which(d015$variable=="pos")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()

ec




d61pos<-d261
d61pos$variable<-as.character(d61pos$variable)
d61pos$variable[which(d61pos$variable=="pos")]<-"pos2"
d61pos$variable<-as.factor(d61pos$variable)
d61pos<-rbind(d61pos,d061)
d61pos<-d61pos[which(d61pos$variable!="neg"),]

table(d61pos$variable)

d61pos <- d61pos %>%
  group_by(variable) %>%
  sample_n(486)

table(d61pos$variable)


p61pos<-ggplot(d61pos, aes(x=(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-50,1000))
p61pos+ theme_minimal()

p61pos<-ggplot(d61pos, aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,10))
p61pos+ theme_minimal()



library(twosamples)
set.seed(31415)
vec1=d61pos$value[which(d61pos$variable=="pos")]
vec2=d61pos$value[which(d61pos$variable=="pos2")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("neg", length(vec1)),rep("pos", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()

ec








d15pos<-d115
d15pos$variable<-as.character(d15pos$variable)
d15pos$variable[which(d15pos$variable=="pos")]<-"pos2"
d15pos$variable<-as.factor(d15pos$variable)
d15pos<-rbind(d15pos,d015)
d15pos<-d15pos[which(d15pos$variable!="neg"),]
table(d15pos$variable)

d15pos <- d15pos %>%
  group_by(variable) %>%
  sample_n(501)

table(d15pos$variable)



p15pos<-ggplot(d15pos, aes(x=(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-50,1000))
p15pos+ theme_minimal()

p15pos<-ggplot(d15pos, aes(x=log(value), fill=variable)) +
  geom_density(alpha=0.4)+
  xlim(c(-2,10))
p15pos+ theme_minimal()



library(twosamples)
set.seed(31415)
vec1=d15pos$value[which(d15pos$variable=="pos")]
vec2=d15pos$value[which(d15pos$variable=="pos2")]
out=wass_test(vec1,vec2)
out
summary(out)
plot(out)

df<-data.frame(cd=c(vec1,vec2), meta=c(rep("pos", length(vec1)),rep("pos2", length(vec2))))
ec<-ggplot(df, aes(x=cd, fill=meta, color=meta)) +
  stat_ecdf(geom = "step")+
  xlim(c(0,1000))+scale_color_manual(values=c("red","blue"))+
  theme_classic()

ec



grid.arrange(p115, p015, ncol=2)
grid.arrange(p261, p061, ncol=2)

grid.arrange(p115, p015, p15pos, p261, p061, p61pos, ncol=3)

p261+theme_minimal()+p061+theme_minimal()+
  p115+theme_minimal()+ p015+theme_minimal()

