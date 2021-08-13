# Gap Analysis


###### This  (through line 36) combines a couple different datasets (trims, changes column names etc.) so probably won't be relevant
# Take CNBH & personal records. Combine and clean for analysis
#Symbiota native download from a polygon search
cnbh<-read.csv('~/Desktop/CURRENT PROJECTS/GAP_analysis/CNBH_allIslands_2020_12_04.csv',stringsAsFactors=F,fileEncoding="latin1")
sjsu<-read.csv('~/Desktop/CURRENT PROJECTS/GAP_analysis/GapAnalysis.csv',stringsAsFactors=F)

setwd("~/Desktop/CURRENT PROJECTS/GAP_analysis")
 ######### match columns and combine datasets
 cnbh<-cnbh[,                  
c("id","catalogNumber", "verbatimElevation",
"decimalLatitude", "decimalLongitude","verbatimCoordinates",
"locality","county","scientificName",
"scientificNameAuthorship","class","order",
"family", "genus", "specificEpithet", 
"infraspecificEpithet","eventDate", "year",
"month", "day", "verbatimEventDate",
"recordedBy", "recordNumber")]

xxx<-rep("NA",dim(sjsu)[1])
sjsu<-cbind(
xxx, sjsu$SJSU_AccessionNumber, sjsu$Elev_m,
sjsu$Latitude, sjsu$Longitude, xxx,
sjsu$General.description, sjsu$County, sjsu$Scientific_Name, 
sjsu$Author, xxx, xxx,
xxx, sjsu$Genus, sjsu$Species,
sjsu$infra.name, sjsu$Date, sjsu$year,
sjsu$month, sjsu$day, sjsu$Date,
sjsu$Collector, sjsu$Col.Num)

colnames(sjsu)<-colnames(cnbh)
all<-rbind(cnbh,sjsu)

############### Clean the dataset #######################


## This is specific to my dataset (I used overlapping datasets)- these aren't herbarium duplicates, they're database duplicates
########### remove duplicates 
d<-duplicated(all$catalogNumber)
all<-all[!d,]

############## drop points with meaningless coords- this removes any specimens that have long/lats that are not in the vicinity of the islands
w<-which(as.numeric(all$decimalLongitude) > (-118) |
	as.numeric(all$decimalLongitude) < (-120.5))
all<-all[-w,]
w<-which(as.numeric(all$decimalLatitude) > (34.1) |
	as.numeric(all$decimalLatitude) < (32.7))
all<-all[-w,]

########## remove obvious georeferencing errors by hand- the 'identify' function allows you to click on a scatterplot (ie distribution map) to identify individual points. This is nice if you have a few coords out in the ocean. If you have a bunch, you'll probably want to do an intersection with island polygons, but I didn't have an island polygon file when I was doing this (the maps I used are line shapefiles, not polygons)
plot(all$decimalLongitude,all$decimalLatitude,asp=1)
identify(all$decimalLongitude,all$decimalLatitude,n=1)
# is singleton?
w<-which(all$decimalLongitude==all$decimalLongitude[975] &
	     all$decimalLatitude==all$decimalLatitude[975])
length(w)
all<-all[-w,]

########## assign records to islands ###
# You'll need to have 8 presence/absence columns (1/0) corresponding to the islands. My approach was to just do a text search using 'grep' since all the specimens I'm working with have the island name somewhere in the locality description. You may need to take another approach, but the point is to end up with an 8 column matrix of 1s & 0s with the number of rows equal to the number of specimens (NOT the number of taxa)

isl.rec<-matrix(nrow=dim(all)[1],ncol=8)
isl.rec[]<-0
isl.names<-c("Anacapa","Santa Cruz","Santa Rosa","San Miguel",
"Santa Catalina","San Clemente","San Nicolas","Santa Barbara")
for(i in 1:8){isl.rec[grep(isl.names[i],all$locality,ignore.case=T),i]<-1}
isl.rec[grep("San Nicholas",all$locality,ignore.case=T),7]<-1
isl.rec[grep("Catalina",all$locality,ignore.case=T),5]<-1

### find the issues- any specimen with a rowSum of zero didn't get assigned to an island, so go back and find/fix these. You'll need to customize this to your dataset
table(rowSums(isl.rec))

to.resolve<-all[rowSums(isl.rec)==0,]
to.resolve$locality

#make some custom replacements
w<-which(rowSums(isl.rec)==0)

isl.rec[w[c(7:27)],3]<-1
isl.rec[w[5],5]<-1
isl.rec[w[1],6]<-1
table(rowSums(isl.rec))

all<-all[rowSums(isl.rec)==1,]
isl.rec<-isl.rec[rowSums(isl.rec)==1,]

################# make an island name column #############
# this converts the 8 column matrix in to a single column with just the island name for each specimen
island<-vector(length=dim(all)[1])
for(i in 1:8){island[isl.rec[,i]==1]<-isl.names[i]}
all<-cbind(all,island)

################## remove records with no epithet #########
# you might also want to do this for taxa identified only to genus
w<-which(all$specificEpithet=="")
all<-all[-w,]
table(all$island)

############ make a binomial vector #######
# make sure you have a column with the full taxon name (binomials or trinomials)
binom<-paste(all$genus,all$specificEpithet,sep="_")
all<-cbind(binom,all)


##################### Done with Cleanup ####################


##########################################################
#################### Spatial #############################
# make a raster
# convert islands to pixels


############ Species maps
# You probably have access to better maps- I just got these off the internet. I used Albers for the mapping so that I could have grid cells in km^2 instead of degrees.
counties<-readOGR("~/Desktop/CURRENT PROJECTS/GAP_analysis/CAcoastline",layer="cnty24k09_1_line")
subco<-counties[counties$NAME%in%c("Ventura", "Santa Barbara","Los Angeles","San Diego","Orange"),]

# crop to decrease rendering time
is.map<-crop(counties, extent(-1.5e+05,2.5e+05,-6e+05,-4e+05))
plot(is.map,col="gray50")
rm(counties)
rm(subco)


####### project long lat coords into Albers (island maps are in Albers)

P4S.latlon <- CRS("+proj=longlat +datum=WGS84")
new_proj<-"+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
sp<-SpatialPointsDataFrame(coords= cbind(as.numeric(all$decimalLongitude),as.numeric(all$decimalLatitude)),
	data=all,proj4string=P4S.latlon)
pointsX<-spTransform(sp,CRS(new_proj))

# attach projected points to 'all'
all<-cbind(all,pointsX@coords)

######## Make some plots
plot(is.map,col="gray50")
points(all$coords.x1,all$coords.x2,cex=.3)

### plotting limits for each island- you can use these to plot one island at a time, since R doesn't have a zoom in/out feature
plot(is.map,col="gray50",axes=T)
abline(v=seq(from=-1.5e+05,to=2.5e+05,by=5e+04))
abline(h=seq(from=-6e+05,to=-4e+05,by=5e+04))

island.plot.lims<-matrix(nrow=10,ncol=4)
colnames(island.plot.lims)<-c("xmn","xmx","ymn","ymx")
rownames(island.plot.lims)<-c(isl.names,"north","south")
island.plot.lims[1,]<-c(50000,60000,-450000,-440000)
island.plot.lims[2,]<-c(0,50000,-460000,-420000)
island.plot.lims[3,]<-c(-30000,10000,-470000,-430000)
island.plot.lims[4,]<-c(-40000,-20000,-450000,-430000)
island.plot.lims[5,]<-c(120000,160000,-530000,-500000)
island.plot.lims[6,]<-c(130000,160000,-590000,-550000)
island.plot.lims[7,]<-c(30000,60000,-540000,-520000)
island.plot.lims[8,]<-c(80000,100000,-510000,-500000)
island.plot.lims[9,]<-c(-50000,70000,-470000,-420000)
island.plot.lims[10,]<-c(30000,160000,-580000,-500000)

write.csv(island.plot.lims,"islandPlotLimits.csv")
####### individual island maps ###########
# i corresponds to the columns in the 8 column matrix
i=7
plot(is.map,
	xlim=island.plot.lims[i,c(1,2)],
	ylim=island.plot.lims[i,c(3,4)],
	col="gray50",asp=1)
points(all$coords.x1,all$coords.x2,cex=.6,col=rgb(0,0,1,0.5),pch=19)

legend("bottomleft",
legend=c(
paste("Species:",length(unique(all$binom[all$island==isl.names[i]]))),
paste("Collections:",length(which(all$island==isl.names[i])))),
col=0)


############### Make Rasters #############
# for collecting intensity & richness grids, make a raster with some arbitrary grid cell size and then find the grid cell for each specimen. 
# 5km grid cells: ncol=50,nrow=34
# 10km: ncol=25,nrow=17
# 1km = ncol= 250, nrow=170
r<-raster(
xmn=-50000, 
xmx=200000, 
ymn=-600000, 
ymx=-430000,
ncol= 250, nrow=170)
r[]<-0

### First, rasterize the islands. This would be a lot easier with a polygon map but I didn't have one so I did some work by hand here- I rasterized the outlines, then selected grid cells by hand from the middle of islands to make up the island rasters. It would be much easier to just rasterize a polygon shapefile.
isl.rast<-rasterize(is.map,r)

# now go through, island by island, and fill in the holes
# (rasterize doesn't quite work here bc it's lines, not polygons)
i=8
plot(isl.rast,
	xlim=island.plot.lims[i,c(1,2)],
	ylim=island.plot.lims[i,c(3,4)],
	col="gray50",asp=1)

plot(is.map,add=T)

#identify cells to add with click()
to.add<-click(isl.rast,n=4,cell=TRUE)

# add selected to the raster
isl.rast[to.add$cell]<-1

# convert to binary
w<-which(is.na(isl.rast[])==FALSE)
isl.rast[w]<-1

############### identify cells for each island ##########
island.pix<-list()

i=5
plot(isl.rast,
	xlim=island.plot.lims[i,c(1,2)],
	ylim=island.plot.lims[i,c(3,4)],
	col="gray50",asp=1)

# select the shape
s<-select(isl.rast,use='pol',draw=T)
# convert to coord
s.coords<-xyFromCell(s,which(s[]==1))
# intersect coords with original rast
s.cells<-cellFromXY(isl.rast,s.coords)
island.pix[[i]]<-s.cells

# now convert to matrix and save

island.cells<-data.frame(matrix(nrow=0,ncol=2))
colnames(island.cells)<-c("island","isl.rastCell")
for(i in 1:8)
{
	x<-island.pix[[i]]
	island.cells<-rbind(island.cells,cbind(rep((isl.names[i])),x))
}

write.csv(island.cells,"islandCells.csv",row.names=F)
writeRaster(isl.rast,"islandRaster.asc",format="ascii")
write.csv(all,"cleanOccurrences.csv",row.names=F)


