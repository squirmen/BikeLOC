list.of.packages <- c("data.table","DBI","doParallel","e1071",
                      "FNN","foreach","foreign","geosphere","ggmap","ggplot2","ggrepel","gmapsdistance",
                      "graphics","Grid2Polygons","gridExtra","Imap","kernlab","maps","maptools","nnet",
                      "osmar","osrm","parallel","plm","plyr","progress","rangeMapper","raster","RDSTK","rgdal",
                      "rgeos","RPostgreSQL","RSQLite","snow","sp","spatstat","stringr","tm","utils","wordcloud")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(nnet)
library(graphics)
library(e1071)
library(plm)
library(kernlab)
library(data.table)
library(stringr)
library(parallel)
library(doParallel)
library(snow)
library(sp)
library(spatstat)
library(maptools)
library(rgeos)
library(ggplot2)
library(plyr)
library(RSQLite)
library(DBI)
library(foreign)
library(rgdal)
library(maps)
library(geosphere)
library(gridExtra)
library(ggmap)
library(wordcloud)
library(tm)
library(ggrepel)
library(gmapsdistance)
library(Imap)
library(utils)
library(osmar)
library(RPostgreSQL)
library(raster)
library(progress)
library(spatstat)
library(RDSTK)
library(osrm)
library(Grid2Polygons)
library(FNN)
library(rangeMapper)
library(foreach)
library(dplyr)

require(compiler)
enableJIT(3)

basedir<-"R:/BikeLOC_Typ"

setwd(basedir)

##Using cell centroids for distance or something else? (default is cell cenntroid, if using other shapefile enter name below)
cell_centroid<-"Name"   #<--name shapefile here

##name of the shapefile for the study area
studyarea<-"wmata_counties"   #<--name shapefile here

##What size grid do you want (in meters)
gridsize<-200

######################
##   BE Variables   ##
######################

#connect to databases
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, user = "postgres", dbname = "osm", host = "localhost")
con_rt <- dbConnect(drv, user = "postgres", dbname = "opbeumDB", host = "localhost")


#----> CREATE GRID CELLS<----

  ## Set-up Inputs/OUTPUTS

##ODs from real station data
  GRID_BIKESHARE_REALODs = fread(paste0("Input/","cabi_OD_01.csv"), sep=",",header = TRUE)
    names(GRID_BIKESHARE_REALODs)[names(GRID_BIKESHARE_REALODs)=="Origin"] <- "TERMINAL_N_O"
    names(GRID_BIKESHARE_REALODs)[names(GRID_BIKESHARE_REALODs)=="Destination"] <- "TERMINAL_N_D"


   #study area shapefile
  place.spr.original<-readOGR(dsn = "Input", layer = paste0(studyarea))
  place.spr.original.proj<-spTransform( place.spr.original, CRS( "+init=epsg:3347" ) ) 
  #place.spr <- spTransform(place.spr.original,CRS("+proj=longlat +datum=WGS84"))
  
  #Bike-share station locations
  bs.stations<-readOGR(dsn = "Input", layer = "Capital_Bike_Share_Locations") 

  #for (z in 51:nrow(place.sp)) { #option shape loop
#place.spr<-place.sp #use for specific area
#place.spr<-place.sp[z,] #use if shape loop

### define SpatialGrid object
place.bb<-bbox(place.spr.original.proj)
place.cs <- c(3.28084, 3.28084)*gridsize  # cell size 300m x 300m
# 1 ft = 3.28084 m
place.cc <- place.bb[, 1] + (place.cs /2)  # cell offset
place.cd <- ceiling(diff(t(place.bb))/place.cs)  # number of cells per direction
place.grd <- GridTopology(cellcentre.offset=place.cc, cellsize=place.cs, cells.dim=place.cd) #create the grid

##make grid topology into SP dataframe
place_sp_grd <- SpatialGridDataFrame(place.grd,
                                     data=data.frame(id=1:prod(place.cd)),
                                     proj4string=CRS(proj4string(place.spr.original.proj)))

##convert SP dataframe to polygons
all_cells<-Grid2Polygons(place_sp_grd)

##clip cells to polygon
all_cells <- gIntersection(place.spr.original.proj, all_cells, byid = TRUE, drop_lower_td = TRUE)

######SWAP PROJECTED WITH WGS 
#transform to WGS
all_cells.proj <- all_cells
all_cells<- spTransform(all_cells,CRS("+proj=longlat +datum=WGS84"))

all_cells_id <- as.data.frame(matrix(0, ncol = 1, nrow = length(all_cells)))
all_cells_id$V1<- 1:nrow(all_cells_id)
colnames(all_cells_id)[1] <- "PageNumber"


# Extract polygon ID's
all_cells_id <- sapply(slot(all_cells, "polygons"), function(x) slot(x, "ID"))

# Create dataframe with correct rownames
all_cells.df <- data.frame( ID=1:length(all_cells), row.names = all_cells_id)

all_cells <- SpatialPolygonsDataFrame(all_cells,all_cells.df)

all_cells <-spTransform(all_cells, CRS("+proj=longlat +datum=WGS84"))
colnames(all_cells@data)[1] <- "PageNumber"

if(cell_centroid=="Name"){
#get the centroid of the cells for routing
origins<- gCentroid(all_cells,byid=TRUE)
origins <- SpatialPointsDataFrame(origins,all_cells.df)
  colnames(origins@data)[1] <- "PageNumber"

}else{
## Set-up Inputs/OUTPUTS
origins<-readOGR(dsn = "CaBi_Archive", layer = "CaBi_Stations") #input origins shapefile
}


#writeOGR(obj=all_cells, dsn=".", layer="all_cells", driver="ESRI Shapefile") 

##
#Build Itersection Query
##
INT_FINAL<- data.frame(ref_count=integer(),
                       lat=numeric(), 
                       lon=numeric(),
                       PageNumber=integer(),
                       stringsAsFactors=FALSE) 

time1<-system.time({
  
  #Chunk grid cells 
  Splitrows=1000
  splitnumb_raw<-(nrow(all_cells)/Splitrows)
  if(nrow(all_cells)<Splitrows){
    Splitrows<-nrow(all_cells)
    splitnumb_l<-1
    leftover<-0
  }else{
    splitnumb_l<-floor(nrow(all_cells)/Splitrows)
    leftover<-splitnumb_raw-floor(splitnumb_l)
  }
  
  #Loops through cell chunks
  for(i in 1:splitnumb_l) {
    if(i==splitnumb_l){
      X1=(i-1)*Splitrows
      X2=Splitrows*i+Splitrows*leftover
    }else{i
      X1=(i-1)*Splitrows
      X2=Splitrows*i
    }
    
    #ith grid cell (row) bounding box and corners
    all_cells_cut<-all_cells[X1:X2,]
    all_cells_cut <-spTransform(all_cells, CRS("+proj=longlat +datum=WGS84"))
    
    bb<-bbox(all_cells_cut)
    
    ##Find intersections & Cul da sacks
    select<-"SELECT 
    ref_count,
    ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
    ST_X(ST_Transform(geom_vertex,4326)) AS lon 
    FROM at_2po_vertex
    WHERE  
    geom_vertex && 
    ST_MakeEnvelope("
    
    bb2<-paste(bb[1,1],",",bb[2,1],",",bb[1,2],",",bb[2,2],",","4283")
    q_tsig <- paste(select,bb2,")")
    
    #Pull SQL result
    tsig_result <- dbGetQuery(con_rt, q_tsig)
    
    #Convert node XY to spatial data frame
    node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                        proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    
    #make point coordinate system same as cells
    proj4string(node.spdf) <- CRS("+proj=longlat +ellps=WGS84") 
    proj4string(all_cells_cut) <- CRS("+proj=longlat +ellps=WGS84") 
    
    #Use just the PageNumber column
    all_cells_cut_red<-all_cells_cut[,1]
    
    #spatial join nodes with cells
    PageNumber = over(node.spdf,all_cells_cut_red)
    tsig_result2<-cbind(tsig_result,PageNumber)
    INT_FINAL=rbind(INT_FINAL,tsig_result2)
    
  }
})

print(time1)

#Add 1 to node directions
INT_FINAL$ref_count=INT_FINAL$ref_count+1

#Split intersections from cul de sacs
INTERS <- plyr::count(INT_FINAL[ which(INT_FINAL$ref_count>1), ], c('PageNumber'))
colnames(INTERS)[2] <- "Intersections"
INTERS3_4 <- plyr::count(INT_FINAL[ which(INT_FINAL$ref_count>2 & INT_FINAL$ref_count<5), ], c('PageNumber'))
colnames(INTERS3_4)[2] <- "Intersections3_4way"
INTERS5 <- plyr::count(INT_FINAL[ which(INT_FINAL$ref_count>4), ], c('PageNumber'))
colnames(INTERS5)[2] <- "Intersections_5plus"
INTCDS <- plyr::count(INT_FINAL[ which(INT_FINAL$ref_count==1), ], c('PageNumber'))
colnames(INTCDS)[2] <- "CulDeSacs"
LINKS <-aggregate(cbind(INT_FINAL$ref_count)~INT_FINAL$PageNumber, data=INT_FINAL, sum, na.rm=TRUE)
names(LINKS) <- c("PageNumber", "LINKS")

#Merge with full cell set
all_cells.df<-as.data.frame(all_cells[,1])
int.temp<-merge(all_cells.df,INTERS,all=T)
int.temp<-merge(int.temp,INTERS3_4,all=T)
int.temp<-merge(int.temp,INTERS5,all=T)
int.temp<-merge(int.temp,INTCDS,all=T)
int.temp<-merge(int.temp,LINKS,all=T)

#Convert NAs to 0s
int.temp[c("Intersections", "CulDeSacs","LINKS")][is.na(int.temp[c("Intersections", "CulDeSacs","LINKS")])] <- 0

int.temp$INTDns<-int.temp$Intersections/2.223948
int.temp$CDSDns<-int.temp$CulDeSacs/2.223948
int.temp$LINKDns<-int.temp$LINKS/2.223948
int.temp$CNR<-int.temp$LINKS/(int.temp$LINKS+int.temp$CulDeSacs)
int.temp$ALPHA<-abs((int.temp$LINKS-(int.temp$LINKS+int.temp$CulDeSacs)+1)/(2*(int.temp$LINKS+int.temp$CulDeSacs)-5))
int.temp$BETA<-int.temp$LINKS/(int.temp$Intersections+int.temp$CulDeSacs)
int.temp$GAMMA<-int.temp$LINKS/(3*((int.temp$Intersections+int.temp$CulDeSacs)-2))
int.temp$CYCL<-int.temp$LINKS-(int.temp$Intersections+int.temp$CulDeSacs)+2

###Get LEHD Data
LEHD_FINAL<- data.frame(tot_jobs=numeric(),
                        n11=numeric(),
                        n21=numeric(),
                        n22=numeric(),
                        n23=numeric(),
                        n31_33=numeric(),
                        n42=numeric(),
                        n44_45=numeric(),
                        n48_49=numeric(),
                        n51=numeric(),
                        n52=numeric(),
                        n53 =numeric(),
                        n54=numeric(),
                        n55=numeric(),
                        n56=numeric(),
                        n61=numeric(),
                        n62=numeric(),
                        n71=numeric(),
                        n72=numeric(),
                        n81=numeric(),
                        n92=numeric(),
                        PageNumber=integer(),
                        BlockAcres=numeric(),
                        Blockinter=numeric(),
                        stringsAsFactors=FALSE) 
time_LEHD<-system.time({
  
  #Chunk grid cells 
  Splitrows=100
  splitnumb_raw<-(nrow(all_cells)/Splitrows)
  if(nrow(all_cells)<Splitrows){
    Splitrows<-nrow(all_cells)
    splitnumb_l<-1
    leftover<-0
  }else{
    splitnumb_l<-floor(nrow(all_cells)/Splitrows)
    leftover<-splitnumb_raw-floor(splitnumb_l)
  }
  
  #Loops through cell chunks
  for(i in 1:splitnumb_l) {
    if(i==splitnumb_l){
      X1=(i-1)*Splitrows
      X2=Splitrows*i+Splitrows*leftover
    }else{i
      X1=(i-1)*Splitrows
      X2=Splitrows*i
    }
    
    #ith grid cell (row) bounding box and corners
    #all_cells <-spTransform(all_cells, CRS("+proj=longlat +datum=WGS84"))
    all_cells_cut<-all_cells[X1:X2,]
    
    bb<-bbox(all_cells_cut)
    
    ##Find employment
    select<-"SELECT 
    C000_y  AS  tot_jobs,
    CNS01_y AS	n11,
    CNS02_y AS	n21,
    CNS03_y AS	n22,
    CNS04_y AS	n23,
    CNS05_y AS	n31_33,
    CNS06_y AS	n42,
    CNS07_y AS	n44_45,
    CNS08_y AS	n48_49,
    CNS09_y AS	n51,
    CNS10_y AS	n52,
    CNS11_y AS	n53, 
    CNS12_y AS	n54,
    CNS13_y AS	n55,
    CNS14_y AS	n56,
    CNS15_y AS	n61,
    CNS16_y AS	n62,
    CNS17_y AS	n71,
    CNS18_y AS	n72,	
    CNS19_y AS	n81,
    CNS20_y AS	n92,
    ST_AsText(geom) AS geom
    
    FROM blocks_lehd_2011_us
    WHERE  
    geom && 
    ST_MakeEnvelope("
    
    bb2<-paste(bb[1,1],",",bb[2,1],",",bb[1,2],",",bb[2,2],",","4283")
    q_tsig <- paste(select,bb2,")")
    
    #Pull SQL result
    tsig_result <- dbGetQuery(con, q_tsig)
    tsig_result$ID<-row.names.data.frame(tsig_result)
    
    #tsig_result1<-tsig_result
    poly.spdf<-WKT2SpatialPolygonsDataFrame(tsig_result, 'geom','ID')
    #plot(poly.spdf)
    
    # #Convert node XY to spatial data frame
    # node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(10,9)], data = tsig_result,
    #                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    # 
    
    #make point coordinate system same as cells
    proj4string(poly.spdf) <- CRS("+proj=longlat +ellps=WGS84") 
    proj4string(all_cells_cut) <- CRS("+proj=longlat +ellps=WGS84") 
    #all_cells_cut.proj <- spTransform(all_cells_cut, CRS( "+init=epsg:3347" ) ) #for area
    poly.spdf.proj <- spTransform( poly.spdf, CRS( "+init=epsg:3347" ) ) #for area
  
    #Use just the PageNumber column
    
    all_cells_cut_red<-all_cells_cut[,1]
    
    #spatial join nodes with cells
    #PageNumber = over(poly.spdf,all_cells_cut_red)
    
    #X1<-gIntersection(poly.spdf, all_cells_cut_red, byid=T,id=all_cells_cut_red$PageNumber)
    #PageNumber<-sapply(X1@polygons, function(x) x@PageNumber)
    
    X1<-raster::intersect(poly.spdf, all_cells_cut_red)
    PageNumber<-as.data.frame(X1$PageNumber)
    ID<-as.data.frame(X1$ID)
    Res<-cbind(PageNumber,ID)
    colnames(Res)<-c("PageNumber","ID")
    Res<-unique(Res)
    
    all_block_area<-data.frame(round(gArea(poly.spdf.proj, byid = TRUE)*0.000247105,2))
    ABA_ID<-as.data.frame(poly.spdf.proj$ID)
    all_block_area<-cbind(all_block_area,ABA_ID)
    colnames(all_block_area) <- c("BlockAcres","ID")
    
    X1.proj <- spTransform( X1, CRS( "+init=epsg:3347" ) ) #for area
    
    all_block_intersect<-as.data.frame(round(gArea(X1.proj, byid = TRUE)*0.000247105,2))
    all_block_intersect<-cbind(all_block_intersect,ID)
    colnames(all_block_intersect) <- c("Blockinter","ID")
    
    blockprop<-merge(all_block_intersect,all_block_area,"ID",all.x=T)
    
    tsig_result2<-merge(blockprop,tsig_result,"ID",all.x=T)
    tsig_result3<-cbind(tsig_result2[c("tot_jobs",	"n11",	"n21",	"n22", "n23",	"n31_33",	"n42",	"n44_45",	"n48_49",	"n51",	"n52",	"n53",	"n54",	"n55",	"n56",	"n61",	"n62",	"n71",	"n72",	"n81",	"n92","Blockinter", "BlockAcres")],Res[1])
    LEHD_FINAL=rbind(LEHD_FINAL,tsig_result3)
    
  }
})

print(time_LEHD)

#get cell area
all_cells_area<-data.frame(round(gArea(all_cells.proj, byid = TRUE)*0.000247105,2))
all_cells_area$PageNumber<- 1:nrow(all_cells_area)
colnames(all_cells_area)[1] <- "CellAcres"

LEHD_FINAL<-merge(all_cells_area,LEHD_FINAL,"PageNumber",all.x=T)

#Convert NAs to 0s
LEHD_FINAL[is.na(LEHD_FINAL)] <- 0

#proportion totals
LEHD_FINAL$acrepercent<-LEHD_FINAL$Blockinter/LEHD_FINAL$BlockAcres
#LEHD_FINAL$acrepercent[LEHD_FINAL$acrepercent > 1] <- 1
LEHD_FINAL$tot_jobs<-round(LEHD_FINAL$tot_jobs*LEHD_FINAL$acrepercent)
LEHD_FINAL$n11<-round(LEHD_FINAL$n11*LEHD_FINAL$acrepercent)
LEHD_FINAL$n21<-round(LEHD_FINAL$n21*LEHD_FINAL$acrepercent)
LEHD_FINAL$n22<-round(LEHD_FINAL$n22*LEHD_FINAL$acrepercent)
LEHD_FINAL$n23<-round(LEHD_FINAL$n23*LEHD_FINAL$acrepercent)
LEHD_FINAL$n31_33<-round(LEHD_FINAL$n31_33*LEHD_FINAL$acrepercent)
LEHD_FINAL$n42<-round(LEHD_FINAL$n42*LEHD_FINAL$acrepercent)
LEHD_FINAL$n44_45<-round(LEHD_FINAL$n44_45*LEHD_FINAL$acrepercent)
LEHD_FINAL$n48_49<-round(LEHD_FINAL$n48_49*LEHD_FINAL$acrepercent)
LEHD_FINAL$n51<-round(LEHD_FINAL$n51*LEHD_FINAL$acrepercent)
LEHD_FINAL$n52<-round(LEHD_FINAL$n52*LEHD_FINAL$acrepercent)
LEHD_FINAL$n53<-round(LEHD_FINAL$n53*LEHD_FINAL$acrepercent)
LEHD_FINAL$n54<-round(LEHD_FINAL$n54*LEHD_FINAL$acrepercent)
LEHD_FINAL$n55<-round(LEHD_FINAL$n55*LEHD_FINAL$acrepercent)
LEHD_FINAL$n56<-round(LEHD_FINAL$n56*LEHD_FINAL$acrepercent)
LEHD_FINAL$n61<-round(LEHD_FINAL$n61*LEHD_FINAL$acrepercent)
LEHD_FINAL$n62<-round(LEHD_FINAL$n62*LEHD_FINAL$acrepercent)
LEHD_FINAL$n71<-round(LEHD_FINAL$n71*LEHD_FINAL$acrepercent)
LEHD_FINAL$n72<-round(LEHD_FINAL$n71*LEHD_FINAL$acrepercent)
LEHD_FINAL$n81<-round(LEHD_FINAL$n81*LEHD_FINAL$acrepercent)
LEHD_FINAL$n92<-round(LEHD_FINAL$n92*LEHD_FINAL$acrepercent)


#Agg varaibles
LEHD_FINAL.agg<-ddply(LEHD_FINAL,~PageNumber,summarise, 
                      tot_jobs=sum(tot_jobs),
                      n11=sum(n11),
                      n21=sum(n21),
                      n22=sum(n22),
                      n23=sum(n23),
                      n31_33=sum(n31_33),
                      n42=sum(n42),
                      n44_45=sum(n44_45),
                      n48_49=sum(n48_49),
                      n51=sum(n51),
                      n52=sum(n52),
                      n53=sum(n53),
                      n54=sum(n54),
                      n55=sum(n55),
                      n56=sum(n56),
                      n61=sum(n61),
                      n62=sum(n62),
                      n71=sum(n71),
                      n72=sum(n72),
                      n81=sum(n81),
                      n92=sum(n92),
                      CellAcres=mean(CellAcres))


##Classify jobs
office_3<-as.data.frame(LEHD_FINAL.agg$n51+LEHD_FINAL.agg$n52+LEHD_FINAL.agg$n53+LEHD_FINAL.agg$n55+LEHD_FINAL.agg$n92)
colnames(office_3)<-"office_3"
retail_3<-as.data.frame(LEHD_FINAL.agg$n44_45+LEHD_FINAL.agg$n72)
colnames(retail_3)<-"retail_3"
service_3<-as.data.frame(LEHD_FINAL.agg$n54+LEHD_FINAL.agg$n56+LEHD_FINAL.agg$n61+LEHD_FINAL.agg$n62+LEHD_FINAL.agg$n71+LEHD_FINAL.agg$n81)
colnames(service_3)<-"service_3"

office_5<-as.data.frame(LEHD_FINAL.agg$n51+LEHD_FINAL.agg$n52+LEHD_FINAL.agg$n53+LEHD_FINAL.agg$n92)
colnames(office_5)<-"office_5"
retail_5<-as.data.frame(LEHD_FINAL.agg$n44_45)
colnames(retail_5)<-"retail_5"
service_5<-as.data.frame(LEHD_FINAL.agg$n54+LEHD_FINAL.agg$n56+LEHD_FINAL.agg$n61+LEHD_FINAL.agg$n62+LEHD_FINAL.agg$n81)
colnames(service_5)<-"service_5"
industrial_5<-as.data.frame(LEHD_FINAL.agg$n11+LEHD_FINAL.agg$n21+LEHD_FINAL.agg$n22+LEHD_FINAL.agg$n23+LEHD_FINAL.agg$n31_33+LEHD_FINAL.agg$n42+LEHD_FINAL.agg$n48_49)
colnames(industrial_5)<-"industrial_5"
entertainment_5<-as.data.frame(LEHD_FINAL.agg$n71+LEHD_FINAL.agg$n72)
colnames(entertainment_5)<-"entertainment_5"

office_8<-as.data.frame(LEHD_FINAL.agg$n51+LEHD_FINAL.agg$n52+LEHD_FINAL.agg$n53+LEHD_FINAL.agg$n55)
colnames(office_8)<-"office_8"
retail_8<-as.data.frame(LEHD_FINAL.agg$n44_45)
colnames(retail_8)<-"retail_8"
service_8<-as.data.frame(LEHD_FINAL.agg$n54+LEHD_FINAL.agg$n56+LEHD_FINAL.agg$n81)
colnames(service_8)<-"service_8"
industrial_8<-as.data.frame(LEHD_FINAL.agg$n11+LEHD_FINAL.agg$n21+LEHD_FINAL.agg$n22+LEHD_FINAL.agg$n23+LEHD_FINAL.agg$n31_33+LEHD_FINAL.agg$n42+LEHD_FINAL.agg$n48_49)
colnames(industrial_8)<-"industrial_8"
entertainment_8<-as.data.frame(LEHD_FINAL.agg$n71+LEHD_FINAL.agg$n72)
colnames(entertainment_8)<-"entertainment_8"
education_8<-as.data.frame(LEHD_FINAL.agg$n61)
colnames(education_8)<-"education_8"
healthcare_8<-as.data.frame(LEHD_FINAL.agg$n62)
colnames(healthcare_8)<-"healthcare_8"
pubadmin_8<-as.data.frame(LEHD_FINAL.agg$n92)
colnames(pubadmin_8)<-"pubadmin_8"

##prep entopy
office_3e<-as.data.frame(office_3/sum(office_3)*log(office_3/sum(office_3)))
office_3e[is.na(office_3e)] <- 0
retail_3e<-as.data.frame(retail_3/sum(retail_3)*log(retail_3/sum(retail_3)))
retail_3e[is.na(retail_3e)] <- 0
service_3e<-as.data.frame(service_3/sum(service_3)*log(service_3/sum(service_3)))
service_3e[is.na(service_3e)] <- 0

office_5e<-as.data.frame(office_5/sum(office_5)*log(office_5/sum(office_5)))
office_5e[is.na(office_5e)] <- 0
retail_5e<-as.data.frame(retail_5/sum(retail_5)*log(retail_5/sum(retail_5)))
retail_5e[is.na(retail_5e)] <- 0
service_5e<-as.data.frame(service_5/sum(service_5)*log(service_5/sum(service_5)))
service_5e[is.na(service_5e)] <- 0
industrial_5e<-as.data.frame(industrial_5/sum(industrial_5)*log(industrial_5/sum(industrial_5)))
industrial_5e[is.na(industrial_5e)] <- 0  
entertainment_5e<-as.data.frame(entertainment_5/sum(entertainment_5)*log(entertainment_5/sum(entertainment_5)))
entertainment_5e[is.na(entertainment_5e)] <- 0 

office_8e<-as.data.frame(office_8/sum(office_8)*log(office_8/sum(office_8)))
office_8e[is.na(office_8e)] <- 0
retail_8e<-as.data.frame(retail_8/sum(retail_8)*log(retail_8/sum(retail_8)))
retail_8e[is.na(retail_8e)] <- 0
service_8e<-as.data.frame(service_8/sum(service_8)*log(service_8/sum(service_8)))
service_8e[is.na(service_8e)] <- 0
industrial_8e<-as.data.frame(industrial_8/sum(industrial_8)*log(industrial_8/sum(industrial_8)))
industrial_8e[is.na(industrial_8e)] <- 0  
entertainment_8e<-as.data.frame(entertainment_8/sum(entertainment_8)*log(entertainment_8/sum(entertainment_8)))
entertainment_8e[is.na(entertainment_8e)] <- 0 
education_8e<-as.data.frame(education_8/sum(education_8)*log(education_8/sum(education_8)))
education_8e[is.na(education_8e)] <- 0 
healthcare_8e<-as.data.frame(healthcare_8/sum(healthcare_8)*log(healthcare_8/sum(healthcare_8)))
healthcare_8e[is.na(healthcare_8e)] <- 0 
pubadmin_8e<-as.data.frame(pubadmin_8/sum(pubadmin_8)*log(pubadmin_8/sum(pubadmin_8)))
pubadmin_8e[is.na(pubadmin_8e)] <- 0 

#Calc Employment variables
emp.temp<-LEHD_FINAL.agg[1]
emp.temp$EmpDens<-LEHD_FINAL.agg$tot_jobs/LEHD_FINAL.agg$CellAcres
emp.temp$EmpRetDens<-retail_5$retail_5/LEHD_FINAL.agg$CellAcres
emp.temp$Emp<-LEHD_FINAL.agg$tot_jobs
emp.temp$EmpRet<-retail_5$retail_5

#Retail, office and serice jobs entopy
#emp.temp$Entr_3Classes<-(-1)*(LEHD_FINAL.agg$retail_jobs/sum(LEHD_FINAL.agg$retail_jobs)*log(LEHD_FINAL.agg$retail_jobs/sum(LEHD_FINAL.agg$retail_jobs))+
#                              LEHD_FINAL.agg$office_jobs/sum(LEHD_FINAL.agg$office_jobs)*log(LEHD_FINAL.agg$office_jobs/sum(LEHD_FINAL.agg$office_jobs))+
#                              LEHD_FINAL.agg$service_jobs/sum(LEHD_FINAL.agg$office_jobs)*log(LEHD_FINAL.agg$service_jobs/sum(LEHD_FINAL.agg$office_jobs))/log(3))


# #Retail, office and serice jobs entopy
# ent.temp<-cbind(ret,off,srv,agr,min,ind)
# colnames(ent.temp) <- c("ret","off","srv","agr","min","ind")  

#3-class entropy
emp.temp$Entr_3C<-(-1)*((retail_3e$retail_3+
                           office_3e$office_3+
                           service_3e$service_3)/log(3))
#5-class entropy
emp.temp$Entr_5C<-(-1)*((retail_5e$retail_5+
                           office_5e$office_5+
                           service_5e$service_5+
                           industrial_5e$industrial_5+
                           entertainment_5e$entertainment_5)/log(5))
#8-class entropy
emp.temp$Entr_8C<-(-1)*((retail_8e$retail_8+
                           office_8e$office_8+
                           service_8e$service_8+
                           industrial_8e$industrial_8+
                           entertainment_8e$entertainment_8+
                           education_8e$education_8+
                           healthcare_8e$healthcare_8+
                           pubadmin_8e$pubadmin_8)/log(8))



###Get census Data
CENSUS_FINAL<- data.frame(pop=numeric(),
                          housing=numeric(),
                          lat=numeric(), 
                          lon=numeric(),
                          PageNumber=integer(),
                          BlockAcres=numeric(),
                          Blockinter=numeric(),
                          stringsAsFactors=FALSE) 
time_census<-system.time({
  #Chunk grid cells 
  Splitrows=100
  splitnumb_raw<-(nrow(all_cells)/Splitrows)
  if(nrow(all_cells)<Splitrows){
    Splitrows<-nrow(all_cells)
    splitnumb_l<-1
    leftover<-0
  }else{
    splitnumb_l<-floor(nrow(all_cells)/Splitrows)
    leftover<-splitnumb_raw-floor(splitnumb_l)
  }
  
  #Loops through cell chunks
  for(i in 1:splitnumb_l) {
    if(i==splitnumb_l){
      X1=(i-1)*Splitrows
      X2=Splitrows*i+Splitrows*leftover
    }else{i
      X1=(i-1)*Splitrows
      X2=Splitrows*i
    }
    
    #ith grid cell (row) bounding box and corners
    all_cells <-spTransform(all_cells, CRS("+proj=longlat +datum=WGS84"))
    all_cells_cut<-all_cells[X1:X2,]
    bb<-bbox(all_cells_cut)
    
    ##Find intersections & Cul da sacks
    select<-"SELECT 
    pop10                   AS  pop,
    housing10 							AS	housing,
    ST_AsText(geom) AS geom
    FROM blocks_census2010_us
    WHERE  
    geom && 
    ST_MakeEnvelope("
    
    bb2<-paste(bb[1,1],",",bb[2,1],",",bb[1,2],",",bb[2,2],",","4283")
    q_tsig <- paste(select,bb2,")")
    
    #Pull SQL result
    tsig_result <- dbGetQuery(con, q_tsig)
    tsig_result$ID<-row.names.data.frame(tsig_result)
    
    
    #tsig_result1<-tsig_result
    poly.spdf<-WKT2SpatialPolygonsDataFrame(tsig_result, 'geom','ID')
    #plot(test.spdf)
    
    # #Convert node XY to spatial data frame
    # node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(4,3)], data = tsig_result,
    #                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    # 
    #make point coordinate system same as cells
    proj4string(poly.spdf) <- CRS("+proj=longlat +ellps=WGS84") 
    proj4string(all_cells_cut) <- CRS("+proj=longlat +ellps=WGS84") 
    poly.spdf.proj <- spTransform( poly.spdf, CRS( "+init=epsg:3347" ) ) #for area
    
    #Use just the PageNumber column
    all_cells_cut_red<-all_cells_cut[,1]
    
    #spatial join nodes with cells
    #PageNumber = over(poly.spdf,all_cells_cut_red)
    
    X1<-raster::intersect(poly.spdf, all_cells_cut_red)
    PageNumber<-as.data.frame(X1$PageNumber)
    ID<-as.data.frame(X1$ID)
    Res<-cbind(PageNumber,ID)
    colnames(Res)<-c("PageNumber","ID")
    Res<-unique(Res)
    
    all_block_area<-data.frame(round(gArea(poly.spdf.proj, byid = TRUE)*0.000247105,2))
    ABA_ID<-as.data.frame(poly.spdf.proj$ID)
    all_block_area<-cbind(all_block_area,ABA_ID)
    colnames(all_block_area) <- c("BlockAcres","ID")
    
    X1.proj <- spTransform( X1, CRS( "+init=epsg:3347" ) ) #for area
    
    all_block_intersect<-as.data.frame(round(gArea(X1.proj, byid = TRUE)*0.000247105,2))
    all_block_intersect<-cbind(all_block_intersect,ID)
    colnames(all_block_intersect) <- c("Blockinter","ID")
    
    blockprop<-merge(all_block_intersect,all_block_area,"ID",all.x=T)
    
    tsig_result2<-merge(blockprop,tsig_result,"ID",all.x=T)
    tsig_result3<-cbind(tsig_result2[c("pop","housing","Blockinter", "BlockAcres")],Res[1])
    CENSUS_FINAL=rbind(CENSUS_FINAL,tsig_result3)
    
  }
})

print(time_census)


#get cell area
all_cells_area<-data.frame(round(gArea(all_cells.proj, byid = TRUE)*0.000247105,2))
all_cells_area$PageNumber<- 1:nrow(all_cells_area)
colnames(all_cells_area)[1] <- "CellAcres"

CENSUS_FINAL<-merge(CENSUS_FINAL,all_cells_area,"PageNumber")

#Convert NAs to 0s
CENSUS_FINAL[is.na(CENSUS_FINAL)] <- 0

#proportion totals
CENSUS_FINAL$acrepercent<-CENSUS_FINAL$Blockinter/CENSUS_FINAL$BlockAcres
#CENSUS_FINAL$acrepercent[CENSUS_FINAL$acrepercent > 1] <- 1
CENSUS_FINAL$pop<-round(CENSUS_FINAL$pop*CENSUS_FINAL$acrepercent)
CENSUS_FINAL$housing<-round(CENSUS_FINAL$housing*CENSUS_FINAL$acrepercent)
#Agg varaibles
CENSUS_FINAL.agg<-ddply(CENSUS_FINAL,~PageNumber,summarise, 
                        pop=sum(pop),
                        housing=sum(housing),
                        CellAcres=mean(CellAcres))

blocks_count<-as.data.frame(table(CENSUS_FINAL$PageNumber))
colnames(blocks_count)<-c("PageNumber","Blocks")

#Calc census variables
cen.temp<-CENSUS_FINAL.agg[1]
cen.temp$PopDns<-CENSUS_FINAL.agg$pop/CENSUS_FINAL.agg$CellAcres
cen.temp$HUDns<-CENSUS_FINAL.agg$housing/CENSUS_FINAL.agg$CellAcres
cen.temp$Pop<-CENSUS_FINAL.agg$pop
cen.temp$HU<-CENSUS_FINAL.agg$housing
cen.temp$Blocks<-blocks_count$Blocks
cen.temp$AreaAcre<-CENSUS_FINAL.agg$CellAcres

##MERGE final data sets
BASEcells<-as.data.frame(all_cells.df[1])
BASEcells<-merge(BASEcells,cen.temp,"PageNumber",all.x=T)
BASEcells<-merge(BASEcells,emp.temp,"PageNumber",all.x=T)
BASEcells<-merge(BASEcells,int.temp,"PageNumber",all.x=T)


##clean 2nd to lasttime
BASEcells[is.na(BASEcells)] <- 0

#Calc last BE variables
BASEcells$ActvDns<-(BASEcells$Pop+BASEcells$Emp)/BASEcells$AreaAcre
RSE<-as.data.frame(retail_5$retail_5+service_5$service_5+entertainment_5$entertainment_5)
colnames(RSE)<-"RSE"
BASEcells$RSEDns<-(RSE$RSE)/BASEcells$AreaAcre
pop1<-as.data.frame(BASEcells$Pop)
colnames(pop1)<-"pop1"
pop1$pop1[pop1$pop1==0]<-1
BASEcells$EPB1<-BASEcells$Emp/pop1$pop1 
alpha<-sum(BASEcells$Emp)/sum(pop1$pop1)
BASEcells$EPB1_n<-(BASEcells$Emp/pop1$pop1)/alpha
BASEcells$EPB2<-retail_5$retail_5/pop1$pop1
beta<-sum(retail_5)/sum(pop1$pop1)
BASEcells$EPB2_n<-(retail_5$retail_5/pop1$pop1)/beta
BASEcells$EPB3<-RSE$RSE/pop1$pop1
gamma<-beta<-sum(RSE)/sum(pop1$pop1)
BASEcells$EPB3_n<-(RSE$RSE/pop1$pop1)/(RSE$RSE)

##last clean
BASEcells[is.na(BASEcells)] <- 0
BASEcells<-as.data.frame(sapply(BASEcells,as.numeric))
#BASEcells<-do.call(data.frame,sapply(BASEcells, function(x) replace(x, is.infinite(x),NA)))

#Join BE data to grid cells 
#all_cells@data = data.frame(all_cells@data, BASEcells[match(all_cells@data[,'PageNumber'], BASEcells[,'PageNumber']),])

#} #end county loop if looping

##-------------Calc netwotk distance to transit----------------##
if (file.exists(paste0("temp/",studyarea,".csv"))) file.remove(paste0("temp/",studyarea,".csv"))

TRANSIT_FINAL<- data.frame(
  originID=numeric(), 
  o_x=numeric(), 
  o_y=numeric(), 
  RailDIS_n=numeric(),
  rail_x=numeric(), 
  rail_y=numeric(),
  BusDIS_n=numeric(),
  bus_x=numeric(), 
  bus_y=numeric(),
  stringsAsFactors=FALSE)

#Chunk grid cells 
Splitrows=1000
splitnumb_raw<-(nrow(origins)/Splitrows)
if(nrow(origins)<Splitrows){
  Splitrows<-nrow(origins)
  splitnumb_l<-1
  leftover<-0
}else{
  splitnumb_l<-floor(nrow(origins)/Splitrows)
  leftover<-splitnumb_raw-floor(splitnumb_l)
}

#Loops through cell chunks
for(j in 1:splitnumb_l) {
  
  if(j==splitnumb_l){
    X1=(j-1)*Splitrows
    X2=Splitrows*j+Splitrows*leftover
  }else{j
    X1=(j-1)*Splitrows
    X2=Splitrows*j
  }  
  
  num.cores <- detectCores() ## detect how many cores on your machine
  cl <- makeCluster(num.cores, type = "SOCK")
  registerDoParallel(cl)
  getDoParWorkers() ## check if all cores are working
  clusterEvalQ(cl,library(RPostgreSQL))
  clusterEvalQ(cl,library(DBI))
  clusterEvalQ(cl,library(RSQLite))
  clusterEvalQ(cl,library(sp))
  clusterEvalQ(cl,library(rgeos))
  
  origins_cut<-origins[X1:X2,]
  
  #Read in bike station locations
  timeTR<-system.time({
    OO<-nrow(origins_cut)

    TRANSIT_FINAL=foreach(i=1:OO,.combine = rbind,.errorhandling='remove') %dopar% {

      drv <- dbDriver("PostgreSQL")
      con <- dbConnect(drv, user = "postgres", dbname = "osm", host = "localhost")
      con_rt <- dbConnect(drv, user = "postgres", dbname = "opbeumDB", host = "localhost")
      
      origins_bb<-bbox(origins_cut[i,])
      #get origin ID and YX
      originID<-origins_cut$PageNumber[i] #change this to the ID column for the loaded shape'
      bbl<-origins_bb*.995  
      bbh<-origins_bb*1.005  
      O_xy<-coordinates(origins_cut[i,])
      
      #transform to planar coordinates
      origins.proj <- spTransform( origins_cut[i,], CRS( "+init=epsg:3347" ) ) 
      
      ##Find destinations in area
      select<-"SELECT 
      route_type,
      ST_Y(geom) AS lat, 
      ST_X(geom) AS lon 
      FROM gtfs_stops_us
      WHERE  
      geom && 
      ST_MakeEnvelope("
      
      bb2<-paste(bbl[1,1],",",bbl[2,1],",",bbh[1,2],",",bbh[2,2],",","4283")
      q_tsig <- paste(select,bb2,")")
      
      #Pull SQL result
      tsig_result <- dbGetQuery(con, q_tsig)
      
      #get only bus stations
      tsig_result_bus <- tsig_result[ which(tsig_result$route_type==3), ]
      
      #get only rail stations
      tsig_result_rail <- tsig_result[ which(tsig_result$route_type<3), ]
      
      ##
      ##get ORIGIN OSM NODE
      ##
      
      bbl<-origins_bb*.9995  
      bbh<-origins_bb*1.0005 
      
      ##Find closest OSM node or the Origin
      select<-"SELECT 
      osm_id,
      ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
      ST_X(ST_Transform(geom_vertex,4326)) AS lon 
      FROM at_2po_vertex
      WHERE  
      geom_vertex && 
      ST_MakeEnvelope("
      
      bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
      q_tsig <- paste(select,bb2,")")
      
      #Pull SQL result
      tsig_result <- dbGetQuery(con_rt, q_tsig)
      
      #Convert node XY to spatial data frame
      node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      #transform to planor coordinates
      node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) ) 
      
      #get closest osm node
      near_node<-gDistance(origins.proj, spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
      node_dist<-min(near_node)
      node_id<-node.spdf$osm_id[which(near_node==min(near_node))]
      
      #get the 'source' node for routing
      select<-"SELECT source
      FROM at_2po_4pgr
      WHERE osm_source_id =" 
      q_node_origin <- paste(select,node_id)
      
      qnode_result_origin <- dbGetQuery(con_rt, q_node_origin)[1,]
      
      #********************
      #Calc RAIL diatance**
      #********************
      if(nrow(tsig_result_rail)>0){
        
        transit.spdf_rail <- SpatialPointsDataFrame(coords = tsig_result_rail[,c(3,2)], data = tsig_result_rail,
                                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
        
        #transform to planar coordinates
        transit.spdf_rail.proj <- spTransform( transit.spdf_rail, CRS( "+init=epsg:3347" ) ) 
        
        #get closest rail location
        near_rail<-gDistance(origins.proj, spgeom2=transit.spdf_rail.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
        rail_dist<-min(near_rail)
        rail_xy<-unique(coordinates(transit.spdf_rail[which(near_rail==min(near_rail)),]))
        rail_row<-which(near_rail==min(near_rail))[1]
        
        ##
        ##get RAIL OSM NODE
        ##
        
        bbl<-rail_xy*.9995  
        bbh<-rail_xy*1.0005  
        
        ##Find closest OSM node
        select<-"SELECT 
        osm_id,
        ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
        ST_X(ST_Transform(geom_vertex,4326)) AS lon 
        FROM at_2po_vertex
        WHERE  
        geom_vertex && 
        ST_MakeEnvelope("
        
        bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
        q_tsig <- paste(select,bb2,")")
        
        
        #Pull SQL result
        tsig_result <- dbGetQuery(con_rt, q_tsig)
        
        if(nrow(tsig_result)>0){ 
          
          #Convert node XY to spatial data frame
          node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          #transform to planor coordinates
          node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) ) 
          
          #get closest osm node
          near_node_rail<-gDistance(transit.spdf_rail.proj[rail_row,], spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
          node_rail_dist<-min(near_node_rail)
          node_rail_id<-node.spdf$osm_id[which(near_node_rail==min(near_node_rail))]
          
          #get the 'source' node for routing
          select<-"SELECT target
          FROM at_2po_4pgr
          WHERE osm_target_id =" 
          q_node_rail <- paste(select,node_rail_id)
          
          qnode_result_rail <- dbGetQuery(con_rt, q_node_rail)[1,]
          
          
          #----------------------
          
          ##get OD dist to rail
          
          q1<-   "SELECT seq, id1 AS node, id2 AS edge,km AS cost
          FROM pgr_astar('
          SELECT id AS id,
          source,
          target,
          cost,
          x1, y1, x2, y2
          FROM at_2po_4pgr as r,
          (SELECT ST_Expand(ST_Extent(geom_way),0.01) as box  FROM at_2po_4pgr as l1  
          WHERE l1.source ="
          o_node<-qnode_result_origin
          q2<-"OR l1.target =" 
          d_node<-qnode_result_rail
          q3<-") as box
          WHERE r.geom_way && box.box',"
          q4<-","
          q5<-", false, false)as r INNER JOIN at_2po_4pgr as g ON r.id2 = g.id ;"
          
          q_railD <- paste0(q1,o_node,q2,d_node,q3,o_node,q4,d_node,q5)
          
          d_result_rail <- dbGetQuery(con_rt, q_railD)
          
          RailDIS_n=sum(d_result_rail$cost)*0.621371
          
          #----------------------------  
        } else{
          RailDIS_n<-9999
          rail_xy<-matrix(0, 1, 2)
        }
        
      } else{
        RailDIS_n<-9999
        rail_xy<-matrix(0, 1, 2)
      }
      
      #********************
      #Calc BUS diatance**
      #********************
      if(nrow(tsig_result_bus)>0){
        
        #Convert node XY to spatial data frame
        transit.spdf_bus <- SpatialPointsDataFrame(coords = tsig_result_bus[,c(3,2)], data = tsig_result_bus,
                                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))    
        #transform to planar coordinates
        transit.spdf_bus.proj <- spTransform( transit.spdf_bus, CRS( "+init=epsg:3347" ) ) 
        
        #get closest bus location
        near_bus<-gDistance(origins.proj, spgeom2=transit.spdf_bus.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
        bus_dist<-min(near_bus)
        bus_xy<-unique(coordinates(transit.spdf_bus[which(near_bus==min(near_bus)),]))
        bus_row<-which(near_bus==min(near_bus))[1]
        
        ##
        ##get ORIGIN OSM NODE
        ##
        
        ##Find closest OSM node or the Origin
        select<-"SELECT 
        osm_id,
        ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
        ST_X(ST_Transform(geom_vertex,4326)) AS lon 
        FROM at_2po_vertex
        WHERE  
        geom_vertex && 
        ST_MakeEnvelope("
        
        bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
        q_tsig <- paste(select,bb2,")")
        
        #Pull SQL result
        tsig_result_bus <- dbGetQuery(con_rt, q_tsig)
        
        if(nrow(tsig_result_bus)>0){  
          
          #Convert node XY to spatial data frame
          node.spdf <- SpatialPointsDataFrame(coords = tsig_result_bus[,c(3,2)], data = tsig_result_bus,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          #transform to planor coordinates
          node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) ) 
          
          #get closest osm node
          near_node<-gDistance(origins.proj, spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
          node_dist<-min(near_node)
          node_id<-node.spdf$osm_id[which(near_node==min(near_node))]
          
          
          #----------------------
          
          ##
          ##get bus OSM NODE
          ##
          
          bbl<-bus_xy*.99995  
          bbh<-bus_xy*1.00005  
          
          ##Find closest OSM node
          select<-"SELECT 
          osm_id,
          ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
          ST_X(ST_Transform(geom_vertex,4326)) AS lon 
          FROM at_2po_vertex
          WHERE  
          geom_vertex && 
          ST_MakeEnvelope("
          
          bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
          q_tsig <- paste(select,bb2,")")
          
          #Pull SQL result
          tsig_result <- dbGetQuery(con_rt, q_tsig)
          
          #Convert node XY to spatial data frame
          node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          #transform to planor coordinates
          node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) ) 
          
          #get closest osm node
          near_node_bus<-gDistance(transit.spdf_bus.proj[bus_row,], spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
          node_bus_dist<-min(near_node_bus)
          node_bus_id<-tsig_result$osm_id[which(near_node_bus==min(near_node_bus))]
          
          #get the 'source' node for routing
          select<-"SELECT target
          FROM at_2po_4pgr
          WHERE osm_source_id =" 
          q_node_bus <- paste(select,node_bus_id)
          
          qnode_result_bus <- dbGetQuery(con_rt, q_node_bus)[1,]
          
          #----------------------------  
          
          ##get OD dist to Bus
          
          q1<-   "SELECT seq, id1 AS node, id2 AS edge,km AS cost
          FROM pgr_astar('
          SELECT id AS id,
          source,
          target,
          cost,
          x1, y1, x2, y2
          FROM at_2po_4pgr as r,
          (SELECT ST_Expand(ST_Extent(geom_way),0.01) as box  FROM at_2po_4pgr as l1  
          WHERE l1.source ="
          o_node<-qnode_result_origin
          q2<-"OR l1.target =" 
          d_node<-qnode_result_bus
          q3<-") as box
          WHERE r.geom_way && box.box',"
          q4<-","
          q5<-", false, false)as r INNER JOIN at_2po_4pgr as g ON r.id2 = g.id ;"
          
          q_busD <- paste0(q1,o_node,q2,d_node,q3,o_node,q4,d_node,q5)
          
          d_result_bus <- dbGetQuery(con_rt, q_busD)
          
          BusDIS_n=sum(d_result_bus$cost)*0.621371
          
          #----------------------------  
          
        } else{
          BusDIS_n<-9999
          bus_xy<-matrix(0, 1, 2)
        }
        
      } else{
        BusDIS_n<-9999
        bus_xy<-matrix(0, 1, 2)
      }
      
      #end distance calcs  
      
      dbDisconnect(con)
      dbDisconnect(con_rt)
      
      #prep column names
      colnames(O_xy) <- c("o_x","o_y")
      colnames(rail_xy) <- c("rail_x","rail_y")
      colnames(bus_xy) <- c("bus_x","bus_y")
      
      #Merge all location results
      cbind.data.frame(originID,O_xy,RailDIS_n,rail_xy,BusDIS_n,bus_xy,row.names = NULL)
      
    }  
    
    write.table(TRANSIT_FINAL, paste0("temp/",studyarea,".csv"), row.names = F, col.names = F, append = T, sep=",",quote=F)
    
    stopCluster(cl)
  })
  print(timeTR)
  print(splitnumb_l-j)
}

TRANSIT_FINAL2 = fread(paste0("temp/",studyarea,".csv"), sep=",",header = F)
colnames(TRANSIT_FINAL2)<-c("PageNumber","o_x","o_y","RailDIS_n","rail_x","rail_y","BusDIS_n","bus_x","bus_y")
#write.table(TRANSIT_FINAL2, paste0(studyarea,"2",".csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)

########################
## End Transit router ##
########################

#Merge transit distances with final dataset
BASEcells<-merge(BASEcells,TRANSIT_FINAL2,all.x=T)

#####CREATE STEVE'S NEW VARIABLE!
BASEcells$ACTIVITY<-BASEcells$Pop +BASEcells$Emp


#Join BE data to grid cells 
all_cells@data = data.frame(all_cells@data, BASEcells[match(all_cells@data[,'PageNumber'], BASEcells[,'PageNumber']),])

#Write outputs (cSV and Shapefile)
writeOGR(obj=all_cells, dsn="Output", layer="OPBEUM_RESULT_GRID", driver="ESRI Shapefile",overwrite_layer=T,check_exists=T) 
write.table(BASEcells, paste0("Output/","OPBEUM_RESULT_GRID.csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)

################
## END OPBEUM ##
################

############################
##Add-on for normalization##
############################

##Gather cells in 1/2 mile radius
all_cells.df<-as.data.frame(all_cells)

##Simpple earth distance formula
earth.dist <- function (long1, lat1, long2, lat2)
{
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

##Get unique combinations,convert to a dataframe, merge with XY data
all_cells.df_red<-all_cells.df[c("PageNumber","o_x","o_y")]
all_cells.df_red_comb<-data.frame((do.call(rbind,combn(all_cells.df_red$PageNumber,2,simplify = FALSE))))
names(all_cells.df_red_comb) <- c("PageNumber_O","PageNumber_D")
all_cells.df_red_comb<-merge(all_cells.df_red_comb,all_cells.df_red,by.x=c("PageNumber_O"),by.y=c("PageNumber"))
  names(all_cells.df_red_comb) <- c("PageNumber_O","PageNumber_D","o_x","o_y")

## get straight line distance b/t cells
all_cells.df$earth_dist <- sapply(1:nrow(all_cells.df),function(i)
  earth.dist(all_cells.df$o_y[i],all_cells.df$o_x[i],all_cells.df$o_y[i],all_cells.df$o_y[i]))

##reduce dataset to .5 miles per pair
all_cells.df$earth_dist<-all_cells.df$earth_dist*0.621371
sInDF_SUB_Comb_OD_red<-all_cells.df[ which(all_cells.df$earth_dist<=.5), ]

##Array raw cell data


####---------------------######################---------------------####
####---------------------##   ADD TYPOLOGY   ##---------------------####
####---------------------######################---------------------####
####################################################
# Steven R. Gehrke
# TRB 2017: Bikeshare Typology
# 07.16.2016
####################################################

####################################################
# Define Directories & Load Libraries
####################################################

# Install packages and load libraries
#install.packages("mclust")
library("mclust")
#citation("mclust")
#install.packages("pastecs")
library("pastecs")

####################################################
# Read, Prepare, and Describe Data
####################################################

# Read and subset data
dat <- read.csv(paste("Input","cabi_opbeum_20160716_v2.csv",sep="/"),header=TRUE)
head(dat)
colnames(dat)
cabi <- dat[,10:42]
head(cabi)
#cabi_unadj <- read.csv(paste(DirData,"cabi_opbeum_20160620_centroid.csv",sep="/"),header=TRUE)
#head(cabi_unadj)

# Examine descriptive statistics
#stat.desc(cabi)
#stat.desc(cabi_unadj)

# Examne correlation matrix
cabix <- cabi[,2:33]
cor(cabix, use="complete.obs", method="kendall")

# Reduce datasets based on unadjusted correlations
cabi2 <- subset(cabi,select=c("ID_CABI","ADJPOP","ADJEMP","ADJACT","ADJEMPRET",
                              "ADJENT5C","ADJEPB1N","ADJEPB2N","ADJEPB3N",
                              "ADJBLKS","ADJINTER","ADJCULDS","ADJLINKS",
                              "ADJCNR","ADJBETA","DISTR","DISTB","NUMBR","NUMBB"))
cabix2 <- cabi2[,2:19]
cor(cabix2, use="complete.obs", method="kendall")

####################################################
# Additional Variable Construction
####################################################

# Retail proportion of total employment
cabi$ADJEMPRETP <- cabi$ADJEMPRET/cabi$ADJEMP
summary(cabi$ADJEMPRETP)

# Presence of WMATA metro station
cabi$PRESR <- ifelse(cabi$NUMBR>=1,1,0)
table(cabi$PRESR)

####################################################
# Latent Class Cluster Analysis
####################################################

### Step 1: Decide on the built environment measures
colnames(cabi)
cabi3 <- subset(cabi, select=c("ADJHOU","ADJACT","ADJEPB1","ADJBETA","DISTR"))
names(cabi3)<-c("HU","ACTIVITY","EPB1","BETA","RailDIS_n")
plot(cabi3, pch=20, cex=0.3)

### Step 2: Estimate latent class cluster analysis model
#?mclust

# Estimate Bayesian information criterion (BIC)
m1bic <- mclustBIC(cabi3)
#m1bic <- mclustBIC(cabi3, G=1:20, x=m1bic)
# summary(m1bic)
# plot(m1bic, G=1:10)

# Estimate latent class cluster analysis (LCCA) model
m1 <- Mclust(cabi3, x=m1bic)
# summary(m1, parameters=TRUE)
# plot(m1, what = "classification")
# table(m1$classification)

# Estimate integrated classification likelihood (ICL) criterion
m1icl <- mclustICL(cabi3)
# summary(m1icl) # same number of classes as BIC ???
# plot(m1icl)

# Perform Lo-Mendell-Rubin test (LRT)
#m1lrt <- mclustBootstrapLRT(cabi3, modelName = "VVV") # long run time
#m1lrt

### Step 3: Examine LCCA model results

# Class assignment and relative probability of class assignment
dat$CLUSTER <- m1$classification
cabix3 <- round(m1$z,5)
cabix3 <- as.data.frame(as.matrix(cabix3))
colnames(cabix3) <- c("ZSCORE1","ZSCORE2","ZSCORE3","ZSCORE4","ZSCORE5")
dat <- cbind(dat,cabix3)
head(dat)

# Calculate average probability of alternative cluster classifications
cl1 <- dat[dat$CLUSTER==1,]
cl2 <- dat[dat$CLUSTER==2,]
cl3 <- dat[dat$CLUSTER==3,]
cl4 <- dat[dat$CLUSTER==4,]
cl5 <- dat[dat$CLUSTER==5,]
colnames(dat)
cl1x <- round(colMeans(cl1[,44:48]),2)
cl2x <- round(colMeans(cl2[,44:48]),2)
cl3x <- round(colMeans(cl3[,44:48]),2)
cl4x <- round(colMeans(cl4[,44:48]),2)
cl5x <- round(colMeans(cl5[,44:48]),2)

m1tab <- do.call("rbind", list(cl1x, cl2x, cl3x, cl4x, cl5x))
m1tab <- as.data.frame(as.matrix(m1tab))
rm(cabix3,cl1,cl2,cl3,cl4,cl5,cl1x,cl2x,cl3x,cl4x,cl5x)
print(m1tab)

####################################################
# Export Dataset
####################################################

# Export dataset for use in BIKELOC proof of concept
head(dat)
write.csv(dat, file=paste("output","cabi_opbeum_lcca_20160718.csv",sep="/"))

####################################################
# End of Script
####################################################



###----------------APPLY TO REGION -----------------
#BASEcells = data.frame(fread(paste0("Output/","OPBEUM_RESULT_GRID.csv"), sep=",",header = TRUE))
  BASEcells_LCAAPPLY_1<-BASEcells[complete.cases(BASEcells),]    
  BASEcells_LCAAPPLY_2<-subset(BASEcells_LCAAPPLY_1, select=c("HU","ACTIVITY","EPB1","BETA","RailDIS_n"))

    
BASEcells_LCA<-predict(m1, BASEcells_LCAAPPLY_2,x=m1bic)

# Class assignment and relative probability of class assignment
BASEcells_LCAAPPLY_2$CLUSTER <- BASEcells_LCA$classification
LCA <- round(BASEcells_LCA$z,5)
LCA <- as.data.frame(as.matrix(LCA))
colnames(LCA) <- c("ZSCORE1","ZSCORE2","ZSCORE3","ZSCORE4","ZSCORE5")
BASEcells_LCAAPPLY_2 <- cbind(BASEcells_LCAAPPLY_2,LCA)

#Subset to results only
BASEcells_LCAAPPLY_2<-BASEcells_LCAAPPLY_2[c("CLUSTER","ZSCORE1","ZSCORE2","ZSCORE3","ZSCORE4","ZSCORE5")]

BASEcells_LCAAPPLY_3<-cbind(BASEcells_LCAAPPLY_1[c(1)],BASEcells_LCAAPPLY_2)

##bring back incomplete cases
BASEcells_LACRESULT<-merge(BASEcells[1],BASEcells_LCAAPPLY_3,"PageNumber",all.x=T)
BASEcells<-merge(BASEcells[1],BASEcells_LCAAPPLY_3,"PageNumber",all.x=T)

#add to shapefile
#all_cells<-readOGR(dsn = "Output", layer = paste0("OPBEUM_RESULT_GRID"))
#names(all_cells)[1]<-"PageNumber"
###clean out data
#all_cells<-all_cells[1]

all_cells@data = data.frame(all_cells@data, BASEcells_LACRESULT[match(all_cells@data[,'PageNumber'], BASEcells_LACRESULT[,'PageNumber']),])
    all_cells.proj<-all_cells
    all_cells.df<-data.frame(all_cells)
    # origins<- gCentroid(all_cells,byid=TRUE)
    # origins <- SpatialPointsDataFrame(origins,all_cells.df)
    # colnames(origins@data)[1] <- "PageNumber"
    # 
####---------------------######################---------------------####
####---------------------##   ADD BIKE LOC   ##---------------------####
####---------------------######################---------------------####

#make point coordinate system same as cells
proj4string(bs.stations) <- CRS("+proj=longlat +ellps=WGS84") 
proj4string(all_cells) <- CRS("+proj=longlat +ellps=WGS84") 

all_cells_PG<-all_cells.proj[1]

BSjoin = over(bs.stations,all_cells[1])
BSjoin<-cbind(data.frame(bs.stations),BSjoin)

##Clean up join with dummy for B/S presence amd dup/rename ID to TERMINAL_N to match OD data
BSjoin$BScell<-BSjoin$PageNumber/BSjoin$PageNumber
#BSjoin<-BSjoin[complete.cases(BSjoin),]
#  BSjoin <- BSjoin %>% mutate( TERMINAL_N = ID )     

#Merge with master BE data
all_cells_bs<-all_cells
all_cells_bs@data = data.frame(all_cells_bs@data, BSjoin[match(all_cells_bs@data[,'PageNumber'], BSjoin[,'PageNumber']),])
all_cells_bs$BScell[is.na(all_cells_bs$BScell)] <- 0

# writeOGR(obj=all_cells_bs, dsn="CaBi_Counties", layer=paste0("OPBEUM_RESULT_GRID_",z), driver="ESRI Shapefile",overwrite_layer=T,check_exists=T) #use for loop
# #writeOGR(obj=all_cells_bs, dsn="", layer="OPBEUM_RESULT_GRID", driver="ESRI Shapefile",overwrite_layer=T,check_exists=T)             #use for specific area
# 
# write.table(all_cells_bs, paste0("CaBi_Counties/OPBEUM_RESULT_GRID.csv",z,".csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)
# #write.table(all_cells_bs, "OPBEUM_RESULT_GRID.csv", row.names = F, col.names = T, append = F, sep=",",quote=F)

######################
##   SVM Location   ##
######################
##Read Base Data & Join XYs
#BASEcells = fread("OPBEUM_RESULT_GRID.csv", sep=",",header = TRUE)
BASEcells<-data.frame(all_cells_bs)
    BASEcells[is.na(BASEcells)]<-0


# CellXY_ALL<-fread("GRID_BIKESHARE_YXPoints.csv", sep=",",header = TRUE)
# CellXY_ALL<-subset(CellXY_ALL, select = c(3:5))
# names(CellXY_ALL)[names(CellXY_ALL)=="PAGENUMBER"] <- "PageNumber"
# names(CellXY_ALL)[names(CellXY_ALL)=="X"] <- "XCOORD"
# names(CellXY_ALL)[names(CellXY_ALL)=="Y"] <- "YCOORD"
# 
# BASEcellsReduced<-merge(BASEcellsReduced,CellXY_ALL,by="PageNumber")
# 
# ##Read and Join Buffer with non-BS cells
# BS_BufferedCells = fread("GRID_BIKESHARE_BUFFER.csv", sep=",",header = TRUE)
# names(BS_BufferedCells)[names(BS_BufferedCells)=="PAGENUMBER"] <- "PageNumber"
# BS_BufferedCells$BSAREA<-1
# BS_BufferedCells_DATA<-merge(BS_BufferedCells,BASEcells,by="PageNumber")
# 
# BS_BufferedCells_APPLY<-merge(BASEcellsReduced,BS_BufferedCells,by="PageNumber", all = TRUE)
# BS_BufferedCells_APPLY$BSAREA[is.na(BS_BufferedCells_APPLY$BSAREA)]<-0
# BS_BufferedCells_APPLY<-subset(BS_BufferedCells_APPLY,BS_BufferedCells_APPLY$BSAREA==0)
  
BASEcells_sample.a<-BASEcells[which(BASEcells$BScell==1),]
BASEcells_sample.b<-BASEcells[sample(nrow(BASEcells), nrow(BASEcells)*.05), ]
BASEcells_sample<-rbind(BASEcells_sample.a,BASEcells_sample.b)

svmfit <- svm(factor(BScell)~
                ZSCORE1+ ZSCORE2+ ZSCORE3+ ZSCORE4+ ZSCORE5,
              data=BASEcells_sample,
              cost=10^10)

# print(svmfit)
# summary(svmfit)

# test with training data
pred <- predict(svmfit, BASEcells_sample)

# Check accuracy:
table(pred, BASEcells_sample$BScell)

###Now predict for the entire study area
BASEcells_NOBS<-BASEcells[which(BASEcells$BScell==0),]

BASEcells_NOBS$pred2 <- predict(svmfit, newdata=BASEcells_NOBS) 
summary(BASEcells_NOBS$pred2)

BASEcells_Result <- BASEcells_NOBS[ which(BASEcells_NOBS$pred2==1), ]
#BASEcells_Result <- subset(BASEcells_Result, select = c("PageNumber","pred2","o_x","o_y"))

##Write final Basic SVM result
#write.table(BASEcells_Result,paste0("Output/","SVM_Result_Class.csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)

######################
## Cluster Analysis ##
######################
origins_accept<-merge(origins,BASEcells_Result[1],all.x=F)
# sp1 <- sp::SpatialPoints(coords = cbind(BASEcells_Result$X,BASEcells_Result$Y))
#spdf <- SpatialPointsDataFrame(origins_accept, BASEcells_Result)

#lnd<-readShapeSpatial("GRID_BIKESHARE_XYPoints.shp")
#class(spdf)
#plot(spdf)
#dev.copy2pdf(file="1_points.pdf")

sSp <- as(SpatialPoints(origins_accept), "ppp")  # convert points to pp class
Dens <- density(sSp, adjust = 0.15)  # create density object
plot(Dens)  # default plot for density
dev.copy2pdf(file="2_Density.pdf")

contour<-contour(density(sSp, adjust = 0.2), nlevels = 10)  # plot as contours - this is where we're heading
plot(contour)
dev.copy2pdf(file="3_Countour.pdf")

Dsg <- as(Dens, "SpatialGridDataFrame")  # convert to spatial grid class
Dim <- as.image.SpatialGridDataFrame(Dsg)  # convert again to an image
Dcl <- contourLines(Dim, nlevels = 7)  # create contour object - change 8 for more/fewer levels
SLDF <- ContourLines2SLDF(Dcl, proj4string=CRS(as.character(NA)))  # convert to SpatialLinesDataFrame
plot(SLDF, col = terrain.colors(5))
dev.copy2pdf(file="4_Terrain.pdf")

Polyclust <- rgeos::gPolygonize(SLDF[4, ])
gas <- rgeos::gArea(Polyclust, byid = T)/10000
Polyclust <- SpatialPolygonsDataFrame(Polyclust, data = data.frame(gas), match.ID = F)
plot(Polyclust)
dev.copy2pdf(file="5_PolyCluster.pdf")

proj4string(Polyclust) <- CRS("+proj=longlat +ellps=WGS84") 
proj4string(origins_accept) <- CRS("+proj=longlat +ellps=WGS84") 

cAg <- aggregate(origins_accept, by = Polyclust, FUN = length)
#lb <- gBoundary(lnd)

plot(Dens, main = "")
#plot(lnd, border = "grey", lwd = 2, add = T)
plot(SLDF, col = terrain.colors(8), add = T)
plot(cAg, col = "red", border = "white", add = T)
graphics::text(coordinates(cAg) + 1000, labels = cAg$CODE)
dev.copy2pdf(file="6_Combined.pdf")


sIn <- origins_accept[cAg, ]  # select the cells inside the clusters
sOut <- bs.stations[!row.names(bs.stations) %in% row.names(sIn), ]  # stations outside the clusters
plot(sIn)  # the more sparsely distributed points - notice the 'holes' of low density
dev.copy2pdf(file="7_CellsInCluster.pdf")
plot(cAg, border = "red", lwd = 3, add = T)
dev.copy2pdf(file="8_CellsInCluster_Boundary.pdf")

nrow(sIn)/nrow(origins_accept)  # proportion of points in cluster

##Write Clustered cells result
write.table(sIn, "Cell_Cluster_Results.csv", row.names = F, col.names = T, append = F, sep=",",quote=F)

######################
##  SVM Regression  ##
######################
BASEcells_Reduced_Acceptable<-BASEcells_Result
BASEcells_Reduced_Acceptable<-BASEcells_Reduced_Acceptable[which(BASEcells_Reduced_Acceptable$ActvDns>0), ]

##Read in location of exisitng stations
# GRID_BIKESHARE_Location = fread("GRID_BIKESHARE_Location.csv", sep=",",header = TRUE)
# names(GRID_BIKESHARE_Location)[names(GRID_BIKESHARE_Location)=="PAGENUMBER"] <- "PageNumber"


##Aggregate Real ODs to Grid (all ovservations)
GRID_BIKESHARE_REALODs_AGG <- plyr::count(GRID_BIKESHARE_REALODs, c("TERMINAL_N_O","TERMINAL_N_D"))
    names(GRID_BIKESHARE_REALODs_AGG)[names(GRID_BIKESHARE_REALODs_AGG)=="freq"] <- "TRIPS"

    test2013<-nrow(GRID_BIKESHARE_REALODs[which(grepl("/2013",GRID_BIKESHARE_REALODs$Start.date)=="TRUE"),])
    test2012<-nrow(GRID_BIKESHARE_REALODs[which(grepl("/2012",GRID_BIKESHARE_REALODs$Start.date)=="TRUE"),])
    test2011<-nrow(GRID_BIKESHARE_REALODs[which(grepl("/2011",GRID_BIKESHARE_REALODs$Start.date)=="TRUE"),])
    test2010<-nrow(GRID_BIKESHARE_REALODs[which(grepl("/2010",GRID_BIKESHARE_REALODs$Start.date)=="TRUE"),])
    
    
##Join FULL GRID data with Bike Station Info
#GRIDCELLS_ALL = fread("MD_Gridcell_Dataset_2016_04_21.csv", sep=",",header = TRUE)

  GRID_BIKESHARE_Location_GRIDInfo<-BASEcells

GRID_BIKESHARE_Location_GRIDInfo_RED_O<-GRID_BIKESHARE_Location_GRIDInfo
  names(GRID_BIKESHARE_Location_GRIDInfo_RED_O)<-paste(names(GRID_BIKESHARE_Location_GRIDInfo_RED_O),"O", sep = "_")
GRID_BIKESHARE_Location_GRIDInfo_RED_D<-GRID_BIKESHARE_Location_GRIDInfo
                  names(GRID_BIKESHARE_Location_GRIDInfo_RED_D)<-paste(names(GRID_BIKESHARE_Location_GRIDInfo_RED_D),"D", sep = "_")

#Clean up data
GRID_BIKESHARE_REALODs_AGG<-GRID_BIKESHARE_REALODs_AGG[which(GRID_BIKESHARE_REALODs_AGG$TERMINAL_N_O>0),]
GRID_BIKESHARE_REALODs_AGG<-GRID_BIKESHARE_REALODs_AGG[which(GRID_BIKESHARE_REALODs_AGG$TERMINAL_N_D>0),]
 
##Join Locations with Real Aggregated ODs - Origins
GRID_BIKESHARE_REALODs_GRIDInfo_O<-merge(GRID_BIKESHARE_REALODs_AGG,GRID_BIKESHARE_Location_GRIDInfo_RED_O,by="TERMINAL_N_O")
##Join Locations with Real Aggregated ODs - Destinations
GRID_BIKESHARE_REALODs_GRIDInfo_OD<-merge(GRID_BIKESHARE_REALODs_GRIDInfo_O,GRID_BIKESHARE_Location_GRIDInfo_RED_D,by="TERMINAL_N_D")

###Get Real OD distances
GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY <- subset(GRID_BIKESHARE_REALODs_GRIDInfo_OD, select = c("PageNumber_O","PageNumber_D","LONGITUDE_O", "LATITUDE_O","LONGITUDE_D", "LATITUDE_D"))
    colnames(GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY)<-c("PageNumber_O", "PageNumber_D", "bs_x_O", "bs_y_O", "bs_x_D", "bs_y_D")

    ##-------------Calc netwotk distance to transit----------------##
    if (file.exists(paste0("temp/",studyarea,"_BSOD.csv"))) file.remove(paste0("temp/",studyarea,"_BSOD.csv"))
    
    
    REALOD_FINAL<- data.frame(
      PageNumber_O=numeric(), 
      PageNumber_D=numeric(), 
      ODDIS_n=numeric(),
      stringsAsFactors=FALSE)
    
    #Chunk grid cells 
    Splitrows=1000
    splitnumb_raw<-(nrow(GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY)/Splitrows)
    if(nrow(GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY)<Splitrows){
      Splitrows<-nrow(GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY)
      splitnumb_l<-1
      leftover<-0
    }else{
      splitnumb_l<-floor(nrow(GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY)/Splitrows)
      leftover<-splitnumb_raw-floor(splitnumb_l)
    }
    
    #Loops through cell chunks
    for(o in 1:splitnumb_l) {
      
      if(o==splitnumb_l){
        X1=(o-1)*Splitrows
        X2=Splitrows*o+Splitrows*leftover
      }else{o
        X1=(o-1)*Splitrows
        X2=Splitrows*o
      }  
      
      num.cores <- detectCores() ## detect how many cores on your machine
      cl <- makeCluster(num.cores, type = "SOCK")
      registerDoParallel(cl)
      getDoParWorkers() ## check if all cores are working
      clusterEvalQ(cl,library(RPostgreSQL))
      clusterEvalQ(cl,library(DBI))
      clusterEvalQ(cl,library(RSQLite))
      clusterEvalQ(cl,library(sp))
      clusterEvalQ(cl,library(rgeos))
      
      GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY_cut<-GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY[X1:X2,]
      
      #Read in bike station locations
      timeTR<-system.time({
        OO<-nrow(GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY_cut)
        
        REALOD_FINAL=foreach(i=1:OO,.combine = rbind,.errorhandling='remove') %dopar% {
          
          drv <- dbDriver("PostgreSQL")
          con <- dbConnect(drv, user = "postgres", dbname = "osm", host = "localhost")
          con_rt <- dbConnect(drv, user = "postgres", dbname = "opbeumDB", host = "localhost")
          
          O_loc<- GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY_cut[i,]
          O_xy<- coordinates(origins[as.numeric(O_loc[1]),])
          
          bbl<- O_xy*.9995  
          bbh<- O_xy*1.0005  
          
          ##
          ##get ORIGIN OSM NODE
          ##
          
          ##Find closest OSM node or the Origin
          select<-"SELECT 
          osm_id,
          ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
          ST_X(ST_Transform(geom_vertex,4326)) AS lon 
          FROM at_2po_vertex
          WHERE  
          geom_vertex && 
          ST_MakeEnvelope("
          
          bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
          q_tsig <- paste(select,bb2,")")
          
          #Pull SQL result
          tsig_result <- dbGetQuery(con_rt, q_tsig)
          
          #Convert node XY to spatial data frame
          node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          #convert cell XY to spatial data frame 
          origin.spdf <- SpatialPointsDataFrame(coords = O_xy, data = O_loc, 
                                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          
          #transform to planor coordinates
          node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) )
          origin.spdf.proj <- spTransform( origin.spdf, CRS( "+init=epsg:3347" ) ) 
          
          
          #get closest osm node
          near_node<-gDistance(origin.spdf.proj, spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
          node_dist<-min(near_node)
          node_id<-node.spdf$osm_id[which(near_node==min(near_node))]
          
          #get the 'source' node for routing
          select<-"SELECT source
          FROM at_2po_4pgr
          WHERE osm_source_id =" 
          q_node_origin <- paste(select,node_id)
          
          qnode_result_origin <- dbGetQuery(con_rt, q_node_origin)[1,]
          
          if(length(qnode_result_origin)<1){
            node_dist<-sort(near_node)
            node_id<-node.spdf$osm_id[which(near_node==near_node[2])]
            
            #get the 'source' node for routing
            select<-"SELECT source
            FROM at_2po_4pgr
            WHERE osm_source_id =" 
            q_node_origin <- paste(select,node_id)
            
            qnode_result_origin <- dbGetQuery(con_rt, q_node_origin)[1,]
          }
          
          
          #----------DESTINATION
          
          D_loc<- GRID_BIKESHARE_REALODs_GRIDInfo_OD_XY_cut[i,]
          D_xy<-  coordinates(origins[as.numeric(O_loc[2]),])
          
          bbl<- D_xy*.9995  
          bbh<- D_xy*1.0005  
          
          ##
          ##get ORIGIN OSM NODE
          ##
          
          ##Find closest OSM node or the Origin
          select<-"SELECT 
          osm_id,
          ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
          ST_X(ST_Transform(geom_vertex,4326)) AS lon 
          FROM at_2po_vertex
          WHERE  
          geom_vertex && 
          ST_MakeEnvelope("
          
          bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
          q_tsig <- paste(select,bb2,")")
          
          #Pull SQL result
          tsig_result <- dbGetQuery(con_rt, q_tsig)
          
          #Convert node XY to spatial data frame
          node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          #convert cell XY to spatial data frame 
          dest.spdf <- SpatialPointsDataFrame(coords = D_xy, data = D_loc, 
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
          
          #transform to planor coordinates
          node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) )
          dest.spdf.proj <- spTransform( dest.spdf, CRS( "+init=epsg:3347" ) ) 
          
          
          #get closest osm node
          near_node<-gDistance(dest.spdf.proj, spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
          node_dist<-min(near_node)
          node_id<-node.spdf$osm_id[which(near_node==min(near_node))]
          
          #get the 'source' node for routing
          select<-"SELECT source
          FROM at_2po_4pgr
          WHERE osm_source_id =" 
          #q_node_dest <- paste(select,node_id)
          q_node_dest <- paste(select,node_id)
          
          qnode_result_dest <- dbGetQuery(con_rt, q_node_dest)[1,]
          
          if(length(qnode_result_dest)<1){
            node_dist<-sort(near_node)
            node_id<-node.spdf$osm_id[which(near_node==near_node[2])]
            
            #get the 'source' node for routing
            select<-"SELECT source
            FROM at_2po_4pgr
            WHERE osm_source_id =" 
            q_node_origin <- paste(select,node_id)
            
            qnode_result_dest <- dbGetQuery(con_rt, q_node_dest)[1,]
          }
          
          
          if(length(qnode_result_origin)<1){
            ODDIS_n=9999
          }
          if(length(qnode_result_dest)<1){
            ODDIS_n=9999
          }else{
            
            #----------------------
            
            ##get OD dist
            
            q1<-   "SELECT seq, id1 AS node, id2 AS edge,km AS cost
            FROM pgr_astar('
            SELECT id AS id,
            source,
            target,
            cost,
            x1, y1, x2, y2
            FROM at_2po_4pgr as r,
            (SELECT ST_Expand(ST_Extent(geom_way),0.01) as box  FROM at_2po_4pgr as l1  
            WHERE l1.source ="
            o_node<-qnode_result_origin
            q2<-"OR l1.target =" 
            d_node<-qnode_result_dest
            q3<-") as box
            WHERE r.geom_way && box.box',"
            q4<-","
            q5<-", false, false)as r INNER JOIN at_2po_4pgr as g ON r.id2 = g.id ;"
            
            q_OD <- paste0(q1,o_node,q2,d_node,q3,o_node,q4,d_node,q5)
            
            d_result_OD <- dbGetQuery(con_rt, q_OD)
            
            ODDIS_n=sum(d_result_OD$cost)*0.621371
            
            #----------------------------  
            
            #end distance calcs  
          }
          dbDisconnect(con)
          dbDisconnect(con_rt)
          
          
          #Merge all location results
          cbind.data.frame(O_loc["PageNumber_O"],O_loc["PageNumber_D"],ODDIS_n)
          
        }  
        
        write.table(REALOD_FINAL, paste0("temp/",studyarea,"_BSOD.csv"), row.names = F, col.names = F, append = T, sep=",",quote=F)
        
        stopCluster(cl)
      })
      print(timeTR)
      print(splitnumb_l-o)
    }
    
    REALOD_FINAL2 = fread(paste0("temp/",studyarea,"_BSOD.csv"), sep=",",header = F)
    colnames(REALOD_FINAL2)<-c("PageNumber_O","PageNumber_D","ODDIS_n")
    #write.table(TRANSIT_FINAL2, paste0(studyarea,"2",".csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)
    
    REALOD_FINAL2_clean<-data.frame(REALOD_FINAL2[ which(REALOD_FINAL2$ODDIS_n<999), ])

  ##get elevation data
    for (e in 1:nrow(REALOD_FINAL2_clean)) {
      
      rodr<-REALOD_FINAL2_clean[e,]
      
      #get eleveation b/t OD pair
      src_el<-(coordinates2statistics(coordinates(origins[as.numeric(rodr[1]),])[2], coordinates(origins[as.numeric(rodr[1]),])[1], "elevation"))
      dst_el<-(coordinates2statistics(coordinates(origins[as.numeric(rodr[2]),])[2], coordinates(origins[as.numeric(rodr[2]),])[1], "elevation"))
      
      rodr$ELEVDELTA<-(dst_el$statistics.elevation.value-src_el$statistics.elevation.value)*3.28084
      
      rodr$SLOPE<- rodr$ELEVDELTA/ rodr$ODDIS_n
      #fix na slope
      rodr$SLOPE[is.na(rodr$SLOPE)] <- 0
      
      write.table(rodr, file=paste0("temp/","REAL_OD_DIST.csv"), row.names=FALSE, col.names=FALSE, append = T, sep=",", quote = F)
    }
    

##Read back complete Real OD distance file
REAL_OD_DIST = fread(paste0("temp/","REAL_OD_DIST.csv"), sep=",",header = TRUE)
names(REAL_OD_DIST) <- c("PageNumber_D", "PageNumber_O", "ODDIS_n", "ELEVDELTA", "SLOPE")
REAL_OD_DIST$SLOPE[is.infinite(REAL_OD_DIST$SLOPE)] <- 0

##Merge Distance into full Real OD Data
GRID_BIKESHARE_REALODs_GRIDInfo_OD<-merge(GRID_BIKESHARE_REALODs_GRIDInfo_OD,REAL_OD_DIST,by=c("PageNumber_D","PageNumber_O"))

#write.table(GRID_BIKESHARE_REALODs_GRIDInfo_OD, "GRID_BIKESHARE_REALODs_GRIDInfo_OD.csv", row.names = F, col.names = T, append = F, sep=",",quote=F)

##Fit basic SVM regression model
model <- ksvm(TRIPS ~ ZSCORE1_O + ZSCORE2_O + ZSCORE3_O + ZSCORE4_O + ZSCORE5_O + 
                      ZSCORE1_D + ZSCORE2_D + ZSCORE3_D + ZSCORE4_D + ZSCORE5_D +
                      ODDIS_n+ SLOPE,
               GRID_BIKESHARE_REALODs_GRIDInfo_OD,
              cost=10^10)


##Convert spatial dataframe to a dataframe
sInDF<-as.data.frame(sIn)
sInDF_SUB <- sInDF[c("PageNumber", "x", "y")]

#export shapes of potental stations
writeOGR(obj=sIn, dsn="output", layer=paste0("studyarea_","Acceptable_BS_Locations"), driver="ESRI Shapefile") 

##-------------FIND Modeled Origin OSM Node----------------##
if (file.exists(paste0("temp/",studyarea,"_MODELORIGIN.csv"))) file.remove(paste0("temp/",studyarea,"_MODELORIGIN.csv"))

MODELORIGIN_FINAL<- data.frame(
  PageNumber=numeric(),
  OSM_source=numeric(),
  stringsAsFactors=FALSE)

#Chunk grid cells 
Splitrows=1000
splitnumb_raw<-(nrow(sInDF_SUB)/Splitrows)
if(nrow(sInDF_SUB)<Splitrows){
  Splitrows<-nrow(sInDF_SUB)
  splitnumb_l<-1
  leftover<-0
}else{
  splitnumb_l<-floor(nrow(sInDF_SUB)/Splitrows)
  leftover<-splitnumb_raw-floor(splitnumb_l)
}

#Loops through cell chunks
for(o in 1:splitnumb_l) {
  
  if(o==splitnumb_l){
    X1=(o-1)*Splitrows
    X2=Splitrows*o+Splitrows*leftover
  }else{o
    X1=(o-1)*Splitrows
    X2=Splitrows*o
  }  
  
  num.cores <- detectCores() ## detect how many cores on your machine
  cl <- makeCluster(num.cores, type = "SOCK")
  registerDoParallel(cl)
  getDoParWorkers() ## check if all cores are working
  clusterEvalQ(cl,library(RPostgreSQL))
  clusterEvalQ(cl,library(DBI))
  clusterEvalQ(cl,library(RSQLite))
  clusterEvalQ(cl,library(sp))
  clusterEvalQ(cl,library(rgeos))
  
  sInDF_SUB_cut<-sInDF_SUB[X1:X2,]
  
  #Read  locations
  timeTR<-system.time({
    OO<-nrow(sInDF_SUB_cut)
    
    MODELORIGIN_FINAL=foreach(i=1:OO,.combine = rbind,.errorhandling='remove') %dopar% {
      drv <- dbDriver("PostgreSQL")
      con <- dbConnect(drv, user = "postgres", dbname = "osm", host = "localhost")
      con_rt <- dbConnect(drv, user = "postgres", dbname = "opbeumDB", host = "localhost")
      
      O_loc<- sInDF_SUB_cut[i,]
      O_xy<- coordinates(origins[as.numeric(O_loc[1]),])
      
      bbl<- O_xy*.9995  
      bbh<- O_xy*1.0005  
      
      ##
      ##get ORIGIN OSM NODE
      ##
      
      ##Find closest OSM node or the Origin
      select<-"SELECT 
      osm_id,
      ST_Y(ST_Transform(geom_vertex,4326)) AS lat, 
      ST_X(ST_Transform(geom_vertex,4326)) AS lon 
      FROM at_2po_vertex
      WHERE  
      geom_vertex && 
      ST_MakeEnvelope("
      
      bb2<-paste(bbl[1],",",bbl[2],",",bbh[1],",",bbh[2],",","4283")
      q_tsig <- paste(select,bb2,")")
      
      #Pull SQL result
      tsig_result <- dbGetQuery(con_rt, q_tsig)
      
      #Convert node XY to spatial data frame
      node.spdf <- SpatialPointsDataFrame(coords = tsig_result[,c(3,2)], data = tsig_result,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      #convert cell XY to spatial data frame 
      origin.spdf <- SpatialPointsDataFrame(coords = O_xy, data = O_loc, 
                                            proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      
      #transform to planor coordinates
      node.spdf.proj <- spTransform( node.spdf, CRS( "+init=epsg:3347" ) )
      origin.spdf.proj <- spTransform( origin.spdf, CRS( "+init=epsg:3347" ) ) 
      
      
      #get closest osm node
      near_node<-gDistance(origin.spdf.proj, spgeom2=node.spdf.proj, byid=T, hausdorff=FALSE, densifyFrac = NULL)
      node_dist<-min(near_node)
      node_id<-node.spdf$osm_id[which(near_node==min(near_node))]
      
      #get the 'source' node for routing
      select<-"SELECT source
      FROM at_2po_4pgr
      WHERE osm_source_id =" 
      q_node_origin <- paste(select,node_id)
      
      qnode_result_origin <- dbGetQuery(con_rt, q_node_origin)[1,]
      
      if(length(qnode_result_origin)<1){
        node_dist<-sort(near_node)
        node_id<-node.spdf$osm_id[which(near_node==near_node[2])]
        
        #get the 'source' node for routing
        select<-"SELECT source
        FROM at_2po_4pgr
        WHERE osm_source_id =" 
        q_node_origin <- paste(select,node_id)
        
        qnode_result_origin <- dbGetQuery(con_rt, q_node_origin)[1,]
      }
      
      ##3rd try
      if(length(qnode_result_origin)<1){
        node_dist<-sort(near_node)
        node_id<-node.spdf$osm_id[which(near_node==near_node[3])]
        
        #get the 'source' node for routing
        select<-"SELECT source
        FROM at_2po_4pgr
        WHERE osm_source_id =" 
        q_node_origin <- paste(select,node_id)
        
        qnode_result_origin <- dbGetQuery(con_rt, q_node_origin)[1,]
      }
      
      if(length(qnode_result_origin)<1){qnode_result_origin=0}
        
      dbDisconnect(con)
      dbDisconnect(con_rt)
      
      OSM_source<-qnode_result_origin
      
      #Merge all location results
      cbind.data.frame(O_loc["PageNumber"],OSM_source)
    }  
    
    write.table(MODELORIGIN_FINAL, paste0("temp/",studyarea,"_MODELORIGIN.csv"), row.names = F, col.names = F, append = T, sep=",",quote=F)
    
    stopCluster(cl)
  })
  print(timeTR)
  print(splitnumb_l-o)
}

MODELORIGIN_FINAL2 = fread(paste0("temp/",studyarea,"_MODELORIGIN.csv"), sep=",",header = F)
colnames(MODELORIGIN_FINAL2)<-c("PageNumber","OSM_source")

sInDF_SUB<-merge(sInDF_SUB,MODELORIGIN_FINAL2,"PageNumber")

#Get Elevation data for potential cells
sInDF_SUB$elevation<-NA
for (v in 1:nrow(sInDF_SUB)) {
  el<-as.data.frame(coordinates2statistics(sInDF_SUB$y[v], sInDF_SUB$x[v], "elevation"))
  el2<-el$statistics.elevation.value*3.28084
  sInDF_SUB$elevation[v]<-el2 
}

##Get unique combinations,convert to a dataframe, merge with XY data
sInDF_SUB_Comb<-data.frame((do.call(rbind,combn(sInDF_SUB$PageNumber,2,simplify = FALSE))))
  names(sInDF_SUB_Comb) <- c("PageNumber_O","PageNumber_D")
sInDF_SUB_Comb_OD <- merge(sInDF_SUB_Comb, sInDF_SUB, by.x = "PageNumber_O", by.y = "PageNumber")
  names(sInDF_SUB_Comb_OD) <- c("PageNumber_O","PageNumber_D", "x_O", "y_O", "OSM_source_O","elevation_O")
sInDF_SUB_Comb_OD <- merge(sInDF_SUB_Comb_OD, sInDF_SUB, by.x = "PageNumber_D", by.y = "PageNumber")
  names(sInDF_SUB_Comb_OD) <- c("PageNumber_O","PageNumber_D", "x_O", "y_O", "OSM_source_O","elevation_O","x_D", "y_D", "OSM_source_D","elevation_D")

  
##Simpple earth distance formula
  earth.dist <- function (long1, lat1, long2, lat2)
  {
    rad <- pi/180
    a1 <- lat1 * rad
    a2 <- long1 * rad
    b1 <- lat2 * rad
    b2 <- long2 * rad
    dlon <- b2 - a2
    dlat <- b1 - a1
    a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
    c <- 2 * atan2(sqrt(a), sqrt(1 - a))
    R <- 6378.145
    d <- R * c
    return(d)
  }

  ##first get straight line distance
  sInDF_SUB_Comb_OD$earth_dist <- sapply(1:nrow(sInDF_SUB_Comb_OD),function(i)
    earth.dist(sInDF_SUB_Comb_OD$y_O[i],sInDF_SUB_Comb_OD$x_O[i],sInDF_SUB_Comb_OD$y_D[i],sInDF_SUB_Comb_OD$x_D[i]))
  ##reduce dataset to 10 miles per pair
  sInDF_SUB_Comb_OD$earth_dist<-sInDF_SUB_Comb_OD$earth_dist*0.621371
  sInDF_SUB_Comb_OD_red<-sInDF_SUB_Comb_OD[ which(sInDF_SUB_Comb_OD$earth_dist<=10), ]



##-------------Calc netwotk distance between cells----------------##
if (file.exists(paste0("temp/",studyarea,"_MODELOD.csv"))) file.remove(paste0("temp/",studyarea,"_MODELOD.csv"))


MODELOD_FINAL<- data.frame(
  PageNumber_O=numeric(),
  PageNumber_D=numeric(),
  n_dist=numeric(),
  stringsAsFactors=FALSE)

#Chunk grid cells 
Splitrows=1000
splitnumb_raw<-(nrow(sInDF_SUB_Comb_OD_red)/Splitrows)
if(nrow(sInDF_SUB_Comb_OD_red)<Splitrows){
  Splitrows<-nrow(sInDF_SUB_Comb_OD_red)
  splitnumb_l<-1
  leftover<-0
}else{
  splitnumb_l<-floor(nrow(sInDF_SUB_Comb_OD_red)/Splitrows)
  leftover<-splitnumb_raw-floor(splitnumb_l)
}

#Loops through cell chunks
for(o in 1:splitnumb_l) {
  
  if(o==splitnumb_l){
    X1=(o-1)*Splitrows
    X2=Splitrows*o+Splitrows*leftover
  }else{o
    X1=(o-1)*Splitrows
    X2=Splitrows*o
  }  
  
  num.cores <- detectCores() ## detect how many cores on your machine
  cl <- makeCluster(num.cores, type = "SOCK")
  registerDoParallel(cl)
  getDoParWorkers() ## check if all cores are working
  clusterEvalQ(cl,library(RPostgreSQL))
  clusterEvalQ(cl,library(DBI))
  clusterEvalQ(cl,library(RSQLite))
  clusterEvalQ(cl,library(sp))
  clusterEvalQ(cl,library(rgeos))
  
  sInDF_SUB_Comb_OD_red_cut<-sInDF_SUB_Comb_OD_red[X1:X2,]
  
  #Read  locations
  timeTR<-system.time({
    OO<-nrow(sInDF_SUB_Comb_OD_red_cut)
    
    MODELOD_FINAL=foreach(i=1:OO,.combine = rbind,.errorhandling='remove') %dopar% {
      
      drv <- dbDriver("PostgreSQL")
      con <- dbConnect(drv, user = "postgres", dbname = "osm", host = "localhost")
      con_rt <- dbConnect(drv, user = "postgres", dbname = "opbeumDB", host = "localhost")
      
      loc<-sInDF_SUB_Comb_OD_red_cut[i,]
      
      OSM_source_O<- loc$OSM_source_O
      OSM_source_D<- loc$OSM_source_D
      
      
      #----------------------
      
      ##get OD dist
      
      q1<-   "SELECT seq, id1 AS node, id2 AS edge,km AS cost
      FROM pgr_astar('
      SELECT id AS id,
      source,
      target,
      cost,
      x1, y1, x2, y2
      FROM at_2po_4pgr as r,
      (SELECT ST_Expand(ST_Extent(geom_way),0.01) as box  FROM at_2po_4pgr as l1  
      WHERE l1.source ="
      o_node<-OSM_source_O
      q2<-"OR l1.target =" 
      d_node<-OSM_source_D
      q3<-") as box
      WHERE r.geom_way && box.box',"
      q4<-","
      q5<-", false, false)as r INNER JOIN at_2po_4pgr as g ON r.id2 = g.id ;"
      
      q_OD <- paste0(q1,o_node,q2,d_node,q3,o_node,q4,d_node,q5)
      
      d_result_OD <- dbGetQuery(con_rt, q_OD)
      
      n_dist=sum(d_result_OD$cost)*0.621371
      
      #----------------------------  
      
      #end distance calcs  
      
      dbDisconnect(con)
      dbDisconnect(con_rt)
      
      
      #Merge all location results
      cbind.data.frame(loc["PageNumber_O"],loc["PageNumber_D"],n_dist)
      
    }  
    
    write.table(MODELOD_FINAL, paste0("temp/",studyarea,"_MODELOD.csv"), row.names = F, col.names = F, append = T, sep=",",quote=F)
    
    stopCluster(cl)
  })
  print(timeTR)
  print(splitnumb_l-o)
}

MODELOD_FINAL2 = fread(paste0("temp/",studyarea,"_MODELOD.csv"), sep=",",header = F)
colnames(MODELOD_FINAL2)<-c("PageNumber_O","PageNumber_D","ODDIS_n")

    #merge back distances with elevation- for the record
    sInDF_SUB_Comb_OD_red<-merge(sInDF_SUB_Comb_OD_red,MODELOD_FINAL2, c("PageNumber_O","PageNumber_D"))

        ### Get rid of pairs over 10 (network) miles
        sInDF_SUB_Comb_OD_red_final<-sInDF_SUB_Comb_OD_red[ which(sInDF_SUB_Comb_OD_red$ODDIS_n<=10), ]

      ##calculate final elevation change and slope
      sInDF_SUB_Comb_OD_red_final$ELEVDELTA<-(sInDF_SUB_Comb_OD_red_final$elevation_D-sInDF_SUB_Comb_OD_red_final$elevation_O)
      sInDF_SUB_Comb_OD_red_final$SLOPE<- sInDF_SUB_Comb_OD_red_final$ELEVDELTA/ sInDF_SUB_Comb_OD_red_final$ODDIS_n

      ##fix na and inf slope
      sInDF_SUB_Comb_OD_red_final$SLOPE[is.na(sInDF_SUB_Comb_OD_red_final$SLOPE)] <- 0
      sInDF_SUB_Comb_OD_red_final$SLOPE[is.infinite(sInDF_SUB_Comb_OD_red_final$SLOPE)] <- 0 

            ##Write final potential OD file - with elevation data
            write.table(sInDF_SUB_Comb_OD_red_final, paste0("temp/",studyarea,"potential_OD_Final.csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)

##extract final OD pagenumber only
 sInDF_SUB_Comb_OD_red_final_ODPageNumber<-sInDF_SUB_Comb_OD_red_final[c("PageNumber_O","PageNumber_D")]     
            
##Join O and D Cell info
BASEcells_Reduced_Acceptable_O<-merge(sInDF_SUB_Comb_OD_red_final_ODPageNumber,BASEcells,by.x=c("PageNumber_O"), by.y=c("PageNumber"))
  names(BASEcells_Reduced_Acceptable_O)[3:length(names(BASEcells_Reduced_Acceptable_O))]<-paste(names(BASEcells_Reduced_Acceptable_O)[3:length(names(BASEcells_Reduced_Acceptable_O))],"O", sep = "_")
BASEcells_Reduced_Acceptable_OD<-merge(BASEcells_Reduced_Acceptable_O,BASEcells,by.x=c("PageNumber_D"), by.y=c("PageNumber"))
  names(BASEcells_Reduced_Acceptable_OD)[(length(names(BASEcells_Reduced_Acceptable_O))+1):length(names(BASEcells_Reduced_Acceptable_OD))]<-paste(names(BASEcells_Reduced_Acceptable_OD)[(length(names(BASEcells_Reduced_Acceptable_O))+1):length(names(BASEcells_Reduced_Acceptable_OD))],"D", sep = "_")

  ##merge in distance etc.   
sInDF_OD_DIST_n_Final<-merge(BASEcells_Reduced_Acceptable_OD,sInDF_SUB_Comb_OD_red_final,by=c("PageNumber_O","PageNumber_D"))

##Predict trips and join to rest of data
predictedY <- predict(model, sInDF_OD_DIST_n_Final)
sInDF_OD_DIST_n_Final_RESULT<-cbind(sInDF_OD_DIST_n_Final,predictedY)

###Calibrated reported trips
cal_base<-min(sInDF_OD_DIST_n_Final_RESULT$predictedY)
if(cal_base<0){
  cal_rate<-abs(cal_base)
  sInDF_OD_DIST_n_Final_RESULT$predictedY<-sInDF_OD_DIST_n_Final_RESULT$predictedY+cal_rate
}
##Write final Basic SVM result
write.table(sInDF_OD_DIST_n_Final_RESULT, paste0("Output/",studyarea,"SVM_Regression_Result_Final.csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)

###Finsih steve's analysis

##adjust##
sInDF_OD_DIST_n_Final_RESULT$Annual_Trips<-sInDF_OD_DIST_n_Final_RESULT$predictedY*.249

radjust<-6286/(sum(GRID_BIKESHARE_REALODs_GRIDInfo_OD$TRIPS)/360)
GRID_BIKESHARE_REALODs_GRIDInfo_OD$TRIPS_adj<-GRID_BIKESHARE_REALODs_GRIDInfo_OD$TRIPS*radjust

MODEL_Trips_byZSCORE_sub<-sInDF_OD_DIST_n_Final_RESULT[c("CLUSTER_O","CLUSTER_D", "Annual_Trips")]                                                                  
MODEL_Trips_byZSCORE <- aggregate(MODEL_Trips_byZSCORE_sub, list(MODEL_Trips_byZSCORE_sub$CLUSTER_O,MODEL_Trips_byZSCORE_sub$CLUSTER_D), sum)
  MODEL_Trips_byZSCORE<-MODEL_Trips_byZSCORE[c("Group.1","Group.2", "Annual_Trips")]     
    names(MODEL_Trips_byZSCORE)<-c("CLUSTER_O","CLUSTER_D", "Annual_Trips")

    write.table(MODEL_Trips_byZSCORE, paste0("Output/",studyarea,"Cluster_Modeled_OD_Result.csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)
    

   REAL_Trips_byZSCORE_sub<-GRID_BIKESHARE_REALODs_GRIDInfo_OD[c("CLUSTER_O","CLUSTER_D", "TRIPS_adj")]                                                                  
   REAL_Trips_byZSCORE <- aggregate(REAL_Trips_byZSCORE_sub, list(REAL_Trips_byZSCORE_sub$CLUSTER_O,REAL_Trips_byZSCORE_sub$CLUSTER_D), sum)
    REAL_Trips_byZSCORE<-REAL_Trips_byZSCORE[c("Group.1","Group.2", "TRIPS_adj")]     
    names(REAL_Trips_byZSCORE)<-c("CLUSTER_O","CLUSTER_D", "Annual_Trips")
    
    write.table(REAL_Trips_byZSCORE, paste0("Output/",studyarea,"Cluster_Real_OD_Result.csv"), row.names = F, col.names = T, append = F, sep=",",quote=F)
    
    x<-sort(GRID_BIKESHARE_REALODs$Start.date)
    
#############
#OD MAPS
############

##Plot
CellXY_MAP<-fread(paste0("Output/",studyarea,"SVM_Regression_Result_Final.csv"), sep=",",header = TRUE, stringsAsFactors=FALSE)

##Add XY data to Real OD trips
# CellXY_MAP$PageNumber_O<-CellXY_MAP$PAGENUMBER
# CellXY_MAP$TERMINAL_N_O<-CellXY_MAP$TERMINAL_N
# CellXY_MAP$oX<-CellXY_MAP$LONGITUDE
# CellXY_MAP$oY<-CellXY_MAP$LATITUDE
CellXY_MAP_O<-subset(CellXY_MAP,select = c("PageNumber_O","TERMINAL_N_O","o_x_D", "o_y_D"))

# CellXY_MAP$PageNumber_D<-CellXY_MAP$PAGENUMBER
# CellXY_MAP$TERMINAL_N_D<-CellXY_MAP$TERMINAL_N
# CellXY_MAP$dX<-CellXY_MAP$LONGITUDE
# CellXY_MAP$dY<-CellXY_MAP$LATITUDE
CellXY_MAP_D<-subset(CellXY_MAP,select = c("PageNumber_D","TERMINAL_N_D","o_x_D", "o_y_D"))

GRID_BIKESHARE_REALODs_AGG_TRIPS<-merge(GRID_BIKESHARE_REALODs_AGG,CellXY_MAP_O,by="TERMINAL_N_O")
GRID_BIKESHARE_REALODs_AGG_TRIPS<-merge(GRID_BIKESHARE_REALODs_AGG_TRIPS,CellXY_MAP_D,by="TERMINAL_N_D")

MAX_Y<-max(GRID_BIKESHARE_REALODs_AGG_TRIPS$oY)
MIN_Y<-min(GRID_BIKESHARE_REALODs_AGG_TRIPS$oY)
MAX_X<-max(GRID_BIKESHARE_REALODs_AGG_TRIPS$oX)
MIN_X<-min(GRID_BIKESHARE_REALODs_AGG_TRIPS$oX)


# Load station points
Stations<-fread("GRID_BIKESHARE_Stationsx.csv", sep=",",header = TRUE, stringsAsFactors=FALSE)
Stations$x <- as.numeric(as.character(Stations$LONGITUDE))
Stations$y <- as.numeric(as.character(Stations$LATITUDE))
coordinates(Stations) <- ~x + y
Stations_DF<-as.data.frame(Stations)

##load background map
b <- bbox(Stations)
MOCOMap = ggmap(get_map(location = b,scale = "auto", color = "bw"))


#####Map Real ODs
R_LAB_MAX <- as.data.frame(Stations[ which(Stations$TRIPS_O>3000), ])
R_LAB_MIN <- as.data.frame(Stations[ which(Stations$TRIPS_O>0 & Stations$TRIPS_O<100), ])

xquiet<- scale_x_continuous("", breaks=NULL)
yquiet<-scale_y_continuous("", breaks=NULL)
quiet<-list(xquiet, yquiet)

# shape<-readShapeSpatial("cb_2015_24_place_500k.shp") 
# shape@data$id <- rownames(shape@data)
# shape.fort <- fortify(shape, region='id') 
MOCOMap+
  #ggplot() + 
  scale_alpha_continuous(range = c(.3, .3)) +
  theme(panel.background = element_rect(fill='white',colour='white')) +
  geom_point(data=Stations_DF, mapping=aes(x=LONGITUDE, y=LATITUDE, size=TRIPS_O),color="skyblue4",alpha = 0.45)+
  geom_segment(aes(x=oX, y=oY,xend=dX, yend=dY, alpha=TRIPS), col = "skyblue4",data = GRID_BIKESHARE_REALODs_AGG_TRIPS)+
  
  geom_label_repel(data=R_LAB_MAX,aes(LONGITUDE, LATITUDE, label = ADDRESS, fill ='tomato4 ', alpha = 0.75), color = 'white',size=1.8,fontface="bold",label.padding = unit(0.05, "lines")) +
  geom_label_repel(data=R_LAB_MIN,aes(LONGITUDE, LATITUDE, label = ADDRESS, fill ='darkgreen', alpha = 0.75), color = 'white',size=1.8,fontface="bold",label.padding = unit(0.05, "lines")) +
  
  ggtitle("Existing Maryland Capital Bikeshare Stations") +
  quiet+coord_equal()

dev.copy2pdf(file="9_RealODs.pdf")

#####Map Modeled ODs

BIKESHARE_ODModel_sInODs_Final_RESULT_XY<-as.data.frame(sInDF_OD_DIST_n_Final_RESULT[ which(sInDF_OD_DIST_n_Final_RESULT$predictedY>100), ])
BIKESHARE_ODModel_sInODs_Final_RESULT_XY<-unique(as.data.frame(BIKESHARE_ODModel_sInODs_Final_RESULT_XY[ which(BIKESHARE_ODModel_sInODs_Final_RESULT_XY$DIST<52800), ]))

CellXY_ALL_O<-CellXY_ALL
CellXY_ALL_D<-CellXY_ALL

CellXY_ALL_O<-rename(CellXY_ALL_O, c("PAGENUMBER"="PageNumber_O", "X"="oX","Y"="oY"))
CellXY_ALL_D<-rename(CellXY_ALL_D, c("PAGENUMBER"="PageNumber_D", "X"="dX","Y"="dY"))

BIKESHARE_ODModel_sInODs_Final_RESULT_XY<-merge(BIKESHARE_ODModel_sInODs_Final_RESULT_XY,CellXY_ALL_O,by="PageNumber_O")
BIKESHARE_ODModel_sInODs_Final_RESULT_XY<-merge(BIKESHARE_ODModel_sInODs_Final_RESULT_XY,CellXY_ALL_D,by="PageNumber_D")

Final_Os<-unique(subset(BIKESHARE_ODModel_sInODs_Final_RESULT_XY, select = c("PageNumber_O")))
Final_Os$USE<-1
CellXY_ALL_Points<-CellXY_ALL
CellXY_ALL_Points$PageNumber_O<-CellXY_ALL_Points$PAGENUMBER
CellXY_ALL_Points<-merge(CellXY_ALL_Points,Final_Os,by="PageNumber_O")
CellXY_ALL_Points<-CellXY_ALL_Points[ which(CellXY_ALL_Points$USE==1), ]


# Spatialdataframe the cells
CellXY_ALL_Points$x <- as.numeric(as.character(CellXY_ALL_Points$X))
CellXY_ALL_Points$y <- as.numeric(as.character(CellXY_ALL_Points$Y))
coordinates(CellXY_ALL_Points) <- ~x + y

BIKESHARE_ODModel_sInODs_Final_RESULT_XY_SDF<-BIKESHARE_ODModel_sInODs_Final_RESULT_XY
BIKESHARE_ODModel_sInODs_Final_RESULT_XY_SDF$x <- as.numeric(as.character(BIKESHARE_ODModel_sInODs_Final_RESULT_XY_SDF$oX))
BIKESHARE_ODModel_sInODs_Final_RESULT_XY_SDF$y <- as.numeric(as.character(BIKESHARE_ODModel_sInODs_Final_RESULT_XY_SDF$oY))
coordinates(BIKESHARE_ODModel_sInODs_Final_RESULT_XY_SDF) <- ~x + y

c <- bbox(CellXY_ALL_Points)
d = matrix(c(1.05, 0 , 1.05, 0), nrow=2, ncol=2) 
e<-c+d
MDMAP = ggmap(get_map(location = e,scale = "auto", color = "bw"))

xquiet<- scale_x_continuous("", breaks=NULL)
yquiet<-scale_y_continuous("", breaks=NULL)
quiet<-list(xquiet, yquiet)
MDMAP+
  #ggplot() + 
  scale_alpha_continuous(range = c(.3, .3)) +
  theme(panel.background = element_rect(fill='white',colour='white')) +
  geom_segment(aes(x=oX, y=oY,xend=dX, yend=dY, alpha=predictedY), col = "pink4",data = BIKESHARE_ODModel_sInODs_Final_RESULT_XY)+
  #geom_point(data=Stations_DF, mapping=aes(x=LONGITUDE, y=LATITUDE, size=TRIPS_O),color="white")+
  #geom_text_repel(data=R_LAB_MAX,aes(LONGITUDE, LATITUDE, label = ADDRESS), color = 'chartreuse1',size=3,fontface="bold") +
  #geom_text_repel(data=R_LAB_MIN,aes(LONGITUDE, LATITUDE, label = ADDRESS), color = 'firebrick1',size=3,fontface="bold") +
  ggtitle("Potential Maryland Bikeshare Stations") +
  quiet+coord_equal()
dev.copy2pdf(file="10_ModelODs.pdf")

##Write final Basic SVM result

write.table(BIKESHARE_ODModel_sInODs_Final_RESULT_XY, "SVM_Regression_Result_Final_XY.csv", row.names = F, col.names = T, append = F, sep=",",quote=F)

xquiet<- scale_x_continuous("", breaks=NULL)
yquiet<-scale_y_continuous("", breaks=NULL)
quiet<-list(xquiet, yquiet)
#MDMAP+
ggplot() + 
  scale_alpha_continuous(range = c(.3, .3)) +
  theme(panel.background = element_rect(fill='black',colour='blue')) +
  geom_segment(aes(x=oX, y=oY,xend=dX, yend=dY, alpha=predictedY), col = "white",data = BIKESHARE_ODModel_sInODs_Final_RESULT_XY)+
  #geom_point(data=Stations_DF, mapping=aes(x=LONGITUDE, y=LATITUDE, size=TRIPS_O),color="white")+
  #geom_text_repel(data=R_LAB_MAX,aes(LONGITUDE, LATITUDE, label = ADDRESS), color = 'chartreuse1',size=3,fontface="bold") +
  #geom_text_repel(data=R_LAB_MIN,aes(LONGITUDE, LATITUDE, label = ADDRESS), color = 'firebrick1',size=3,fontface="bold") +
  #ggtitle("Potential Maryland Bikeshare Stations") +
  quiet+coord_equal()
dev.copy2pdf(file="graphic.pdf")
