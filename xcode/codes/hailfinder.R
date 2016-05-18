# install.packages("devtools")
# library(devtools)
# if(devtools::find_rtools()) devtools::install_github("itsrainingdata/ccdr")
# library(ccdr) ## manually build this package and load

setwd("~/ccdr/xcode/codes") ## change if necessary

## test on established graphs from `bnlearn` package
# install.packages("bnlearn")
library(bnlearn)

# try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("graph")
library(graph)

# install.packages("pcalg")
# try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("RBGL")
# library(RBGL)
library(pcalg)

# library(igraph)

source("rmvDAG_fix.R") ## generate random data with points allowed to be fixed (in a specified range)
source("convert.R")
source("compare.R") ## compare graphs. different from `compareGraphs`
source("maintest.R") ## main test
source("swaptest.R") ## outdated method if in each test we just swap columns on the same sample

data(hailfinder)
res = empty.graph(names(hailfinder))
modelstring(res) = paste("[N07muVerMo][SubjVertMo][QGVertMotion][SatContMoist][RaoContMoist]",
                         "[VISCloudCov][IRCloudCover][AMInstabMt][WndHodograph][MorningBound][LoLevMoistAd][Date]",
                         "[MorningCIN][LIfr12ZDENSd][AMDewptCalPl][LatestCIN][LLIW]",
                         "[CombVerMo|N07muVerMo:SubjVertMo:QGVertMotion][CombMoisture|SatContMoist:RaoContMoist]",
                         "[CombClouds|VISCloudCov:IRCloudCover][Scenario|Date][CurPropConv|LatestCIN:LLIW]",
                         "[AreaMesoALS|CombVerMo][ScenRelAMCIN|Scenario][ScenRelAMIns|Scenario][ScenRel34|Scenario]",
                         "[ScnRelPlFcst|Scenario][Dewpoints|Scenario][LowLLapse|Scenario][MeanRH|Scenario]",
                         "[MidLLapse|Scenario][MvmtFeatures|Scenario][RHRatio|Scenario][SfcWndShfDis|Scenario]",
                         "[SynForcng|Scenario][TempDis|Scenario][WindAloft|Scenario][WindFieldMt|Scenario]",
                         "[WindFieldPln|Scenario][AreaMoDryAir|AreaMesoALS:CombMoisture]",
                         "[AMCINInScen|ScenRelAMCIN:MorningCIN][AMInsWliScen|ScenRelAMIns:LIfr12ZDENSd:AMDewptCalPl]",
                         "[CldShadeOth|AreaMesoALS:AreaMoDryAir:CombClouds][InsInMt|CldShadeOth:AMInstabMt]",
                         "[OutflowFrMt|InsInMt:WndHodograph][CldShadeConv|InsInMt:WndHodograph][MountainFcst|InsInMt]",
                         "[Boundaries|WndHodograph:OutflowFrMt:MorningBound][N34StarFcst|ScenRel34:PlainsFcst]",
                         "[CompPlFcst|AreaMesoALS:CldShadeOth:Boundaries:CldShadeConv][CapChange|CompPlFcst]",
                         "[InsChange|CompPlFcst:LoLevMoistAd][CapInScen|CapChange:AMCINInScen]",
                         "[InsSclInScen|InsChange:AMInsWliScen][R5Fcst|MountainFcst:N34StarFcst]",
                         "[PlainsFcst|CapInScen:InsSclInScen:CurPropConv:ScnRelPlFcst]",
                         sep = "")
# there are too many nodes for plot(), use graphviz.plot().
g.res <- as.graphNEL(res)

pp <- numNodes(g.res)
nodes(g.res) <- as.character(1:pp)
vfix <- rep(sample(1:pp), 10)
test <- maintest(g.res, vfix = vfix, N = 50)

## focus on reversed edges
## get edges with fewer true estimates and many more reversed estimates
redge0 <- redge(test)
redge0
## and true edges
tedge0 <- tedge(test)
tedge0

revnodes <- c(1, 2, 3, 4, 5) # change this to whatever nodes we want to add intervention
vfix.rev <- rep(revnodes[sample(length(revnodes))], 50)
# test.rev <- maintest(g.res, vfix.rev, N = 50)
test.rev <- swaptest(g.res, vfix.rev, N = 50, originaldata = test$data, originalvfix = test$vfix)

colMeans(test.rev$metric)
apply(test.rev$metric, 2, sd)
## see if any changes
summarynewtest(test.rev, redge0)
summarynewtest(test.rev, tedge0)
