

require(xts)
require(quantmod)
#Download data
symbols <- c("VLO","BA","HD","CL","GSK","WBA","GS","MCD","PEP","XOM","TM","PFE","CVX","SNE","SAP","F","GE","UN","NVS","IBM","MTU","KMB","MAR","MDLZ","CAJ","HMC","CMCSA","AXP","WMT","AIG","R","HPQ","TXN","GD","LMT","DD","RTN","NOC","CAT","K","SNY","WFC","MMM","BAC","JPM","NAV","COP","CSCO","TOT","KO","CVS","PG","XRX","TWX")
symbols <- getSymbols.yahoo(symbols, from = "2008-11-01", env = globalenv())
symData <- xts()
for(i in seq_along(symbols)) {
  symbol <- symbols[i]
  symData <- merge(symData, get(symbol)[,paste(symbol, "Close", sep=".")])
}

# Partial correlation network
part_cor_mat <- matrix(data=0L, ncol=NCOL(symData),nrow=NCOL(symData))
colnames(part_cor_mat)<-symbols
rownames(part_cor_mat)<-symbols
for(z in 1:NCOL(symData))
{
  for (x in 1:NCOL(symData))
  {
    for (y in 1:NCOL(symData))
    {
      if(x!=y && x!=z && y!=z)
      {
        corxy = cor(coredata(Cl(symData))[,x], coredata(Cl(symData))[,y])
        corxz = cor(coredata(Cl(symData))[,x], coredata(Cl(symData))[,z])
        coryz = cor(coredata(Cl(symData))[,y], coredata(Cl(symData))[,z])
        corxyz = (corxy - corxz*coryz)/sqrt((1-corxz^2)*(1-coryz^2))
        part_cor_mat[z,x] <- part_cor_mat[z,x] + (corxy - corxyz)
        part_cor_mat[z,y] <- part_cor_mat[z,y] + (corxy - corxyz)
      }
    }
  }
}
threshold_factor_part <- 1
part_threshold_value <- mean(part_cor_mat[part_cor_mat>0]) + sd(part_cor_mat[part_cor_mat>0])*threshold_factor_part
part_cor_mat[ abs(part_cor_mat) < part_threshold_value]<- 0

library(igraph)
png(filename = "Partial_correlation_network.png", width = 10000, height = 10000, res = 200)
test.net1 = graph_from_adjacency_matrix(part_cor_mat,
                                        mode = "directed", weighted = TRUE) 
E(test.net1)$weight<-t(part_cor_mat)[abs(t(part_cor_mat))>part_threshold_value]
E(test.net1)$color <- (abs(E(test.net1)$weight)-min(part_cor_mat[part_cor_mat>0]))/(max(part_cor_mat)-min(part_cor_mat[part_cor_mat>0])) * 50 + 1
V(test.net1)$color <- ifelse(degree(test.net1, v=V(test.net1), mode="out")>degree(test.net1, v=V(test.net1), mode="in"), yes=adjustcolor("green", alpha.f = .6), no=adjustcolor("red", alpha.f = .6))
V(test.net1)$size <- 3 + degree(test.net1, v=V(test.net1), mode="out")/8
V(test.net1)$label.cex <- degree(test.net1, v=V(test.net1), mode="out")/(degree(test.net1, v=V(test.net1), mode="out")-1)
par(bg = "white")

plot(test.net1, rescale=T, frame=T,edge.curved=T, layout=layout_with_fr(test.net1),
     vertex.label.color= "black", edge.label.color= "black", vertex.frame.color="white",
     edge.width=7,
     edge.arrow.size = 3,
     edge.label=round(E(test.net1)$weight, 2),
     palette=palette(heat.colors(55, alpha=1)))
dev.off()
# http://finance.yahoo.com/chart/JPM#eyJjb21wYXJpc29ucyI6IkJBQyxHUyIsImNvbXBhcmlzb25zQ29sb3JzIjoiIzFhYzU2NywjZjAxMjZmIiwiY29tcGFyaXNvbnNHaG9zdGluZyI6IjAsMCIsImNvbXBhcmlzb25zV2lkdGhzIjoiMSwxIiwibXVsdGlDb2xvckxpbmUiOmZhbHNlLCJib2xsaW5nZXJVcHBlckNvbG9yIjoiI2UyMDA4MSIsImJvbGxpbmdlckxvd2VyQ29sb3IiOiIjOTU1MmZmIiwibWZpTGluZUNvbG9yIjoiIzQ1ZTNmZiIsIm1hY2REaXZlcmdlbmNlQ29sb3IiOiIjZmY3YjEyIiwibWFjZE1hY2RDb2xvciI6IiM3ODdkODIiLCJtYWNkU2lnbmFsQ29sb3IiOiIjMDAwMDAwIiwicnNpTGluZUNvbG9yIjoiI2ZmYjcwMCIsInN0b2NoS0xpbmVDb2xvciI6IiNmZmI3MDAiLCJzdG9jaERMaW5lQ29sb3IiOiIjNDVlM2ZmIiwicmFuZ2UiOiIxeSJ9



#create correlation matrix
cor_mat<-matrix(ncol=NCOL(symData),nrow=NCOL(symData))
colnames(cor_mat)<-symbols
rownames(cor_mat)<-symbols

for (i in 1:NCOL(symData))
{
  for (j in 1:NCOL(symData))
  {
    if(i!=j)
    {
      cor_mat[i,j]<-cor(coredata(Cl(symData))[,i], coredata(Cl(symData))[,j])
    }
    else
    {
      cor_mat[i,j]<-1
    }
  }
}

threshold_value <- 0.8
cor_mat[ lower.tri(cor_mat, diag=TRUE) ]<- 0
cor_mat[ abs(cor_mat) < threshold_value]<- 0

library(igraph)
png(filename = "Correlation_network.png", width = 10000, height = 10000, res = 200)
test.net1 = graph_from_adjacency_matrix(cor_mat,
                                        mode = "undirected", weighted = TRUE)
E(test.net1)$weight<-t(cor_mat)[abs(t(cor_mat))>threshold_value]
E(test.net1)$color <- ifelse(test=E(test.net1)$weight*100>0, yes=E(test.net1)$weight*100, no=E(test.net1)$weight*100+100)#100 - abs(E(test.net1)$weight)*100
V(test.net1)$color <- "black"
V(test.net1)$size <- 3
# L = igraph_layout_mds(test.net1); #Scaling factor s
# autocurve.edges(test.net1)
par(bg = 'black')
plot(test.net1, rescale=T, frame=T,edge.curved=T,# edge.curved=seq(-0.5, 0.5, length = ecount(test.net1)),
     vertex.label.color= "white", edge.label.color= "white", vertex.frame.color="white",
     edge.width=2,#abs(E(test.net1)$weight)*10,
     edge.label=round(E(test.net1)$weight, 2),
     palette=palette(heat.colors(100, alpha=1))) #rainbow(100)))
dev.off()

# sed '1d;$d' nasdaqlisted.txt | awk -F'|' '{print $1}' | tr '\n' '\",\"'