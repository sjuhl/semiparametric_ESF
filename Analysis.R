##########################
# Analysis Script
##########################

# clean environment
rm(list=ls())

# read in simulation putput
load("./SimOut/MC_Out.RData")

# how long did the simulations take
time.taken

# session info (for Monte Carlos)
session_info

# split SAR and SEM DGPs
sar_out <- sim_out[!input$SEM,]
sem_out <- sim_out[input$SEM,]
sar_input <- input[!input$SEM,!(colnames(input) %in% "SEM")]
sem_input <- input[input$SEM,!(colnames(input) %in% "SEM")]


# residual autocorrelation
cols <- colnames(sim_out)[grepl("moran_+",colnames(sim_out))]
pts <- 21:25
lins <- 1:length(cols)
layout(matrix(c(1,2,3),ncol=3,byrow=F))
par(oma=c(3,3.5,2,3.5),mar=c(.3,.3,.3,.3))
for(id in unique(input$W_id)){
  sel <- input$SEM==TRUE & input$W_id==id
  plot(0,xlim=c(0,.75),ylim=c(-.35,8),type="n",axes=F,ann=F)
  abline(h=0,lwd=2,col="red")
  if (id==unique(input$W_id)[1]) axis(2)
  axis(1,at=unique(input$p))
  for(c in cols){
    lines(x=unique(input$p),y=c(median(sim_out[sel & input$p==0,c])
                                ,median(sim_out[sel & input$p==.25,c])
                                ,median(sim_out[sel & input$p==.5,c])
                                ,median(sim_out[sel & input$p==.75,c]))
          ,type="b",pch=pts[match(c,cols)],lty=lins[match(c,cols)])
  }
}

