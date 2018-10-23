require(shiny)
require(ggplot2)
source('load_pops.r',local = FALSE)

shinyServer(
  function(input, output) 
  {    
    
    HCR.VBGF<-function(Linf,k,t0,ages)
    {
      Lengths_exp<-Linf*(1-exp(-k*(ages-t0)))
      return(Lengths_exp)
    }
    
    HCR.VBGF.fit<-function(p)
    {
      exp.lts<-VBGF(p[1],p[2],p[3],obs[,1])
      return(exp.lts)
    }
    
    Lt2Age<-function(Linf,k,t0,lt)
    {
      age.lt<-(t0-(log(1-(lt/Linf))/k))
      return(age.lt)
    }
    
    HCR.LtWt.fit<-function(p,Lts)
    {
      exp.wts<-p[1]*(Lts)^p[2]
      return(exp.wts)
    }
    
    HCR.Mat.fit<-function(p,Lts)
    {
      exp.mat<-1/(1+exp(p[1]*(Lts-p[2])))
      return(exp.mat)
    }
    
    HCR.Log.Sel.fit<-function(p,met.vec)
    {
      Exp.sel<-(1/(1+exp(-log(19)*(met.vec-p[1])/p[2])))
      return(Exp.sel)
    }
    
    Numbers<-function(R_in,M,F_yr,Sel,ages,a_i)
    {
      M_step<-M*a_i
      F_step<-F_yr*a_i
      Z_step<-M_step+(F_step*Sel)
      Z<-M+(F_yr*Sel)
      nums<-matrix(data= ages*0, nrow= length(ages), ncol= 2) #Empty matrix for numbers at age per sex
      nums[1,1:2]<-R_in/2		#Initialize female & male columns
      for (i in 1:2)		#Loop to calculate numbers for years 1 to max-1
      {
        for(j in 2:(length(ages)-1))
        {
          nums[j,i]<-nums[j-1,i]*exp(-(Z_step[j-1]))
        }
      }
      
      for (ii in 1:2)	#Loop to calculate numbners in plus group
      {
        nums[length(ages),ii]<-(nums[length(ages)-1,ii]*exp(-(Z_step[length(ages)-1])))/(1-exp(-Z[length(ages)]))
      }
      return(nums)
    }
    
    SpawnB<-function(R_in,Maturity,Weight,M,F_yr,Sel,ages,a_i)
    {
      numbers<-Numbers(R_in,M,F_yr,Sel,ages,a_i)
      SBPRF<-sum(Maturity*Weight*numbers[,1])
      return(SBPRF)
    }
    
    BHRecruit.max<- function(R0,Rt,steep,Weight,Maturity,M,F_yr,Sel,ages,a_i,sigmaR)
    {
      SBPRF0<-SpawnB(R0,Maturity,Weight,M,0,Sel,ages,a_i)
      SBPRF<-SpawnB(Rt,Maturity,Weight,M,F_yr,Sel,ages,a_i)
      RF<-(R0*((4*steep*SBPRF)-(SBPRF0*(1-steep))))/(SBPRF*(5*steep-1))* exp(rnorm(1,0,sigmaR)-sigmaR^2/2)
      return(RF)
    }
    
####### END FUNCTIONS ########
    
    Pop.out.rockfish<- reactive({
      species.in<-rockfish.in
      ages=c(0:species.in$maxage)  
      Sel= HCR.Log.Sel.fit(c(species.in$Sel.a,species.in$Sel.b),ages)
      Lt.out<-HCR.VBGF(species.in$Linf,species.in$K,species.in$t0,ages)
      Maturity<-Maturity<-HCR.Mat.fit(c(species.in$Mat.a,species.in$Mat.b),Lt.out)
      Weight<-HCR.LtWt.fit(c(species.in$Ltwt.a,species.in$Ltwt.b),Lt.out)
      
      Num_all<-Num_F<-Num_M<-matrix(NA,length(ages),input$years.rockfish+1)
      SB.ts<-Rec.ts<-rep(NA,input$years.rockfish+1)
      
      Num<-Numbers(R_in=species.in$R0,species.in$M,F_yr=0,Sel,ages,species.in$a_i)
      Num_all[,1]<-rowSums(Num)
      Num_F[,1]<-Num[,1]
      Num_M[,2]<-Num[,2]
      SB.ts[1]<-SpawnB(R_in=species.in$R0,Maturity,Weight,species.in$M,F_yr=0,Sel,ages,species.in$a_i)
      Rec.ts[1]<-BHRecruit.max(R0=species.in$R0,Rt=species.in$R0,species.in$steep,Weight,Maturity,species.in$M,F_yr=0,Sel,ages,species.in$a_i,sigmaR=0)
      for(i in 1:input$years.rockfish)
      {
        Num<-Numbers(R_in=Rec.ts[i],species.in$M,F_yr=input$F.in.rockfish,Sel,ages,species.in$a_i)
        Num_all[,i+1]<-rowSums(Num)
        Num_F[,i+1]<-Num[,1]
        Num_M[,i+1]<-Num[,2]
        SB.ts[i+1]<-SpawnB(R_in=Rec.ts[i],Maturity,Weight,species.in$M,F_yr=input$F.in.rockfish,Sel,ages,species.in$a_i)
        Rec.ts[i+1]<-BHRecruit.max(R0=species.in$R0,Rt=Rec.ts[i],species.in$steep,Weight,Maturity,species.in$M,F_yr=input$F.in.rockfish,Sel,ages,species.in$a_i,sigmaR=0.25)
      }
      Pop.out.rockfish<-list(Num_all,Num_F,Num_M,SB.ts,Rec.ts)
      names(Pop.out.rockfish)<-c("Numbers_all","Numbers_F","Numbers_M","SB","Rec")
      return(Pop.out.rockfish)
    })
    
     
    Pop.out.flatfish<- reactive({
      species.in<-flatfish.in
      ages=c(0:species.in$maxage)  
      Sel= HCR.Log.Sel.fit(c(species.in$Sel.a,species.in$Sel.b),ages)
      Lt.out<-HCR.VBGF(species.in$Linf,species.in$K,species.in$t0,ages)
      Maturity<-Maturity<-HCR.Mat.fit(c(species.in$Mat.a,species.in$Mat.b),Lt.out)
      Weight<-HCR.LtWt.fit(c(species.in$Ltwt.a,species.in$Ltwt.b),Lt.out)
      
      Num_all<-Num_F<-Num_M<-matrix(NA,length(ages),input$years.flatfish+1)
      SB.ts<-Rec.ts<-rep(NA,input$years.flatfish+1)
      
      Num<-Numbers(R_in=species.in$R0,species.in$M,F_yr=0,Sel,ages,species.in$a_i)
      Num_all[,1]<-rowSums(Num)
      Num_F[,1]<-Num[,1]
      Num_M[,2]<-Num[,2]
      SB.ts[1]<-SpawnB(R_in=species.in$R0,Maturity,Weight,species.in$M,F_yr=0,Sel,ages,species.in$a_i)
      Rec.ts[1]<-BHRecruit.max(R0=species.in$R0,Rt=species.in$R0,species.in$steep,Weight,Maturity,species.in$M,F_yr=0,Sel,ages,species.in$a_i,sigmaR=0)
      for(i in 1:input$years.flatfish)
      {
        Num<-Numbers(R_in=Rec.ts[i],species.in$M,F_yr=input$F.in.flatfish,Sel,ages,species.in$a_i)
        Num_all[,i+1]<-rowSums(Num)
        Num_F[,i+1]<-Num[,1]
        Num_M[,i+1]<-Num[,2]
        SB.ts[i+1]<-SpawnB(R_in=Rec.ts[i],Maturity,Weight,species.in$M,F_yr=input$F.in.flatfish,Sel,ages,species.in$a_i)
        Rec.ts[i+1]<-BHRecruit.max(R0=species.in$R0,Rt=Rec.ts[i],species.in$steep,Weight,Maturity,species.in$M,F_yr=input$F.in.flatfish,Sel,ages,species.in$a_i,sigmaR=0.5)
      }
      Pop.out.flatfish<-list(Num_all,Num_F,Num_M,SB.ts,Rec.ts)
      names(Pop.out.flatfish)<-c("Numbers_all","Numbers_F","Numbers_M","SB","Rec")
      return(Pop.out.flatfish)
    })
    
    Pop.out.shark<- reactive({
      species.in<-shark.in
      ages=c(0:species.in$maxage)  
      Sel= HCR.Log.Sel.fit(c(species.in$Sel.a,species.in$Sel.b),ages)
      Lt.out<-HCR.VBGF(species.in$Linf,species.in$K,species.in$t0,ages)
      Maturity<-Maturity<-HCR.Mat.fit(c(species.in$Mat.a,species.in$Mat.b),Lt.out)
      Weight<-HCR.LtWt.fit(c(species.in$Ltwt.a,species.in$Ltwt.b),Lt.out)
      
      Num_all<-Num_F<-Num_M<-matrix(NA,length(ages),input$years.shark+1)
      SB.ts<-Rec.ts<-rep(NA,input$years.shark+1)
      
      Num<-Numbers(R_in=species.in$R0,species.in$M,F_yr=0,Sel,ages,species.in$a_i)
      Num_all[,1]<-rowSums(Num)
      Num_F[,1]<-Num[,1]
      Num_M[,2]<-Num[,2]
      SB.ts[1]<-SpawnB(R_in=species.in$R0,Maturity,Weight,species.in$M,F_yr=0,Sel,ages,species.in$a_i)
      Rec.ts[1]<-BHRecruit.max(R0=species.in$R0,Rt=species.in$R0,species.in$steep,Weight,Maturity,species.in$M,F_yr=0,Sel,ages,species.in$a_i,sigmaR=0)
      for(i in 1:input$years.shark)
      {
        Num<-Numbers(R_in=Rec.ts[i],species.in$M,F_yr=input$F.in.shark,Sel,ages,species.in$a_i)
        Num_all[,i+1]<-rowSums(Num)
        Num_F[,i+1]<-Num[,1]
        Num_M[,i+1]<-Num[,2]
        SB.ts[i+1]<-SpawnB(R_in=Rec.ts[i],Maturity,Weight,species.in$M,F_yr=input$F.in.shark,Sel,ages,species.in$a_i)
        Rec.ts[i+1]<-BHRecruit.max(R0=species.in$R0,Rt=Rec.ts[i],species.in$steep,Weight,Maturity,species.in$M,F_yr=input$F.in.shark,Sel,ages,species.in$a_i,sigmaR=0.1)
      }
      Pop.out.shark<-list(Num_all,Num_F,Num_M,SB.ts,Rec.ts)
      names(Pop.out.shark)<-c("Numbers_all","Numbers_F","Numbers_M","SB","Rec")
      return(Pop.out.shark)
    })
    
    Pop.out.custom<- reactive({
      species.in<-shark.in
      species.in$maxage<-input$Max_age_cs
      species.in$Linf<-input$Linf_cs
      species.in$K<-input$k_cs
      species.in$t0<-input$t0_cs
      species.in$R0<-input$R0_cs
      species.in$M<-input$M_cs
      species.in$a_i<-input$a_i_cs
      species.in$steep<-input$h_cs
      species.in$Ltwt.a<-input$LW_a_cs
      species.in$Ltwt.b<-input$LW_b_cs
      species.in$Mat.a<-input$Mat_a_cs
      species.in$Mat.b<-input$Mat_L50_cs
      species.in$Sel.a<-input$Sel_a_cs
      species.in$Sel.b<-input$Sel_b_cs
      ages=c(0:species.in$maxage)  
      Sel= HCR.Log.Sel.fit(c(species.in$Sel.a,species.in$Sel.b),ages)
      Lt.out<-HCR.VBGF(species.in$Linf,species.in$K,species.in$t0,ages)
      Maturity<-Maturity<-HCR.Mat.fit(c(species.in$Mat.a,species.in$Mat.b),Lt.out)
      Weight<-HCR.LtWt.fit(c(species.in$Ltwt.a,species.in$Ltwt.b),Lt.out)
      
      Num_all<-Num_F<-Num_M<-matrix(NA,length(ages),input$years.custom+1)
      SB.ts<-Rec.ts<-rep(NA,input$years.custom+1)
      
      Num<-Numbers(R_in=species.in$R0,species.in$M,F_yr=0,Sel,ages,species.in$a_i)
      Num_all[,1]<-rowSums(Num)
      Num_F[,1]<-Num[,1]
      Num_M[,2]<-Num[,2]
      SB.ts[1]<-SpawnB(R_in=species.in$R0,Maturity,Weight,species.in$M,F_yr=0,Sel,ages,species.in$a_i)
      Rec.ts[1]<-BHRecruit.max(R0=species.in$R0,Rt=species.in$R0,species.in$steep,Weight,Maturity,species.in$M,F_yr=0,Sel,ages,species.in$a_i,sigmaR=0)
      for(i in 1:input$years.custom)
      {
        Num<-Numbers(R_in=Rec.ts[i],species.in$M,F_yr=input$F.in.custom,Sel,ages,species.in$a_i)
        Num_all[,i+1]<-rowSums(Num)
        Num_F[,i+1]<-Num[,1]
        Num_M[,i+1]<-Num[,2]
        SB.ts[i+1]<-SpawnB(R_in=Rec.ts[i],Maturity,Weight,species.in$M,F_yr=input$F.in.custom,Sel,ages,species.in$a_i)
        Rec.ts[i+1]<-BHRecruit.max(R0=species.in$R0,Rt=Rec.ts[i],species.in$steep,Weight,Maturity,species.in$M,F_yr=input$F.in.custom,Sel,ages,species.in$a_i,sigmaR=0.1)
      }
      Pop.out.custom<-list(Num_all,Num_F,Num_M,SB.ts,Rec.ts)
      names(Pop.out.custom)<-c("Numbers_all","Numbers_F","Numbers_M","SB","Rec")
      return(Pop.out.custom)
    })
    
       output$rockfish_pop <- renderPlot({
       Pop.out.rockfish<-Pop.out.rockfish()
       for(i in 1:length(Pop.out.rockfish$SB)){if(Pop.out.rockfish$SB[i]<0){Pop.out.rockfish$SB[i:length(Pop.out.rockfish$SB)]=0}}
       SB.gg<-as.data.frame(cbind(c(1:(input$years.rockfish+1)),Pop.out.rockfish$SB))
       names(SB.gg)<-c("Year","SB")
       extict<-SB.gg$SB==0
       rockfish_pop<-ggplot(data=SB.gg,aes(x=Year,y=SB))+
         geom_line(col="orange",lwd=1.5)+
         ylim(0, max(SB.gg$SB)*1.1)+
         geom_line(data=SB.gg[extict,],col="black",lwd=1.5)+
         annotate("text", -Inf, Inf, label = paste("Scale_0=",round(SB.gg$SB[1],0)), hjust = 0, vjust = 1,size=8,color="orange")+
         annotate("text", -Inf, Inf, label = paste("Scale_end=",round(SB.gg$SB[length(SB.gg$SB)],0)), hjust = -1.1, vjust = 1,size=8,color="orange")+
         annotate("text", -Inf, Inf, label = paste("Status=",round(SB.gg$SB[length(SB.gg$SB)]/SB.gg$SB[1],2)), hjust = -3.5, vjust = 1,size=8,color="orange")
       print(rockfish_pop)
     })

       output$flatfish_pop <- renderPlot({
       Pop.out.flatfish<-Pop.out.flatfish()
       for(i in 1:length(Pop.out.flatfish$SB)){if(Pop.out.flatfish$SB[i]<0){Pop.out.flatfish$SB[i:length(Pop.out.flatfish$SB)]=0}}
       SB.gg<-as.data.frame(cbind(c(1:(input$years.flatfish+1)),Pop.out.flatfish$SB))
       names(SB.gg)<-c("Year","SB")
       extict<-SB.gg$SB==0
       flatfish_pop<-ggplot(data=SB.gg,aes(x=Year,y=SB))+
         geom_line(col="green",lwd=1.5)+
         ylim(0, max(SB.gg$SB)*1.1)+
         geom_line(data=SB.gg[extict,],col="black",lwd=1.5)+
         annotate("text", -Inf, Inf, label = paste("Scale_0=",round(SB.gg$SB[1],0)), hjust = 0, vjust = 1,size=8,color="green")+
         annotate("text", -Inf, Inf, label = paste("Scale_end=",round(SB.gg$SB[length(SB.gg$SB)],0)), hjust = -1.1, vjust = 1,size=8,color="green")+
         annotate("text", -Inf, Inf, label = paste("Status=",round(SB.gg$SB[length(SB.gg$SB)]/SB.gg$SB[1],2)), hjust = -3.5, vjust = 1,size=8,color="green")
       print(flatfish_pop)
       })
     
     output$shark_pop <- renderPlot({
       Pop.out.shark<-Pop.out.shark()
       for(i in 1:length(Pop.out.shark$SB)){if(Pop.out.shark$SB[i]<0){Pop.out.shark$SB[i:length(Pop.out.shark$SB)]=0}}
       SB.gg<-as.data.frame(cbind(c(1:(input$years.shark+1)),Pop.out.shark$SB))
       names(SB.gg)<-c("Year","SB")
       extict<-SB.gg$SB==0
       shark_pop<-ggplot(data=SB.gg,aes(x=Year,y=SB))+
         geom_line(col="blue",lwd=1.5)+
         ylim(0, max(SB.gg$SB)*1.1)+
         geom_line(data=SB.gg[extict,],col="black",lwd=1.5)+
         annotate("text", -Inf, Inf, label = paste("Scale_0=",round(SB.gg$SB[1],0)), hjust = 0, vjust = 1,size=8,color="blue")+
         annotate("text", -Inf, Inf, label = paste("Scale_end=",round(SB.gg$SB[length(SB.gg$SB)],0)), hjust = -1.1, vjust = 1,size=8,color="blue")+
         annotate("text", -Inf, Inf, label = paste("Status=",round(SB.gg$SB[length(SB.gg$SB)]/SB.gg$SB[1],2)), hjust = -3.5, vjust = 1,size=8,color="blue")
       print(shark_pop)
     })
     
     output$custom_pop <- renderPlot({
       Pop.out.custom<-Pop.out.custom()
       for(i in 1:length(Pop.out.custom$SB)){if(Pop.out.custom$SB[i]<0){Pop.out.custom$SB[i:length(Pop.out.custom$SB)]=0}}
       SB.gg<-as.data.frame(cbind(c(1:(input$years.custom+1)),Pop.out.custom$SB))
       names(SB.gg)<-c("Year","SB")
       extict<-SB.gg$SB==0
       custom_pop<-ggplot(data=SB.gg,aes(x=Year,y=SB))+
         geom_line(col="purple",lwd=1.5)+
         ylim(0, max(SB.gg$SB)*1.1)+
         geom_line(data=SB.gg[extict,],col="black",lwd=1.5)+
         annotate("text", -Inf, Inf, label = paste("Scale_0=",round(SB.gg$SB[1],0)), hjust = 0, vjust = 1,size=8,color="purple")+
         annotate("text", -Inf, Inf, label = paste("Scale_end=",round(SB.gg$SB[length(SB.gg$SB)],0)), hjust = -1.1, vjust = 1,size=8,color="purple")+
         annotate("text", -Inf, Inf, label = paste("Status=",round(SB.gg$SB[length(SB.gg$SB)]/SB.gg$SB[1],2)), hjust = -3.5, vjust = 1,size=8,color="purple")
       print(custom_pop)
     })
     
  }
)

