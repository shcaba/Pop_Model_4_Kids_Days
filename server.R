require(shiny)
require(ggplot2)

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
    
    BHRecruit.max<- function(R0,steep,Weight,Maturity,M,F_yr,Sel,ages,a_i)
    {
      SBPRF0<-SpawnB(R0,Maturity,Weight,M,0,Sel,ages,a_i)
      SBPRF<-SpawnB(R0,Maturity,Weight,M,F_yr,Sel,ages,a_i)
      RF<-(R0*((4*steep*SBPRF)-(SBPRF0*(1-steep))))/(SBPRF*(5*steep-1))* exp(rnorm(1,0,0.5)-0.5^2/2)
      return(RF)
    }
    
    
####### END FUNCTIONS ########
    
    
 Pop_out<- reactive({
   Lt.out<-HCR.VBGF(Linf,K,t0,0:maxage)
   Maturity<-HCR.Mat.fit(c(-0.5,40),Lt.out)
   Weight<-HCR.LtWt.fit(c(0.0000118,3.094),Lt.out)
   
   Num_all<-Num_F<-Num_M<-matrix(NA,length(ages),years+1)
   SB.ts<-Rec.ts<-rep(NA,years+1)
   
   Num<-Numbers(R_in=R0,M,F_yr=0,Sel,ages,a_i)
   Num_all[,1]<-rowSums(Num)
   Num_F[,1]<-Num[,1]
   Num_M[,2]<-Num[,2]
   SB.ts[1]<-SpawnB(R_in=R0,Maturity,Weight,M,F_yr=0,Sel,ages,a_i)
   Rec.ts[1]<-BHRecruit.max(R0,steep,Weight,Maturity,M,F_yr=0,Sel,ages,a_i)

   for(i in 1:years)
   {
     Num<-Numbers(R_in=Rec.ts[i],M,F_yr=F.in,Sel,ages,a_i)
     Num_all[,i+1]<-rowSums(Num)
     Num_F[,i+1]<-Num[,1]
     Num_M[,i+1]<-Num[,2]
     SB.ts[i+1]<-SpawnB(R_in=Rec.ts[i],Maturity,Weight,M,F_yr=F.in,Sel,ages,a_i)
     Rec.ts[i+1]<-BHRecruit.max(R0,steep,Weight,Maturity,M,F_yr=F.in,Sel,ages,a_i)
   }
    Pop.out<-list(Num_all,Num_F,Num_M,SB.ts,Rec.ts)
    names(Pop.out)<-c("Numbers_all","Numbers_F","Numbers_M","SB","Rec")
    Pop.out
   })
     
 output$Mplot <- renderPlot({
   M_vals_all<-M_vals_all()
   M_methods<-c("Then_Amax 1","Then_Amax 2","Then_Amax 3","Hamel_Amax","AnC","Then_VBGF","Jensen_VBGF 1","Jensen_VBGF 2","Pauly_lt","Gislason","Chen-Wat","Roff","Jensen_Amat","Pauly_wt","PnW","Lorenzen","GSI")
  # plot M
   if(all(is.na(M_vals_all))){ymax<-0.5}
   if(!(all(is.na(M_vals_all)))){ymax<-ceiling((max(M_vals_all,na.rm=TRUE)*1.1*10))/10}
   par(mar=c(8,4,2,6),xpd =TRUE)
   plot(M_vals_all, col = "black",bg=c("blue","blue","blue","blue","green","green","green","green","yellow","yellow","orange","red","red","black","black","black","purple"),xlab=" ",ylab="Natural mortality",ylim=c(0,ymax),pch=22,cex=1.5,axes=F)
   box()
   axis(1,at=1:length(M_vals_all),labels=M_methods,las=3)
   axis(2)
   legend(x="topright",legend=c("Amax","VBGF","VBGF:Temp","VBGF;Amat","Amat","Weight","GSI"),pch=22,col="black",pt.bg=c("blue","green","yellow","orange","red","black","purple"),bty="n",horiz=FALSE,cex=1,inset=c(-0.125,0))
   M_table<-data.frame(cbind(M_methods,M_vals_all))
   colnames(M_table)<-c("Method","M")
  # if(all(is.na(M_vals()))){return(NULL)}
   output$downloadMs <- downloadHandler(
    filename = function() {paste0("M_values", '.csv') },
    content = function(file) {write.csv(M_table, file=file)}
  )
 })


#Plot Composite M
 output$Mcomposite<- renderPlot({    
   if(all(is.na(M_vals_all()))){return(NULL)}
   else{
   M.wts<-c(input$Then_Amax_1,input$Then_Amax_2,input$Then_Amax_3,input$Hamel_Amax,input$AnC,input$Then_VBGF,input$Jensen_VBGF_1,input$Jensen_VBGF_2,input$Pauly_lt,input$Gislason,input$Chen_Wat,input$Roff,input$Jensen_Amat,input$Pauly_wt,input$PnW,input$Lorenzen,input$Gonosoma)
   #remove NAs
   if(any(is.na(M_vals_all()))){
     NA.ind<-attributes(na.omit(M_vals_all()))$na.action
     M.sub<-M_vals_all()[-NA.ind]
     M.wts.sub<-M.wts[-NA.ind]
   }
   else{
     M.sub<-M_vals_all()
     M.wts.sub<-M.wts
   }
   #remove 0 weight
   M.sub.n0<-M.sub[M.wts.sub>0]
   M.wts.sub.n0<-M.wts.sub[M.wts.sub>0]
   M.wts.sub.stand<-M.wts.sub.n0/sum(M.wts.sub.n0)
   M.densum<-density(M.sub.n0,weights=M.wts.sub.stand,from=0,cut=0)
   #Approximate the denisty function
   f<- approxfun(M.densum$x, M.densum$y, yleft=0, yright=0)
   #Standardize densities
   pdf_counts<-round(100000*(M.densum$y/sum(M.densum$y)))
   #Expand densities to samples
   pdf.samples<-unlist(mapply(rep,M.densum$x,pdf_counts))
   #Calculate the cdf
   cdf.out<-ecdf(pdf.samples)
   #Plot the density function
   M.densum.plot<- data.frame(x = M.densum$x, y = M.densum$y)
   Mcomposite.densityplot<- ggplot(data=M.densum.plot,aes(x,y,fill="blue"))+
     geom_line(col="black")+
     labs(x="Natural Mortality",y="Density")+ 
     geom_area(fill="gray")+ 
     #scale_x_continuous(limits=c(0,quantile(M.densum$x,0.99)))+
     geom_vline(xintercept = quantile(cdf.out,0.5),color="darkblue",size=1.2)
   print(Mcomposite.densityplot)
   output$downloadMcompositedensityplot <- downloadHandler(
   filename = function() { paste0('Mcomposite_densityplot',Sys.time(), '.png')},
   content = function(file) {
     png(file, type='cairo',width=800,height=720)
     print(Mcomposite.densityplot)
     dev.off()},contentType = 'image/png') 
   output$downloadMcompositedist <- downloadHandler(
     filename = function() {  paste0("Mcomposite_samples",Sys.time(),".DMP") },
     content = function(file) {save(pdf.samples,file=file)}) 
    }
   })
  }
)

