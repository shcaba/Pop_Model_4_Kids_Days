rockfish.ex<-flatfish.ex<-shark.ex<-tilapia

tilapia.2<-tilapia
tilapia.2$M<-0.05

rockfish.ex<-list("Rockfish","growth_VB",0.128749,61.2916,0.240421,65,0.0118,3.094,0.0521,0,1000,"pmat_w",40.5,65,0.2619,0.0217,"knife_edge",40.5,"srrBH")
#names(rockfish.ex<-c("species","K","Linf",))

rockfish.ex$species<-"Rockfish"
rockfish.ex$growthFun<-"growth_VB"
rockfish.ex$K<-0.128749 
rockfish.ex$Linf<-61.2916
rockfish.ex$t0<--0.240421
rockfish.ex$amax<-65
rockfish.ex$LWa<-0.0118
rockfish.ex$LWb<-3.094
rockfish.ex$M<-0.0521
rockfish.ex$F<-0
rockfish.ex$N0<-1000
rockfish.ex$matFun<-"pmat_w"
rockfish.ex$Lmat<- 40.5
rockfish.ex$fec_a<- -1000 #-0.2619
rockfish.ex$fec_b<- 0.55 #0.0217
#rockfish.ex$nspawn<-
rockfish.ex$selectFun<-"knife_edge"
#rockfish.ex$select_p1<-
#rockfish.ex$select_p2<-
rockfish.ex$knife_edge_size<-40.5
#rockfish.ex$mesh_size<-
#rockfish.ex$select_dist<-
#rockfish.ex$mesh_size1<-
rockfish.ex$srrFun<-"srrFecBH"
#rockfish.ex$srrFecBH_a<-
#rockfish.ex$srrFecBH_b<-
rockfish.ex$rmax<-1000
rockfish.ex$beta<-500
rockfish.ex$SB<-500
rockfish.ex


res <- cohortSim(rockfish.ex, t_incr=1)
plot(Neggst ~ Wt, res) # Number of eggs as a function of ind. weight
plot(res$t, res$Neggst*res$pmat*res$Nt) # total population fecundity by age
plot(pcap ~ Lt, res, t="l")
plot(Lt ~ t, res, t="l")
plot(Wt ~ t, res, t="l")

plot(Bt ~ t, res, t="l")
lines(SBt ~ t, res, col=2)

plot(Bt ~ Lt, res, t="l")
lines(SBt ~ Lt, res, col=2)

plot(Yt ~ t, res, t="l")

plot(Nt ~ t, res, t="l", log="y")
lines(Nt.noF ~ t, res, col=2, lty=2)



nyears<-100
#params<-rockfish.ex
Ft <- rep(1.5, nyears)
env_at <- runif(nyears, min=0.5, max=1.5)
env_bt <- rep(1, nyears); env_bt[20:35] <- 0.5
rockfish.sim<- stockSim(Ft=Ft, params=rockfish.ex, nyears=nyears, env_at=env_at, env_bt=env_bt)
plot(rockfish.ex$Bt, t="l")

tilapia.sim<- stockSim(Ft=Ft, params=tilapia, nyears=nyears)
plot(tilapia.sim$Bt, t="l")
tilapia2.sim<- stockSim(Ft=Ft, params=tilapia.2, nyears=nyears)
plot(tilapia2.sim$Bt, t="l")
tilapia.sim<- stockSim(Ft=Ft, params=tilapia, nyears=nyears)
plot(tilapia.sim$Bt, t="l")



years<-100
F.in<-0.05
maxage<-65
Linf=61
K=0.13
t0=0
R0=1000
M = 0.05
ages=c(0:maxage)
a_i=1
steep=0.6
Ltwt.a<-0.0000118
Ltwt.b<-3.094
Mat.a<--0.5
Mat.b<-40
Sel.a<-10
Sel.b<-0.5
  
Sel= HCR.Log.Sel.fit(c(Sel.a,Sel.b),0:65)
Lt.out<-HCR.VBGF(Linf,K,t0,0:maxage)
Maturity<-HCR.Mat.fit(c(Mat.a,Mat.b),Lt.out)
Weight<-HCR.LtWt.fit(c(Ltwt.a,Ltwt.b),Lt.out)

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

