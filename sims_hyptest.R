## Simulations for CP hypothesis testing paper

library(dplyr)
library(ggplot2)
library(Rcpp)
library(parallel)
library(ggpubr)
library(igraph)
library(reticulate)
library(RColorBrewer)


source("sims_functions.R")

################## Simulations

myColors <- brewer.pal(6,"Spectral")
myColors[3] <- "darkblue"
myColors[5] <- "darkgreen"


#############################################  
################# ER vs. CL ################# 
############################################# 

# Setting 1 (a): ER vs. CL (generated with CP), vary rho

n.iters=100
delta.seq = seq(0, 0.4, length=11)
p11 = p12 = p22 = 1

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(delta.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(delta.seq)), 
             delta=0, rho=0, reject=0, time=0)

n=1000
cnt=1

for(delta in delta.seq){
  for(iter in 1:n.iters){
    k = 0.1*n
    thetaC <- runif(k, 0.2+delta/2, 0.4+delta/2)
    thetaP = runif(n-k, 0.2-delta/2, 0.4-delta/2)
    
    theta = c(thetaC, thetaP)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 3] <- delta
    df[cnt:(cnt+3), 4] <- rho

    
    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(delta)
  #save(df, file="sims_er_cl_hyptest_rho_042425.RData")
}

#load("sims_er_cl_hyptest_rho_042425.RData")


df_plot <- df %>% group_by(delta, method) %>% summarise(reject=mean(reject), rho=mean(rho), time=mean(time))
df_plot <- df_plot[df_plot$method=="Proposed", ]


p1 <- ggplot(df_plot, aes(x=rho, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_er_cl_rho_042425.pdf", device="pdf", width=6, height=4, units="in")


# Setting 1 (b): ER vs. CL, vary n


n.iters=100
n.seq = seq(1000, 2000, 200)

delta = 0.22

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(n.seq)), 
             n=0, rho=0, reject=0, time=0)

cnt=1

for(n in n.seq){
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- n
    
    k = 0.1*n
    thetaC <- runif(k, 0.2+delta/2, 0.4+delta/2)
    thetaP = runif(n-k, 0.2-delta/2, 0.4-delta/2)
    
    theta = c(thetaC, thetaP)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 4] <- rho

    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(n)
  #save(df, file="sims_er_cl_hyptest_n_042425.RData")
}

#load("sims_er_cl_hyptest_n_042425.RData")

df <- df[!rowSums(is.na(df)), ]


df_plot <- df %>% group_by(n, method) %>% summarise(reject=mean(reject), rho=mean(rho))
df_plot <- df_plot[df_plot$method=="Proposed", ]

p1 <- ggplot(df_plot, aes(x=n, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  ylim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_er_cl_n_042425.pdf", device="pdf", width=6, height=4, units="in")


# Setting 1 (c): ER vs. CL, not disparate enough theta

# theta[k]*theta[n]
# theta[k+1]*theta[k+2]

n.iters=100
delta.seq = seq(0, 0.4, length=11)
n=1000

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(delta.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(delta.seq)), 
             delta=0, rho=0, reject=0, time=0, thetakn=0, thetak1k2=0)

cnt=1

for(delta in delta.seq){
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- delta
    
    k = 0.1*n
    theta <- runif(n, 0.3-delta/2, 0.3+delta/2)
    theta <- sort(theta, decreasing = TRUE)
    
    df[cnt+3, 7] = theta[k]   * theta[n]
    df[cnt+3, 8] = theta[k+1] * theta[k+2]
    
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 4] <- rho
    
    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(delta)
  #save(df, file="sims_er_cl_hyptest_null_n_042425.RData")
}

load("sims_er_cl_hyptest_null_n_042425.RData")

df <- df[!is.na(df$reject), ]
df_plot <- df %>% group_by(delta, method) %>% summarise(reject=mean(reject), rho=mean(rho),
                                                        thetakn = mean(thetakn),
                                                        thetak1k2 = mean(thetak1k2))
df_plot <- df_plot[df_plot$method=="Proposed", ]

p1 <- ggplot(df_plot, aes(x=rho, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  ylim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_er_cl_null_042425.pdf", device="pdf", width=6, height=4, units="in")

# Table of theta[k] * theta[n] and theta[k+1] * theta[k+2]
df_plot

round(df_plot$thetakn,2)
round(df_plot$thetak1k2,2)


#############################################  
################# ER vs. SBM ################
############################################# 

# Setting 2 (a): ER vs. SBM (generated with CP), vary p

n.iters=100
n=1000

p12.seq = seq(0.01, 0.05, length=10)
p22=0.01

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(p12.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(p12.seq)), 
             rho=0, reject=0, time=0)

cnt=1

for(p12 in p12.seq){
  p11=2*p12
  
  Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
  P <- generateP(n, p11, p12, p22, prop=0.1)
  rho <- rho.fun(P, Cstar)
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- rho
    
    A <- generateA(n, p11, p12, p22, prop=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    # start = proc.time()[3]
    # df[cnt, 4]   <- boot_boyd(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+1, 4]   <- boot_deg(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+1, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+2, 4]   <- boot_ER(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+2, 5] = end - start
    
    start = proc.time()[3]
    df[cnt+3, 4] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 5] = end - start
    
    cnt = cnt + 4
  }
  
  print(p12)
  #save(df, file="sims_sbm_hyptest_rho_040924.RData")
}


load("sims_sbm_hyptest_rho_040924.RData")


df_plot <- df %>% group_by(rho, method) %>% summarise(reject=mean(reject), time=mean(time))

df_plot = df_plot[df_plot$method=="Proposed", ]

p1 <- ggplot(df_plot, aes(x=rho, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  xlab("rho")+
  theme_bw()+
  theme(text = element_text(size = 20))
p1

#ggsave(file = "sims_sbm_rho_0409224.pdf", device="pdf", width=6, height=4, units="in")



# Setting 2 (b): ER vs. SBM (generated with CP), vary n

n.iters=100
n.seq = seq(1000, 2500, 250)


df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(n.seq)), 
             n=0, reject=0, time=0)

cnt=1

for(n in n.seq){
  
  p12 = 0.02
  p11 = 0.05
  p22 = 0.01
  
  Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
  P <- generateP(n, p11, p12, p22, prop=0.1)
  rho <- rho.fun(P, Cstar)
  print(rho)
  
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- n
    
    A <- generateA(n, p11, p12, p22, prop=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    # start = proc.time()[3]
    # df[cnt, 4]   <- boot_boyd(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+1, 4]   <- boot_deg(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+1, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+2, 4]   <- boot_ER(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+2, 5] = end - start
    
    start = proc.time()[3]
    df[cnt+3, 4] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 5] = end - start
    
    cnt = cnt + 4
    #print(iter)
  }
  
  print(n)
  #save(df, file="sims_sbm_hyptest_n_040924.RData")
}

load("sims_sbm_hyptest_n_040924.RData")

df_plot <- df %>% group_by(n, method) %>% summarise(reject=mean(reject), time=mean(time))
df_plot <- df_plot[df_plot$method=="Proposed", ]

p2 <- ggplot(df_plot, aes(x=n, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  theme_bw()+
  theme(text = element_text(size = 20))
p2

#ggsave(file = "sims_sbm_n_040924.pdf", device="pdf", width=6, height=4, units="in")




# Setting 2 (c): ER vs. SBM (not generated with CP)

n.iters=100
n=1000

p12.seq = seq(0.01, 0.10, length=5)
p22=0.01

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(p12.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(p12.seq)), 
             rho=0, reject=0, time=0)

cnt=1

for(p12 in p12.seq){
  p11=p22
  
  for(iter in 1:n.iters){
    
    A <- generateA(n, p11, p12, p22, prop=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    P <- generateP(n, p11, p12, p22, prop=0.1)
    df[cnt:(cnt+3), 3] <- rho.fun(P, C)
    
    # start = proc.time()[3]
    # df[cnt, 4]   <- boot_boyd(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+1, 4]   <- boot_deg(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+1, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+2, 4]   <- boot_ER(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+2, 5] = end - start
    
    start = proc.time()[3]
    df[cnt+3, 4] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 5] = end - start
    
    cnt = cnt + 4
    #print(iter)
  }
  df[(cnt  - 4*n.iters):(cnt - 1) , 3] <- mean(df$rho[(cnt - 4*n.iters):(cnt - 1)])
  print(p12)
  #save(df, file="sims_sbm_hyptest_null_040924.RData")
}

load("sims_sbm_hyptest_null_040924.RData")

df_plot <- df %>% group_by(rho, method) %>% summarise(reject=mean(reject), time=mean(time))


p3 <- ggplot(df_plot, aes(x=rho, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  ylim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 20))
p3



#ggsave(file = "sims_dis_040924.pdf", device="pdf", width=6, height=4, units="in")


#############################################  
################# ER vs. DCBM ############### 
#############################################

# Setting 3 (a): ER vs. DCBM (generated with CP), vary p

n.iters=100
p.seq = seq(0.05, 0.25, length=10)


df <- tibble(iter = rep(rep(1:n.iters,each=4), length(p.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(p.seq)), 
             p12=0, rho=0, reject=0, time=0)

n=1000
cnt=1

for(p12 in p.seq){
  
  p11 = 2*p12
  p22 = 0.05
  
  
  for(iter in 1:n.iters){
    
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)
    
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 3] <- p12
    df[cnt:(cnt+3), 4] <- rho
    
    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(p12)
  #save(df, file="sims_er_dcbm_hyptest_rho_042425.RData")
}

#load("sims_er_dcbm_hyptest_rho_042425.RData")


df_plot <- df %>% group_by(p12, method) %>% summarise(reject=mean(reject), rho=mean(rho), time=mean(time))
df_plot <- df_plot[df_plot$method=="Proposed", ]


p1 <- ggplot(df_plot, aes(x=rho, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_er_dcbm_rho_042425.pdf", device="pdf", width=6, height=4, units="in")


# Setting 3 (b): ER vs. DCBM, wtih CP, vary n


n.iters=100
n.seq = seq(1000, 2000, 200)


df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(n.seq)), 
             n=0, rho=0, reject=0, time=0)

p11=0.24
p12=0.12
p22=0.05
cnt=1

for(n in n.seq){
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- n
    
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)

    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 4] <- rho
    
    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(n)
  #  save(df, file="sims_er_dcbm_hyptest_n_042425.RData")
}

#load("sims_er_dcbm_hyptest_n_042425.RData")

df <- df[!rowSums(is.na(df)), ]


df_plot <- df %>% group_by(n, method) %>% summarise(reject=mean(reject), rho=mean(rho))
df_plot <- df_plot[df_plot$method=="Proposed", ]

p1 <- ggplot(df_plot, aes(x=n, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  ylim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_er_dcbm_n_042425.pdf", device="pdf", width=6, height=4, units="in")


# Setting 3 (c): ER vs. DCBM, wtih disassortative, vary n


n.iters=100
n.seq = seq(1000, 5000, 1000)


df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(n.seq)), 
             n=0, rho=0, reject=0, time=0)

p11=0.05
p12=0.10
p22=0.05
cnt=1

for(n in n.seq){
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- n
    
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 4] <- rho
    
    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "ER")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(n)
  #save(df, file="sims_er_dcbm_hyptest_null_n_042425.RData")
}

load("sims_er_dcbm_hyptest_null_n_042425.RData")


df_plot <- df %>% group_by(n, method) %>% summarise(reject=mean(reject), rho=mean(rho))
df_plot <- df_plot[df_plot$method=="Proposed", ]

p1 <- ggplot(df_plot, aes(x=n, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  ylim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_er_dcbm_null_042425.pdf", device="pdf", width=6, height=4, units="in")


#############################################  
################# CL vs. DCBM ###############
############################################# 

# Setting 4 (a): CL vs. DCBM (generated with CP), vary p

n.iters=100
p.seq = seq(0.05, 0.25, length=10)


df <- tibble(iter = rep(rep(1:n.iters,each=4), length(p.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(p.seq)), 
             p12=0, rho=0, reject=0, time=0)

n=1000
cnt=1

for(p12 in p.seq){
  
  p11 = 2*p12
  p22 = 0.05
  
  
  for(iter in 1:n.iters){
    
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)
    #theta[1:k] <- theta[1:k]/sum(theta[1:k])*k
    #theta[(k+1):n] <- theta[(k+1):n]/sum(theta[(k+1):n])*(n-k)
    
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 3] <- p12
    df[cnt:(cnt+3), 4] <- rho
    
    # start = proc.time()[3]
    # df[cnt, 4]   <- boot_boyd(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+1, 4]   <- boot_deg(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+1, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+2, 4]   <- boot_ER(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+2, 5] = end - start
    
    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "CL")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(p12)
  #save(df, file="sims_dcbm_hyptest_rho_040924.RData")
}

#load("sims_dcbm_hyptest_rho_040924.RData")


df_plot <- df %>% group_by(p12, method) %>% summarise(reject=mean(reject), rho=mean(rho), time=mean(time))
df_plot <- df_plot[df_plot$method=="Proposed", ]


p1 <- ggplot(df_plot, aes(x=rho, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_dcbm_rho_040924.pdf", device="pdf", width=6, height=4, units="in")


# Setting 4 (b): CL vs. DCBM, wtih CP, vary n


n.iters=100
n.seq = seq(1000, 2000, 200)


df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(n.seq)), 
             n=0, rho=0, reject=0, time=0)

p11=0.24
p12=0.12
p22=0.05
cnt=1

for(n in n.seq){
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- n
    
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)
    #theta[1:k] <- theta[1:k]/sum(theta[1:k])*k
    #theta[(k+1):n] <- theta[(k+1):n]/sum(theta[(k+1):n])*(n-k)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 4] <- rho
    
    # start = proc.time()[3]
    # df[cnt, 4]   <- boot_boyd(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+1, 4]   <- boot_deg(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+1, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+2, 4]   <- boot_ER(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+2, 5] = end - start
    
    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "CL")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(n)
  #  save(df, file="sims_dcbm_hyptest_n_040924.RData")
}

#load("sims_dcbm_hyptest_n_040924.RData")

df <- df[!rowSums(is.na(df)), ]


df_plot <- df %>% group_by(n, method) %>% summarise(reject=mean(reject), rho=mean(rho))
df_plot <- df_plot[df_plot$method=="Proposed", ]

p1 <- ggplot(df_plot, aes(x=n, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  ylim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_dcbm_n_040924.pdf", device="pdf", width=6, height=4, units="in")


# Setting 4 (c): CL vs. DCBM, wtih disassortative, vary n


n.iters=100
n.seq = seq(1000, 5000, 1000)


df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("Boyd", "Holme", "ER", "Proposed"), n.iters*length(n.seq)), 
             n=0, rho=0, reject=0, time=0)

p11=0.05
p12=0.10
p22=0.05
cnt=1

for(n in n.seq){
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- n
    
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)
    #theta[1:k] <- theta[1:k]/sum(theta[1:k])*k
    #theta[(k+1):n] <- theta[(k+1):n]/sum(theta[(k+1):n])*(n-k)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    C <- borgattiCpp(A)
    tstat <- obj.fun(A, C)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    df[cnt:(cnt+3), 4] <- rho
    
    # start = proc.time()[3]
    # df[cnt, 4]   <- boot_boyd(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+1, 4]   <- boot_deg(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+1, 5] = end - start
    # 
    # start = proc.time()[3]
    # df[cnt+2, 4]   <- boot_ER(A, tstat, 200, 6)
    # end = proc.time()[3]
    # df[cnt+2, 5] = end - start
    
    start = proc.time()[3]
    df[cnt+3, 5] <- anal_rejc(A, C, "CL")$reject
    end = proc.time()[3]
    df[cnt+3, 6] = end - start
    
    cnt = cnt + 4
  }
  
  print(n)
  #save(df, file="sims_dcbm_hyptest_null_n_040924.RData")
}

load("sims_dcbm_hyptest_null_n_040924.RData")


df_plot <- df %>% group_by(n, method) %>% summarise(reject=mean(reject), rho=mean(rho))
df_plot <- df_plot[df_plot$method=="Proposed", ]

p1 <- ggplot(df_plot, aes(x=n, y=reject))+
  geom_point(size=3, color=myColors[1])+
  geom_line(linewidth=1.2, color=myColors[1])+
  ylab("Rejection rate")+
  ylim(0,1)+
  theme_bw()+
  theme(text = element_text(size = 20))
p1


#ggsave(file = "sims_dcbm_null_040924.pdf", device="pdf", width=6, height=4, units="in")
