## CP identification simulations
setwd("~/Documents/Research/Srijan/CP_clean/Code")


library(dplyr)
library(ggplot2)
library(Rcpp)
library(parallel)
library(ggpubr)
library(igraph)
library(reticulate)
library(RColorBrewer)

#py_install("cpnet")

source("sims_functions.R")
source_python('cp_lip.py')
source_python('cp_rombach.py')
source_python('cp_surprise.py')
source_python('cp_cur.py')
source_python('cp_be.py')

################## Simulations

myColors <- brewer.pal(6,"Spectral")
myColors[3] <- "darkblue"
myColors[5] <- "darkgreen"

# Setting 0 (a): CL with CP, vary p

n.iters=100
n=1000

delta.seq = seq(0, 0.2, length=11)
p11 = p12 = p22 = 1

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(delta.seq)), 
             method=rep(c("BE", "Rombach", "Curcuringu", "BayesSBM"), n.iters*length(delta.seq)), 
             rho=0, class=0, time=0, delta=0)

cnt=1

for(delta in delta.seq){
  
  for(iter in 1:n.iters){
    k = 0.1*n
    thetaC <- runif(k, 0.1+delta/2, 0.2+delta/2)
    thetaP = runif(n-k, 0.1-delta/2, 0.2-delta/2)
    
    theta = c(thetaC, thetaP)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    
    
    df[cnt:(cnt+3), 3] <- rho
    df[cnt:(cnt+3), 6] <- delta
    
    # Borgatti and Everett
    start = proc.time()[3]
    C <- borgattiCpp(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    
    # Rombach
    start = proc.time()[3]
    C <- C_rombach(A, k)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start 
    
    # Curcuringu
    start = proc.time()[3]
    C <- C_cur(A)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start  
    
    # BayesSBM
    start = proc.time()[3]
    C <- C_SBM(A, 100, 50)
    end = proc.time()[3]
    df[cnt+3, 4] <- class_acc(C, Cstar)
    df[cnt+3, 5] <- end - start  
    
    cnt = cnt + 4
    save(df, file = "sims_cl_class_rho_042325.RData")
  }
  
  print(delta)
}


load("sims_cl_class_rho_042325.RData")

df_plot <- df %>% group_by(delta, method) %>% summarise(rho = mean(rho), accuracy = mean(class), time = mean(time))

df_plot$method[df_plot$method=="BE"] <- "Proposed"

df_plot$method <- factor(df_plot$method, levels = c("Proposed", "BayesSBM", "Curcuringu", "Rombach"))


p1 <- ggplot(df_plot, aes(x=rho, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()


p2 <- ggplot(df_plot, aes(x=rho, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylab("Time (sec)")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

#ggsave(file = "sims_cl_class_rho_042325.pdf", device="pdf", width=6, height=4, units="in")


# Setting 0 (b): CL with CP, vary n

n.iters=100
n.seq = seq(500, 2000, 250)
delta = 0.16

p11 = p12 = p22 = 1

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("BE", "Rombach", "Curcuringu", "BayesSBM"), n.iters*length(n.seq)), 
             n=0, class=0, time=0, rho=0)

cnt=1

for(n in n.seq){
  
  for(iter in 1:n.iters){
    k = 0.1*n
    thetaC <- runif(k, 0.1+delta/2, 0.2+delta/2)
    thetaP = runif(n-k, 0.1-delta/2, 0.2-delta/2)
    theta = c(thetaC, thetaP)
    
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    
    
    df[cnt:(cnt+3), 3] <- n
    df[cnt:(cnt+3), 6] <- rho
    
    # Borgatti and Everett
    start = proc.time()[3]
    C <- borgattiCpp(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    
    # Rombach
    start = proc.time()[3]
    C <- C_rombach(A, k)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start 
    
    # Curcuringu
    start = proc.time()[3]
    C <- C_cur(A)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start  
    
    # BayesSBM
    start = proc.time()[3]
    C <- C_SBM(A, 100, 50)
    end = proc.time()[3]
    df[cnt+3, 4] <- class_acc(C, Cstar)
    df[cnt+3, 5] <- end - start  
    
    cnt = cnt + 4
    save(df, file = "sims_cl_class_n_042325.RData")
  }
  
  print(n)
}


load("sims_cl_class_n_042325.RData")

df_plot <- df %>% group_by(n, method) %>% summarise(accuracy = mean(class), time = mean(time), rho=mean(rho))

df_plot$method[df_plot$method=="BE"] <- "Proposed"

df_plot$method <- factor(df_plot$method, levels = c("Proposed", "BayesSBM", "Curcuringu", "Rombach"))


p1 <- ggplot(df_plot, aes(x=n, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

p2 <- ggplot(df_plot, aes(x=n, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylab("Time (sec)")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

#ggsave(file = "sims_cl_class_n_042325.pdf", device="pdf", width=6, height=4, units="in")



#######################################################

# Setting 1 (a): SBM with CP, vary p

n.iters=100
n=1000

p12.seq = seq(0.002, 0.026, length=11)
p22=0.001

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(p12.seq)), 
             method=rep(c("BE", "Rombach", "Curcuringu", "BayesSBM"), n.iters*length(p12.seq)), 
             rho=0, class=0, time=0)

cnt=1

for(p12 in p12.seq){
  p11=2*p12
  
  k = 0.1*n
  Cstar = c(rep(1,k), rep(0, n-k)) 
  P <- generateP(n, p11, p12, p22, prop=k/n)
  rho <- rho.fun(P, Cstar)
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- rho
    
    A <- generateA(n, p11, p12, p22, prop=k/n)
    
    # Borgatti and Everett
    start = proc.time()[3]
    C <- borgattiCpp(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    
    # Rombach
    start = proc.time()[3]
    C <- C_rombach(A, k)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start 
    
    # Curcuringu
    start = proc.time()[3]
    C <- C_cur(A)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start  
    
    # BayesSBM
    start = proc.time()[3]
    C <- C_SBM(A, 100, 50)
    end = proc.time()[3]
    df[cnt+3, 4] <- class_acc(C, Cstar)
    df[cnt+3, 5] <- end - start  
    
    cnt = cnt + 4
    #save(df, file = "sims_sbm_class_rho_041824.RData")
  }
  
  print(p12)
}


load("sims_sbm_class_rho_041824.RData")

df_plot <- df %>% group_by(rho, method) %>% summarise(accuracy = mean(class), time = mean(time))

df_plot$method[df_plot$method=="BE"] <- "Proposed"

df_plot$method <- factor(df_plot$method, levels = c("Proposed", "BayesSBM", "Curcuringu", "Rombach"))


p1 <- ggplot(df_plot, aes(x=rho, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()
p1


p2 <- ggplot(df_plot, aes(x=rho, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylab("Time (sec)")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

#ggsave(file = "sims_sbm_class_rho_041824.pdf", device="pdf", width=6, height=4, units="in")






# Setting 1 (b): SBM with CP, vary n

n.iters=100
n.seq=seq(500, 2500, 250)

p22 = 0.001
p12 = 0.0075
p11 = 2*p12

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("BE", "Rombach", "Curcuringu", "BayesSBM"), n.iters*length(n.seq)), 
             n=0, class=0, time=0)

cnt=1

for(n in n.seq){
  
  k = 0.1*n
  Cstar = c(rep(1,k), rep(0, n-k)) 
  P <- generateP(n, p11, p12, p22, prop=k/n)
  rho <- rho.fun(P, Cstar)
  
  for(iter in 1:n.iters){
    df[cnt:(cnt+3), 3] <- n
    
    A <- generateA(n, p11, p12, p22, prop=k/n)
    
    # Borgatti and Everett
    start = proc.time()[3]
    C <- borgattiCpp(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    
    # Rombach
    start = proc.time()[3]
    C <- C_rombach(A, k)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start 
    
    # Curcuringu
    start = proc.time()[3]
    C <- C_cur(A)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start  
    
    # BayesSBM
    start = proc.time()[3]
    C <- C_SBM(A, 100, 50)
    end = proc.time()[3]
    df[cnt+3, 4] <- class_acc(C, Cstar)
    df[cnt+3, 5] <- end - start  
    
    cnt = cnt + 4
    save(df, file = "sims_sbm_class_n_041824.RData")
  }
  
  print(n)
}

load("sims_sbm_class_n_041824.RData")

df_plot <- df %>% group_by(n, method) %>% summarise(accuracy = mean(class), time = mean(time))

df_plot$method[df_plot$method=="BE"] <- "Proposed"

df_plot$method <- factor(df_plot$method, levels = c("Proposed", "BayesSBM", "Curcuringu", "Rombach"))


p1 <- ggplot(df_plot, aes(x=n, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()
p1


p2 <- ggplot(df_plot, aes(x=n, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylab("Time (sec)")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

#ggsave(file = "sims_sbm_class_n_041824.pdf", device="pdf", width=6, height=4, units="in")



# Setting 2 (a): DCBM with CP, vary p

n.iters=100
n=1000

p12.seq = seq(0.05, 0.15, length=11)
p22=0.05

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(p12.seq)), 
             method=rep(c("BE", "Rombach", "Curcuringu", "BayesSBM"), n.iters*length(p12.seq)), 
             rho=0, class=0, time=0, p12=0)

cnt=1

for(p12 in p12.seq){
  p11=2*p12
  
  for(iter in 1:n.iters){
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)
    #theta[1:k] <- theta[1:k]/sum(theta[1:k])*k
    #theta[(k+1):n] <- theta[(k+1):n]/sum(theta[(k+1):n])*(n-k)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    
    
    df[cnt:(cnt+3), 3] <- rho
    df[cnt:(cnt+3), 6] <- p12
    
    # Borgatti and Everett
    start = proc.time()[3]
    C <- borgattiCpp(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    
    # Rombach
    start = proc.time()[3]
    C <- C_rombach(A, k)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start 
    
    # Curcuringu
    start = proc.time()[3]
    C <- C_cur(A)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start  
    
    # BayesSBM
    start = proc.time()[3]
    C <- C_SBM(A, 100, 50)
    end = proc.time()[3]
    df[cnt+3, 4] <- class_acc(C, Cstar)
    df[cnt+3, 5] <- end - start  
    
    cnt = cnt + 4
    #save(df, file = "sims_dcbm_class_rho_041824.RData")
  }
  
  print(p12)
}


load("sims_dcbm_class_rho_041824.RData")

df_plot <- df %>% group_by(p12, method) %>% summarise(rho = mean(rho), accuracy = mean(class), time = mean(time))

df_plot$method[df_plot$method=="BE"] <- "Proposed"

df_plot$method <- factor(df_plot$method, levels = c("Proposed", "BayesSBM", "Curcuringu", "Rombach"))


p1 <- ggplot(df_plot, aes(x=rho, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()
p1


p2 <- ggplot(df_plot, aes(x=rho, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylab("Time (sec)")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

#ggsave(file = "sims_dcbm_class_rho_041824.pdf", device="pdf", width=6, height=4, units="in")


# Setting 2 (b): DCBM with CP, vary n

n.iters=100
n.seq = seq(500, 2000, 250)

p22 = 0.05
p12 = 0.10
p11 = 2*p12

df <- tibble(iter = rep(rep(1:n.iters,each=4), length(n.seq)), 
             method=rep(c("BE", "Rombach", "Curcuringu", "BayesSBM"), n.iters*length(n.seq)), 
             n=0, class=0, time=0, rho=0)

cnt=1

for(n in n.seq){
  
  for(iter in 1:n.iters){
    k = 0.1*n
    theta <- runif(n, 0.6, 0.8)
    # theta[1:k] <- theta[1:k]/sum(theta[1:k])*k
    # theta[(k+1):n] <- theta[(k+1):n]/sum(theta[(k+1):n])*(n-k)
    A <- generateDCBM(theta, p11, p12, p22, q=0.1)
    
    Cstar = c(rep(1,0.1*n), rep(0, 0.9*n)) 
    
    P <- generatePDCBM(theta, p11, p12, p22, prop=0.1)
    rho <- rho.fun(P, Cstar)
    
    
    df[cnt:(cnt+3), 3] <- n
    df[cnt:(cnt+3), 6] <- rho
    
    # Borgatti and Everett
    start = proc.time()[3]
    C <- borgattiCpp(A)
    end = proc.time()[3]
    df[cnt, 4] <- class_acc(C, Cstar)
    df[cnt, 5] <- end - start
    
    # Rombach
    start = proc.time()[3]
    C <- C_rombach(A, k)
    end = proc.time()[3]
    df[cnt+1, 4] <- class_acc(C, Cstar)
    df[cnt+1, 5] <- end - start 
    
    # Curcuringu
    start = proc.time()[3]
    C <- C_cur(A)
    end = proc.time()[3]
    df[cnt+2, 4] <- class_acc(C, Cstar)
    df[cnt+2, 5] <- end - start  
    
    # BayesSBM
    start = proc.time()[3]
    C <- C_SBM(A, 100, 50)
    end = proc.time()[3]
    df[cnt+3, 4] <- class_acc(C, Cstar)
    df[cnt+3, 5] <- end - start  
    
    cnt = cnt + 4
    #save(df, file = "sims_dcbm_class_n_041824.RData")
  }
  
  print(n)
}


load("sims_dcbm_class_n_041824.RData")

df_plot <- df %>% group_by(n, method) %>% summarise(accuracy = mean(class), time = mean(time), rho=mean(rho))

df_plot$method[df_plot$method=="BE"] <- "Proposed"

df_plot$method <- factor(df_plot$method, levels = c("Proposed", "BayesSBM", "Curcuringu", "Rombach"))


p1 <- ggplot(df_plot, aes(x=n, y=accuracy, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylim(0, 1)+
  ylab("Detection Accuracy")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

p2 <- ggplot(df_plot, aes(x=n, y=time, color=method, linetype=method))+
  geom_line(linewidth=1.2)+
  ylab("Time (sec)")+
  scale_colour_manual(name = "Method",values = myColors[c(1,2,5,6)])+
  labs(color="Method", linetype="Method")+
  theme(text = element_text(size = 16))+
  theme_bw()

ggarrange(p1, p2, nrow = 1, legend = "bottom", common.legend = TRUE)

#ggsave(file = "sims_dcbm_class_n_041824.pdf", device="pdf", width=6, height=4, units="in")











