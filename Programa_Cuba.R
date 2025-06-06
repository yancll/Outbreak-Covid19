#Librerias######################################################################
library(readr)
library(timeDate)
library(tidyverse)
library(foreign)
library(openxlsx)
library(readxl)
library(tidyverse)
library(srvyr)
library(timeDate)
library(tidyverse)
library(foreign)
library(openxlsx)
library(readxl)
library(tidyverse)
library(srvyr)
library(sp)
library(maptools)
library(RColorBrewer)
library(KernSmooth)
library(classInt)
library(tidyverse)
library(reshape2)
library(srvyr)
library(viridis)
library(rainbow)
library(plotly)
library(ggsankey)
library(timeDate)
library(tidyverse)
library(foreign)
library(openxlsx)
library(readxl)
library(tidyverse)
library(srvyr)
library(sp)
library(maptools)
library(RColorBrewer)
library(KernSmooth)
library(classInt)
library(tidyverse)
library(reshape2)
library(srvyr)
library(viridis)
#library(rainbow)
library(plotly)
library(haven)
library(ggsci)
library(vcd)
library(rstan)
library(bayesplot)
library(tidyr)
library(dplyr)
library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(EnvStats)
library(matrixStats)
library(scales)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)
library(svglite)
library(ggplot2)
library("stringr")    
#Definir directorio de trabajo##################################################
setwd("F:/Tesis de master/Programa")
#Leer funciones auxiliares######################################################
source("R function/Leer_transformar_datos.R")
#Obrener los datos##############################################################
data_Cuba_Week=get_data_week()
data_Cuba_varian=get_variant(data_Cuba_Week )
data_Cuba_varian$omicron[117]=1


# 
# vector [1500]IS; // SI in reverse order
# vector [1500]IM; // SI in reverse order
# vector [1500]IR; // SI in reverse order
# vector [1500]ICV; // SI in reverse order
# vector [4] lar;
# for (i in 1:1500){
#   IS[i]=gamma_cdf(L_sup_dia[i], 5.3*(5.3/0.26)/0.26,(5.3/0.26)/0.26)-gamma_cdf(L_inf_dia[i], 5.3*(5.3/0.26)/0.26,(5.3/0.26)/0.26);
#   IM[i]=gamma_cdf(L_sup_dia[i], 21.33660,1/1.03174)-gamma_cdf(L_inf_dia[i], 21.33660,1/1.03174);
#   IR[i]=gamma_cdf(L_sup_dia[i], 18.62067,1/1.18892)-gamma_cdf(L_inf_dia[i], 18.62067,1/1.18892);
#   ICV[i]=gamma_cdf(L_sup_dia[i], 2.061856,0.003436426)-gamma_cdf(L_inf_dia[i], 2.061856,0.003436426);
#   if (IS[i]==0){ lar[1]+=1;}
#   if (IM[i]==0){ lar[2]+=1;}
#   if (IR[i]==0){ lar[3]+=1;}
#   if (ICV[i]==0){ lar[4]+=1;}
# }
# }

L_inf_dia=c(0,1:1499+0.5)
L_sup_dia=c(1:1500+0.5)
l_sup_7=(1:1500)*7

#IM=pgamma(L_sup_dia, shape=21.33660,rate =1/1.03174)-pgamma(L_inf_dia, shape=18.33660,rate =1/1.03174)


IS=pgamma(L_sup_dia, shape=5.3*(5.3/0.26)/0.26,rate=(5.3/0.26)/0.26)-pgamma(L_inf_dia, shape=5.3*(5.3/0.26)/0.26,rate = (5.3/0.26)/0.26)
IM=pgamma(L_sup_dia, shape=21.33660,rate =1/1.03174)-pgamma(L_inf_dia, shape=21.33660,rate =1/1.03174)
IR=pgamma(l_sup_7, shape=18.62067,rate =1/1.18892)
IC=pgamma(l_sup_7,shape= 2.061856,rate =0.003436426)
IV=pgamma(l_sup_7,shape= 2.061856,rate =0.003436426)

# IV=pgamma(l_sup_7,shape= 2.382732276,rate =0.003694158)
IR=IR[1500:1]
IC=IC[1500:1]
IV=IV[1500:1]
lar=c(sum(IS==0),sum(IM==0),sum(IR==0),sum(IV==0))

data_Cuba_varian[83:90,7]=1
CFR=data_Cuba_Week$total_deaths/data_Cuba_Week$total_cases
CFR[1:6]=CFR[7]
# CFR=(min(CFR)-CFR)/(min(CFR)-max(CFR))
# CFR=1/(1+exp(-CFR))

momento1=rep(0,length(CFR))
momento1[1:16]=1

momento2=rep(0,length(CFR))
momento2[17:37]=1


momento3=rep(0,length(CFR))
momento3[38:93]=1
momento4=rep(0,length(CFR))
momento4[94:length(CFR)]=1

cantidad=117
momento1=momento1[1:cantidad]
momento2=momento2[1:cantidad]
momento3=momento3[1:cantidad]
momento4=momento4[1:cantidad]
data_Cuba_varian=data_Cuba_varian[1:cantidad,]
data_Cuba_Week=data_Cuba_Week[1:cantidad,]
CFR=CFR[1:cantidad]
CFR=rep(0,length(CFR))
# nuevos_v=c(rep(0,7),data_Cuba_Week$new_vaccinations)
# data_Cuba_Week$new_vaccinations=nuevos_v[1:cantidad]
# data_Cuba_Week$total_vaccinations=cumsum(nuevos_v[1:cantidad])

data_stan=list(n=dim(data_Cuba_Week)[1],
               m=dim(data_Cuba_Week)[2],
               N=11000000,
               Casos_datos=data_Cuba_Week$new_cases,
               Muertes_datos=data_Cuba_Week$new_deaths,
               Nuevos_V=data_Cuba_Week$new_vaccinations,
               Total_V=data_Cuba_Week$total_vaccinations,
               Variantes=data_Cuba_varian[,-1],
               nv=dim(data_Cuba_varian)[2]-1,
               M=diag(dim(data_Cuba_Week)[1]),
               IS=IS,
               IM=IM,
               IR=IR,
               IC=IC,
               IV=IV,
               lar=lar,
               CFR=CFR,
               momento1=momento1,
               momento2=momento2,
               momento3=momento3,
               momento4=momento4,
               Casos_d_acum=data_Cuba_Week$total_cases)
options(mc.cores = parallel::detectCores())
set_cppo(mode='fast')
rstan_options(auto_write = TRUE)


model_stan <- stan_model('Stan function/Stan_Modelo_Cuba_8.stan', verbose = TRUE)
inicial=get_posterior_mean(fit)
ss=100
fit <- sampling(model_stan, data =data_stan,chains=4,warmup=ss*1000,iter=ss*2000,control = list(max_treedepth = 20,adapt_delta = 0.95),
                sample_file=paste0('C:/Users/Yan/Desktop/Trabajo/Resultados/','distribucion normal primer modelo solo', cantidad,' semanas prub2','.csv'),
                diagnostic_file=paste0('C:/Users/Yan/Desktop/Trabajo/Resultados/','distribucion normal primer modelo solo', cantidad,' semanas prub diagnostic2','.csv'))


save(fit,data_stan,
     file=paste0('C:/Users/Yan/Desktop/Trabajo/Resultados/','distribucion normal primer modelo solon', cantidad,' semanas prub','.Rdata'))

resul=as.data.frame(monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE)))

resul=cbind(1:dim(resul)[1],resul)
rrr=resul[str_detect(row.names(resul),'dt'),2][1:cantidad]
sum((rrr-mean(data_stan$Muertes_datos))^2)/sum((data_stan$Muertes_datos-mean(data_stan$Muertes_datos))^2)  
1-sum((data_stan$Muertes_datos - rrr)^2)/sum((data_stan$Muertes_datos-mean(data_stan$Muertes_datos))^2)


resul=as.data.frame(resul)
save(fit,data_stan,resul,
     file=paste0('C:/Users/Yan/Desktop/Trabajo/Resultados/','distribucion normal primer modelo solo', cantidad,' semanas prubn1','.Rdata'))
aaa=summary(fit)
save(resul,aaa,
     file=paste0('C:/Users/Yan/Desktop/Trabajo\Resultados/','distribucion normal primer modelo solo', cantidad,' resul','.Rdata'))

