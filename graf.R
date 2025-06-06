
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(tibble)


### se comienza por las variantes

variantes=data_stan$Variantes
variantes$Seman=1:117

variantes1=as.data.frame(1:117)
colnames(variantes1)='Semana'
variantes1$Frecuencia=variantes[,1]

for (i in 2:15) {
  temporal=variantes[,c(16,i)]
  colnames(temporal)=colnames(variantes1)
  variantes1=rbind(variantes1,temporal)
  
}

variantes1$Variante=rep(colnames(variantes)[-16],rep(117,15))
# variantes1$Semana=as.character(variantes1$Semana)

round(variantes1$Frecuencia*1000)

salida=as.data.frame(rep(variantes1$Semana,round(variantes1$Frecuencia*1000)))
colnames(salida)='Semana'
salida$Variante=rep(variantes1$Variante,round(variantes1$Frecuencia*1000))


p_varantes <- salida %>%
  ggplot( aes(x=Variante,y=Semana,fill=Variante, color=Variante)) +
  geom_violin(width=1.5, size=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none"
  ) +
  coord_flip() + # This switch X and Y axis and allows to get the horizontal version
  xlab("Variantes") +
  ylab("Semana")

p_varantes

resul=as.data.frame(aaa$summary)
posterior=as.array(fit)
color_scheme_set("red")
mcmc_areas(posterior,
           pars = c('lp__'),
           prob = 0.95)+
  scale_y_discrete(
    labels = 'lp_',
  )


# parametros de las variantes del rt negativos y wan
mcmc_areas(posterior,
           pars = c('alpa_j_M[4]','alpa_j_M[5]','alpa_j_M[6]','alpa_j_M[8]','alpa_j_M[12]','alpa_j_M[14]','alpa_j_M[15]'),
           prob = 0.9)+
  scale_y_discrete(
    labels = colnames(variantes)[c(4,5,6,8,12,14,15)]
  )

# parametros de las variantes del rt Positivos y wan
mcmc_areas(posterior,
           pars = c('alpa_j_M[1]','alpa_j_M[2]','alpa_j_M[3]','alpa_j_M[7]',
                    'alpa_j_M[9]','alpa_j_M[10]','alpa_j_M[11]','alpa_j_M[13]','alpa_j_M[14]'),
           prob = 0.9)+
  scale_y_discrete(
    labels = colnames(variantes)[c(1,2,3,7,9,10,11,13,14)]
  )


# parametros de los momentos RT
mcmc_areas(posterior,
           pars = c('momento_Rt[1]','momento_Rt[2]','momento_Rt[3]','momento_Rt[4]'),
           prob = 0.9)+
  scale_y_discrete(
    labels = c('Momento 1','Momento 2','Momento 3','Momento 4')
  )


# Distribucion de R0
mcmc_areas(posterior,
           pars = c('R0'),
           prob = 0.9)+
  scale_y_discrete(
    labels = c('R0')
  )
resul$CV= resul$sd/resul$mean*100
resul$CV[2]





# parametros de las variantes del IFR 
color_scheme_set("red")
mcmc_areas(posterior,
           pars = c('alpa_IFR_j_M[1]','alpa_IFR_j_M[2]','alpa_IFR_j_M[3]','alpa_IFR_j_M[4]','alpa_IFR_j_M[5]','alpa_IFR_j_M[6]','alpa_IFR_j_M[7]','alpa_IFR_j_M[8]',
                    'alpa_IFR_j_M[9]','alpa_IFR_j_M[10]','alpa_IFR_j_M[11]','alpa_IFR_j_M[12]','alpa_IFR_j_M[13]','alpa_IFR_j_M[14]','alpa_IFR_j_M[15]'),
           prob = 0.9)+
  scale_y_discrete(
    labels = colnames(variantes)
  )

# parametros de los momentos IFR
color_scheme_set("red")
mcmc_areas(posterior,
           pars = c('momento_IFRt[1]','momento_IFRt[2]','momento_IFRt[3]','momento_IFRt[4]'),
           prob = 0.9)+
  scale_y_discrete(
    labels = c('Momento 1','Momento 2','Momento 3','Momento 4')
  )

# parametros otros dos IFR
color_scheme_set("red")
mcmc_areas(posterior,
           pars = c('alpa_IFR_c'),
           prob = 0.9)+
  scale_y_discrete(
    labels = c('Activos')
  )

# parametros otros dos IFR
color_scheme_set("red")
mcmc_areas(posterior,
           pars = c('alpa_IFR_s'),
           prob = 0.9)+
  scale_y_discrete(
    labels = c('Vacunados')
  )

#IFR0
color_scheme_set("red")
mcmc_areas(posterior,
           pars = c('IFR0'),
           prob = 0.9)+
  scale_y_discrete(
    labels = c('IFR0')
  )

# st

color_scheme_set("green")
mcmc_areas(posterior,
           pars = c('alpa_error_Ct','alpa_error'),
           prob = 0.9)+
  scale_y_discrete(
    labels = c('Infección','Vacunación')
  )

# st

color_scheme_set("green")
mcmc_areas(posterior,
           pars = c('vac_susep'),
           prob = 0.9)+
  scale_y_discrete(
    labels = c('Vacunados Susceptibles')
  )


library(tibble)

#susceptibles

df=as.data.frame(aaa$summary[c(1:1268)[str_detect(row.names(resul),'St')],])
hill_bayes_ribbons <- tibble(
  Semana = 1:117,
  Susceptibles = df$mean,
  pred_lower = df$`25%`,
  pred_upper = df$`75%`,
  fitted_lower = df$`2.5%`,
  fitted_upper = df$`97.5%`)

colourcodes <- as.character(color_scheme_get("green")) 
i=3
susep=ggplot(hill_bayes_ribbons , aes(x=Semana, y=Susceptibles)) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.5, fill=colourcodes[i]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.2, fill=colourcodes[i]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Susceptibles), colour=colourcodes[i],
            size=1)+
  theme_classic()+
  geom_vline(xintercept = 16,alpha=0.4)+
  geom_vline(xintercept = 37,alpha=0.4)+
  geom_vline(xintercept = 93,alpha=0.4)+
  geom_text(x=6,y=0.7,label=c('M 1'),alpha=1,color='red',size=5)+
  geom_text(x=25,y=0.7,label=c('M 2'),alpha=1,color='red',size=5)+
  geom_text(x=60,y=0.7,label=c('M 3'),alpha=1,color='red',size=5)+
  geom_text(x=110,y=0.7,label=c('M 4'),alpha=1,color='red',size=5)+
  transition_reveal(Semana)


# Save at gif:
anim_save("D:/Tesis de master/Programa/graficos/anima/susceptibles.gif", width = 10, height = 4)




df=as.data.frame(aaa$summary[902:1018,])
# df=as.data.frame(aaa$summary[902:942,])

# df=as.data.frame(aaa$summary[902:918,])
hill_bayes_ribbons <- tibble(
  Semana = 1:117,
  Rt = df$mean,
  pred_lower = df$`25%`,
  pred_upper = df$`75%`,
  fitted_lower = df$`2.5%`,
  fitted_upper = df$`97.5%`)

colourcodes <- as.character(color_scheme_get("blue")) 
i=4
ggplot(hill_bayes_ribbons , aes(x=Semana, y=Rt)) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.5, fill=colourcodes[i]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.2, fill=colourcodes[i]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Rt), colour=colourcodes[i],
            size=1)+
  theme_classic()+
  geom_vline(xintercept = 16,alpha=0.4)+
  geom_vline(xintercept = 37,alpha=0.4)+
  geom_vline(xintercept = 93,alpha=0.4)+
  geom_text(x=6,y=0.7,label=c('M 1'),alpha=1,color='red',size=5)+
  geom_text(x=25,y=0.7,label=c('M 2'),alpha=1,color='red',size=5)+
  geom_text(x=60,y=0.7,label=c('M 3'),alpha=1,color='red',size=5)+
  geom_text(x=110,y=0.7,label=c('M 4'),alpha=1,color='red',size=5)+
transition_reveal(Semana)


# Save at gif:
anim_save("D:/Tesis de master/Programa/graficos/anima/Rt.gif", width = 10, height = 4)




df=as.data.frame(aaa$summary[1019:1135,])

# df=as.data.frame(aaa$summary[902:918,])
hill_bayes_ribbons <- tibble(
  Semana = 1:117,
  Rt = df$mean,
  pred_lower = df$`25%`,
  pred_upper = df$`75%`,
  fitted_lower = df$`2.5%`,
  fitted_upper = df$`97.5%`)

colourcodes <- as.character(color_scheme_get("red")) 
i=4
ggplot(hill_bayes_ribbons , aes(x=Semana, y=Rt)) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.5, fill=colourcodes[i]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.2, fill=colourcodes[i]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Rt), colour=colourcodes[i],
            size=1)+
  theme_classic()+
  geom_vline(xintercept = 16,alpha=0.4)+
  geom_vline(xintercept = 37,alpha=0.4)+
  geom_vline(xintercept = 93,alpha=0.4)+
  geom_text(x=6,y=0.004,label=c('M 1'),alpha=1,color='red',size=5)+
  geom_text(x=25,y=0.004,label=c('M 2'),alpha=1,color='red',size=5)+
  geom_text(x=60,y=0.004,label=c('M 3'),alpha=1,color='red',size=5)+
  geom_text(x=110,y=0.004,label=c('M 4'),alpha=1,color='red',size=5)+
  ylab('IFRt')+
  transition_reveal(Semana)


# Save at gif:
anim_save("D:/Tesis de master/Programa/graficos/anima/IFRt.gif", width = 10, height = 4)







df=as.data.frame(aaa$summary[317:433,])


# df=as.data.frame(aaa$summary[317:332,])
hill_bayes_ribbons <- tibble(
  Semana = 1:117,
  Rt = df$mean,
  pred_lower = df$`25%`,
  pred_upper = df$`75%`,
  fitted_lower = df$`2.5%`,
  fitted_upper = df$`97.5%`,
  Decesos=data_Cuba_Week$new_deaths)

colourcodes <- as.character(color_scheme_get("gray")) 
i=4
ggplot(hill_bayes_ribbons , aes(x=Semana, y=Rt)) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.9, fill=colourcodes[i]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.7, fill=colourcodes[i]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Rt), colour=colourcodes[i],
            size=1)+
  geom_bar(aes(x=Semana, y=Decesos),stat = 'identity',alpha=0.2)+
  theme_classic()+
  geom_vline(xintercept = 16,alpha=0.4)+
  geom_vline(xintercept = 37,alpha=0.4)+
  geom_vline(xintercept = 93,alpha=0.4)+
  geom_text(x=6,y=200,label=c('M 1'),alpha=1,color='red',size=5)+
  geom_text(x=25,y=200,label=c('M 2'),alpha=1,color='red',size=5)+
  geom_text(x=60,y=200,label=c('M 3'),alpha=1,color='red',size=5)+
  geom_text(x=110,y=200,label=c('M 4'),alpha=1,color='red',size=5)+
  ylab('Muertes por Semana')+
  transition_reveal(Semana)


# Save at gif:
anim_save("D:/Tesis de master/Programa/graficos/anima/nuevas muertes.gif", width = 10, height = 4)




a=200
df=as.data.frame(aaa$summary[a:(a+116),])

hill_bayes_ribbons <- tibble(
  Semana = 1:117,
  Rt = df$mean,
  pred_lower = df$`25%`,
  pred_upper = df$`75%`,
  fitted_lower = df$`2.5%`,
  fitted_upper = df$`97.5%`,
  Casos=data_Cuba_Week$new_cases)

colourcodes <- as.character(color_scheme_get("purple")) 
i=4
ggplot(hill_bayes_ribbons , aes(x=Semana, y=Rt)) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.9, fill=colourcodes[i]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.7, fill=colourcodes[i]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Rt), colour=colourcodes[i],
            size=1)+
  geom_bar(aes(x=Semana, y=Casos),stat = 'identity',alpha=0.2,fill=colourcodes[i])+
  theme_classic()+
  geom_vline(xintercept = 16,alpha=0.4)+
  geom_vline(xintercept = 37,alpha=0.4)+
  geom_vline(xintercept = 93,alpha=0.4)+
  geom_text(x=6,y=40000,label=c('M 1'),alpha=1,color='red',size=5)+
  geom_text(x=25,y=40000,label=c('M 2'),alpha=1,color='red',size=5)+
  geom_text(x=60,y=40000,label=c('M 3'),alpha=1,color='red',size=5)+
  geom_text(x=110,y=40000,label=c('M 4'),alpha=1,color='red',size=5)+
  ylab('Nuevos Casos Por Semana')+
transition_reveal(Semana)


# Save at gif:
anim_save("D:/Tesis de master/Programa/graficos/anima/nuevas casos.gif", width = 10, height = 4)




#para la evaluacion


df=as.data.frame(aaa$summary[317:433,])


# df=as.data.frame(aaa$summary[317:332,])
temporal <- as.data.frame(tibble(
  Semana = 1:117,
  Rt = df$mean,
  pred_lower = df$`25%`,
  pred_upper = df$`75%`,
  fitted_lower = df$`2.5%`,
  fitted_upper = df$`97.5%`,
  Decesos=data_Cuba_Week$new_deaths))


#### r cuadrado
Rcuad=1-sum((temporal$Decesos-temporal$Rt)^2)/sum((temporal$Decesos-mean(temporal$Decesos))^2)

shapiro.test(temporal$Decesos-temporal$Rt)

hist(temporal$Decesos-temporal$Rt)
Residuos=temporal$Decesos-temporal$Rt

shapiro.test(rnorm(20))
ggplot(data = temporal, aes(x = seq_along(Residuos), 
                             y = temporal$Decesos-temporal$Rt)) +
  geom_point(aes(color = Residuos)) +
  scale_color_gradient2(low = "blue3", mid = "grey", high = "red") +
  geom_line(size = 0.3) + 
  labs(title = "Distribución de los residuos", x = "index", y = "residuo")+ 
  geom_hline(yintercept = 0) + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
theme_classic()+xlab('')+ylab('')


require(nortest)  # Se debe haber instalado nortest
ad.test(temporal$Decesos-temporal$Rt)

resul=as.data.frame(aaa$summary)

rhats <- resul$Rhat
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)

temporal=aaa$summary[rhats>1.1,]

### 11 parameters con r hat mayor que 1.1


# tama;o de muestreo efectivo
ratios_cp <- neff_ratio(fit)
print(ratios_cp)


mcmc_neff(ratios_cp, size = 0.1)
temporal=aaa$summary[ratios_cp>1.1,]



posterior_cp <- as.array(fit)

library("bayesplot")
library("ggplot2")
library("rstan")      

np_cp <- nuts_params(fit)
color_scheme_set("darkgray")
mcmc_parcoord(posterior_cp, np = np_cp)






# variables que me quedaban





df=as.data.frame(aaa$summary[317:433,])


# df=as.data.frame(aaa$summary[317:332,])
hill_bayes_ribbons <- tibble(
  Semana = 1:117,
  Rt = df$mean,
  pred_lower = df$`25%`,
  pred_upper = df$`75%`,
  fitted_lower = df$`2.5%`,
  fitted_upper = df$`97.5%`,
  Decesos=data_Cuba_Week$new_deaths,
  Casos=data_Cuba_Week$new_cases,
  vacunados=data_Cuba_Week$total_vaccinations)
colourcodes1 <- as.character(color_scheme_get("purple")) 
i=4
colourcodes2 <- as.character(color_scheme_get("green")) 
i=4
colourcodes <- as.character(color_scheme_get("gray")) 
i=4
ggplot(hill_bayes_ribbons , aes(x=Semana, y=Rt)) +
  # geom_bar(aes(x=Semana, y=Casos),stat = 'identity',alpha=0.2,fill=colourcodes1[i])+
  # geom_bar(aes(x=Semana, y=Decesos),stat = 'identity',alpha=0.2,fill=colourcodes[i])+
  geom_bar(aes(x=Semana, y=vacunados),stat = 'identity',alpha=0.8,fill=colourcodes2[3])+
  theme_classic()+
  geom_vline(xintercept = 16,alpha=0.4)+
  geom_vline(xintercept = 37,alpha=0.4)+
  geom_vline(xintercept = 93,alpha=0.4)+
  geom_text(x=6,y=4000000,label=c('M 1'),alpha=1,color='red',size=5)+
  geom_text(x=25,y=4000000,label=c('M 2'),alpha=1,color='red',size=5)+
  geom_text(x=60,y=4000000,label=c('M 3'),alpha=1,color='red',size=5)+
  geom_text(x=110,y=4000000,label=c('M 4'),alpha=1,color='red',size=5)+
  ylab('Total de vacunados por Semana')



a=200
df=as.data.frame(aaa$summary[a:(a+116),])

hill_bayes_ribbons <- tibble(
  Semana = 1:117,
  Rt = df$mean,
  pred_lower = df$`25%`,
  pred_upper = df$`75%`,
  fitted_lower = df$`2.5%`,
  fitted_upper = df$`97.5%`,
  Casos=data_Cuba_Week$new_cases)

colourcodes <- as.character(color_scheme_get("purple")) 
i=4
ggplot(hill_bayes_ribbons , aes(x=Semana, y=Rt)) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.9, fill=colourcodes[i]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.7, fill=colourcodes[i]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Rt), colour=colourcodes[i],
            size=1)+
  geom_bar(aes(x=Semana, y=Casos),stat = 'identity',alpha=0.2,fill=colourcodes[i])+
  theme_classic()+
  geom_vline(xintercept = 16,alpha=0.4)+
  geom_vline(xintercept = 37,alpha=0.4)+
  geom_vline(xintercept = 93,alpha=0.4)+
  geom_text(x=6,y=40000,label=c('M 1'),alpha=1,color='red',size=5)+
  geom_text(x=25,y=40000,label=c('M 2'),alpha=1,color='red',size=5)+
  geom_text(x=60,y=40000,label=c('M 3'),alpha=1,color='red',size=5)+
  geom_text(x=110,y=40000,label=c('M 4'),alpha=1,color='red',size=5)+
  ylab('Nuevos Casos Por Semana')




np_cp <- nuts_params(fit)

lp_ncp <- log_posterior(fit)
np_ncp <- nuts_params(fit)
color_scheme_set("mix-brightblue-gray")
posterior=as.array(fit)
color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior, pars = "R0", np = np_cp) +
  xlab("Iteraciones despés del calentamiento")

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior, pars = "IFR0", np = np_cp) +
  xlab("Iteraciones despés del calentamiento")


color_scheme_set("darkgray")
mcmc_parcoord(posterior, np = np_cp)

compare_cp_ncp <- function(cp_plot, ncol = 1, ...) {
  bayesplot_grid(
    cp_plot, 
    grid_args = list(ncol = ncol),
    subtitles = c("Centered parameterization"),
    ...
  )
}

compare_cp_ncp(
  mcmc_acf(posterior, pars = "R0", lags = 10)
)

mcmc_acf(posterior, pars = c('alpa_j_M[1]','alpa_j_M[2]'), lags = 10)
mcmc_acf(posterior, pars = c('IFR0','R0'), lags = 10)
