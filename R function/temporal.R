resul=cbind(1:dim(resul)[1],resul)



plot(resul[c(1:1268)[str_detect(row.names(resul),'St')],1])

plot(resul[str_detect(row.names(resul),'dt'),2][1:cantidad])
plot(resul[str_detect(row.names(resul),'Ct'),2][1:cantidad])



plot(resul[907:1023,2])
plot(resul[1019:1135,2])

resul[1268,] 
resul[547:600,]
plot(resul[str_detect(row.names(resul),'Rt'),2][1:cantidad])
plot(resul[str_detect(row.names(resul),'IFR'),2][1:cantidad])


plot(resul[str_detect(row.names(resul),'St'),2])
plot(resul[str_detect(row.names(resul),'St'),2])
plot(resul[str_detect(row.names(resul),'St'),2])


csvfiles=dir(system.file('misc',package = 'rstan'),
             pattern = 'C:/Users/Mercedes Mena/Desktop/Resultados/distribucion normal primer modelo solo117 semanas prub diagnostic_1-7.csv',
             full.names = TRUE)
dir("../..", pattern = "^[a-lr]", full.names = TRUE, ignore.case = TRUE)

fit=read_stan_csv('C:/Users/Mercedes Mena/Desktop/Resultados/distribucion normal primer modelo solo117 semanas prub diagnostic_1-7.csv')


csvfiles <- dir('C:/Users/Mercedes Mena/Desktop/Resultados/',
                pattern = 'distribucion normal primer modelo solo117 semanas prub2_[2-4].csv', full.names = TRUE)
fit <- read_stan_csv(csvfiles)



