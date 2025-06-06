

get_data_day=function(){
  #los datos estan tomados de The Oxford Martin Programme on Global Development
  #https://ourworldindata.org/covid-vaccinations
  #https://github.com/owid/covid-19-data/tree/master/public/data 
  #la data general de todos los paises son las mismas bd
  Global_data <- read_csv("Data/owid-covid-data (2).csv")
  cuba_data=Global_data[Global_data$location=='Cuba',]
  n=dim(cuba_data)[1]
  semanas=floor(((n-3)/7))
  cuba_data$Semanas=rep(c(1:(semanas+2)),c(3,rep(7,semanas),n-3-semanas*7))
  cuba_data=cuba_data[-816,]
  return(cuba_data)
}


get_data_week=function(){
  Global_data <- as.data.frame(read_csv("Data/owid-covid-data (2).csv"))
  cuba_data=Global_data[Global_data$location=='Cuba',]
  n=dim(cuba_data)[1]
  semanas=floor(((n-3)/7))
  cuba_data$Semanas=rep(c(1:(semanas+2)),c(3,rep(7,semanas),n-3-semanas*7))
  cuba_data=cuba_data[-816,]
  for (i in 5:dim(cuba_data)[2]) {
    if(is.numeric(cuba_data[,i])){
      cuba_data[is.na(cuba_data[,i]),i]=0
    }
  }
  
  cont1=0
  cont2=0
  for (i in 1:dim(cuba_data)[1]) {
    if(cuba_data$people_fully_vaccinated[i]!=0){
      cont1=1
      cont2=cuba_data$people_fully_vaccinated[i]
    }
    if(cont1==1){
      if(cuba_data$people_fully_vaccinated[i]==0){
        cuba_data$people_fully_vaccinated[i]=cont2
      }
    
    }
    
  }
  
  cuba_data1= cuba_data%>% group_by(Semanas)
  cuba_data_week=cuba_data1%>%  summarise(
    iso_code=unique(iso_code),
    continent=unique(continent),  
    location=unique(location),
    date0=date[1],
    date1=date[length(date) ],
    total_day=n(),
    new_cases=sum(new_cases),
    total_cases=total_cases[length(total_cases)],
    new_deaths=sum(new_deaths),
    total_deaths=total_deaths[length(total_deaths)],
    total_vaccinations=people_fully_vaccinated[length(people_fully_vaccinated)]
    
  )
  cuba_data_week$new_vaccinations=0
  for (i in 2:dim(cuba_data_week)[1]) {
    cuba_data_week$new_vaccinations[i]=cuba_data_week$total_vaccinations[i]-cuba_data_week$total_vaccinations[i-1]
  }
  
  return(cuba_data_week)
}


get_variant=function(data_Cuba_Week){
  variantes <- read_excel("Data/variantes.xlsx", sheet = "Hoja1",
                                       col_types = c( "date","text", "text", "text"))
  colnames(variantes)=tolower(colnames(variantes))
  for (i in 1:dim(variantes)[1]) {
    variantes[i,3]=tolower(variantes[i,3])
  }
  variantes=as.data.frame(variantes)
  variantes$ajustada=variantes$variantes
  variantes$ajustada[variantes$ajustada%in%c('b.1.1.519','gamma','lambda','pm4','pm8','b.1.1.519')]='otra'
  salida=matrix(0,dim(data_Cuba_Week)[1],1+length(levels(as.factor(variantes$ajustada))))
  variantes_nombres=levels(as.factor(variantes$ajustada))
  colnames(salida)=c('Semana',variantes_nombres)
  salida[,1]=data_Cuba_Week$Semanas
  for (i in 1:dim(salida)[1]) {
    for (j in 2:dim(salida)[2]) {
      salida[i,j]=sum((variantes$date>=data_Cuba_Week$date0[i])*(variantes$date<=data_Cuba_Week$date1[i])*(variantes$ajustada%in%variantes_nombres[j-1]))
    }
    
  }
  salida[1:42,colnames(salida)%in%"wuhan"]=1
  salida[44,colnames(salida)%in%"wuhan"]=0
  sumas=colSums(t(salida[,2:dim(salida)[2]]))
  sumas[sumas==0]=1
  for (i in 2:dim(salida)[2]) {
    salida[,i]=salida[,i]/sumas
    
  }
  salida=as.data.frame(salida)
  return(salida)
}






