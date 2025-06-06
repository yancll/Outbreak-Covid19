data {
  int<lower=0> n;
  int<lower=0> m;
  int<lower=0> N;
  int<lower=0> nv;
  vector [n]Casos_datos;
  int Casos_d_acum[n];
  int Muertes_datos[n];
  vector [n]Nuevos_V;
  vector [1500]IS; // SI in reverse order
  vector [1500]IM; // SI in reverse order
  vector [1500]IR; // SI in reverse order
  vector [1500]IV; // SI in reverse order
  vector [1500]IC;
  vector [4] lar;
  vector [n] Total_V;
  matrix [n,nv] Variantes;
  vector [n] CFR;
  vector [n] momento1;
  vector [n] momento2;
  vector [n] momento3;
  vector [n] momento4;
}


parameters {
  real<lower=0> tau;
  real<lower=0> R0;
  real<lower=0.008,upper=0.016> IFR0;
  real<lower=0> kapa;
  real<lower=0> kapa_IFR0;
  real<lower=0,upper=1> alpa_error;
  real<lower=0> jota_alpa_error;
  // real salva_Rt;
  real<lower=0.5,upper=1> alpa_error_Ct;
  real<lower=0> jota_alpa_error_Ct;
  real<lower=0> fi;
  real<lower=0>j_salva;
  real<lower=0> c_inicial;
  // real alpa_s;
  // real alpa_c;
  // real<lower=0> jota_s;
  // real<lower=0> jota_c;
  real alpa_IFR_s;
  real alpa_IFR_c;
  real alpa_IFR_CFR;
  real <lower=0>jota_IFR_CFR;
  real<lower=0> jota_IFR_s;
  real<lower=0> jota_IFR_c;
  vector [nv] alpa_j_M;
  vector  [nv]alpa_IFR_j_M;
  real<lower=0.5,upper=1> vac_susep;
  real<lower=0,upper=1> j_vac_susep;
  
  vector  [4] momento_Rt;
  vector<lower=0>  [4]j_Rt;
  vector  [4] momento_IFRt;
  vector<lower=0>  [4]j_IFRt;
  
  vector[5] weekly_effect_R0;
  real<lower=0, upper=1> weekly_rho_R0;
  real<lower=0, upper=1> weekly_rho1_R0;
  real<lower=0> weekly_sd_R0;
  
  vector[5] weekly_effect_IFR;
  real<lower=0, upper=1> weekly_rho_IFR;
  real<lower=0, upper=1> weekly_rho1_IFR;
  real<lower=0> weekly_sd_IFR;
  
}



transformed parameters {
  vector<lower=0,upper=1> [n] St;
  vector<lower=0>  [n]Ct;
  vector<lower=0>  [n]dt;
  vector<lower=0,upper=3000000>  [n]Ctacum;
  vector<lower=0,upper=11000>  [n]dtacum;
  vector  [n]Ct_R_ac;
  // vector  [n]Ct_R_ac_VC;
  vector  [n]Vt_R_ac;
  vector<lower=0,upper=4> [n]Rt;
  vector<lower=0,upper=0.02> [n]IFR;
  vector<lower=0> [n]Rec;
  vector<lower=0> [7]yy1;
  vector<lower=0> [7]yy2;
  real <lower=0> momento_dich;
  real<lower=0> borde1; 
  real<lower=0> borde2; 
  real<lower=0> borde3; 
  real<lower=0> borde4; 
   
  St[1]=1;
  Ct[1]=c_inicial;
  dt[1]=0.0001;
  Ctacum[1]=Ct[1];
  dtacum[1]=0;
  Ct_R_ac[1]=0;
  Vt_R_ac[1]=0;
  // Ct_R_ac_VC[1]=0;
  // Ct_R_ac_VC[1]=0;
  Rec[1]=0;
  momento_dich=1;
  
  Rt[1]=R0*2*St[1]*inv_logit(+dot_product(-alpa_j_M,Variantes[1,])-momento_Rt[1]*momento1[1]
  -momento_Rt[2]*momento2[1]-momento_Rt[3]*momento3[1]-momento_Rt[4]*momento4[1]-weekly_effect_R0[1]);
 
 // Rt[1]=R0*2*St[1]*inv_logit(+dot_product(-alpa_j_M,Variantes[1,])-salva_Rt*R0-momento_Rt[1]*momento1[1]
 //  -momento_Rt[2]*momento2[1]-momento_Rt[3]*momento3[1]-momento_Rt[4]*momento4[1]-weekly_effect_R0[1]);
 
  
  IFR[1]=IFR0*2*inv_logit(-weekly_effect_IFR[1]-alpa_IFR_CFR*CFR[1]-alpa_IFR_s*Total_V[1]-momento_IFRt[1]*momento1[1]-momento_IFRt[2]*momento2[1]-
  momento_IFRt[3]*momento3[1]-momento_IFRt[4]*momento4[1]-alpa_IFR_c*(Ct[1])+dot_product(-alpa_IFR_j_M,Variantes[1,]));
  
  
  
  for (i in 2:n){
    if(momento1[i]==1){momento_dich=1;}
    if(momento2[i]==1){momento_dich=2;}
    if(momento3[i]==1){momento_dich=3;}
    if(momento4[i]==1){momento_dich=4;}
      Vt_R_ac[i]=dot_product(vac_susep*Nuevos_V[1:(i-1)],IC[(1500-i+2):1500]);
      Ct_R_ac[i]=dot_product(Ct[1:(i-1)],IC[(1500-i+2):1500]);
      // Ct_R_ac_VC[i]=dot_product(Ct[1:(i-1)],IC[(1500-i+2):1500]);
      
      St[i]=(N+Ct_R_ac[i]+alpa_error*Vt_R_ac[i]-vac_susep*Total_V[i-1]-Ctacum[i-1])/N;
      
      if(momento1[i]==1){
        Rt[i]=R0*2*St[i]*inv_logit(dot_product(-alpa_j_M,Variantes[i,])-momento_Rt[1]*momento1[i]-momento_Rt[2]*momento2[i]-
      momento_Rt[3]*momento3[i]-momento_Rt[4]*momento4[i]-
      weekly_effect_R0[1]);
      }
      if(momento2[i]==1){
        Rt[i]=R0*2*St[i]*inv_logit(dot_product(-alpa_j_M,Variantes[i,])-momento_Rt[1]*momento1[i]-momento_Rt[2]*momento2[i]-
      momento_Rt[3]*momento3[i]-momento_Rt[4]*momento4[i]-
      weekly_effect_R0[2]);
      }
      if(momento3[i]==1){
        Rt[i]=R0*2*St[i]*inv_logit(dot_product(-alpa_j_M,Variantes[i,])-momento_Rt[1]*momento1[i]-momento_Rt[2]*momento2[i]-
      momento_Rt[3]*momento3[i]-momento_Rt[4]*momento4[i]-
      weekly_effect_R0[3]);
      }
      if(momento4[i]==1){
        Rt[i]=R0*2*St[i]*inv_logit(dot_product(-alpa_j_M,Variantes[i,])-momento_Rt[1]*momento1[i]-momento_Rt[2]*momento2[i]-
      momento_Rt[3]*momento3[i]-momento_Rt[4]*momento4[i]-
      weekly_effect_R0[4]);
      }
      
      
      // if(i<n){Rt[i]=R0*2*St[i]*inv_logit(dot_product(-alpa_j_M,Variantes[i,])-momento_Rt[1]*momento1[i]-momento_Rt[2]*momento2[i]-
      // momento_Rt[3]*momento3[i]-momento_Rt[4]*momento4[i]-salva_Rt*((Casos_datos[i]+Casos_datos[1+i])/(Casos_datos[i]+Casos_datos[i-1]))-
      // weekly_effect_R0[i]);}
      
      
      // if(i==n){Rt[i]=R0*2*St[i]*inv_logit(dot_product(-alpa_j_M,Variantes[i,])-momento_Rt[1]*momento1[i]-momento_Rt[2]*momento2[i]-
      // momento_Rt[3]*momento3[i]-momento_Rt[4]*momento4[i]-salva_Rt*((Casos_datos[i])/(Casos_datos[i-1]))-weekly_effect_R0[i]);}
      // 
      
      yy1=rep_vector(0,7);
      yy2=rep_vector(0,7);
      for (j in 1:(i-1)) {
        if((7*(j-1)+7)<(1500-lar[1]+7)){
          yy1[1]+=sum(IS[(7*(j-1)+1):(7*(j-1)+7)])*Ct[i-j]/7;
          yy1[2]+=sum(IS[(7*(j-1)+2):(7*(j-1)+8)])*Ct[i-j]/7;
          yy1[3]+=sum(IS[(7*(j-1)+3):(7*(j-1)+9)])*Ct[i-j]/7;
          yy1[4]+=sum(IS[(7*(j-1)+4):(7*(j-1)+10)])*Ct[i-j]/7;
          yy1[5]+=sum(IS[(7*(j-1)+5):(7*(j-1)+11)])*Ct[i-j]/7;
          yy1[6]+=sum(IS[(7*(j-1)+6):(7*(j-1)+12)])*Ct[i-j]/7;
          yy1[7]+=sum(IS[(7*(j-1)+7):(7*(j-1)+13)])*Ct[i-j]/7;
        }
        
        if((7*(j-1)+7)<(1500-lar[2]+7)){
          yy2[1]+=sum(IM[(7*(j-1)+1):(7*(j-1)+7)])*Ct[i-j]/7;
          yy2[2]+=sum(IM[(7*(j-1)+2):(7*(j-1)+8)])*Ct[i-j]/7;
          yy2[3]+=sum(IM[(7*(j-1)+3):(7*(j-1)+9)])*Ct[i-j]/7;
          yy2[4]+=sum(IM[(7*(j-1)+4):(7*(j-1)+10)])*Ct[i-j]/7;
          yy2[5]+=sum(IM[(7*(j-1)+5):(7*(j-1)+11)])*Ct[i-j]/7;
          yy2[6]+=sum(IM[(7*(j-1)+6):(7*(j-1)+12)])*Ct[i-j]/7;
          yy2[7]+=sum(IM[(7*(j-1)+7):(7*(j-1)+13)])*Ct[i-j]/7;
        }
      }
        yy1[1]=yy1[1]*Rt[i];
        yy1[2]=yy1[1]*IS[1]+yy1[2];
        yy1[2]=yy1[2]*Rt[i];
        yy1[3]=yy1[2]*IS[1]+yy1[1]*IS[2]+yy1[3];
        yy1[3]=yy1[3]*Rt[i];
        yy1[4]=yy1[3]*IS[1]+yy1[2]*IS[2]+yy1[1]*IS[3]+yy1[4];
        yy1[4]=yy1[4]*Rt[i];
        yy1[5]=yy1[4]*IS[1]+yy1[3]*IS[2]+yy1[2]*IS[3]+yy1[1]*IS[4]+yy1[5];
        yy1[5]=yy1[5]*Rt[i];
        yy1[6]=yy1[5]*IS[1]+yy1[4]*IS[2]+yy1[3]*IS[3]+yy1[2]*IS[4]+yy1[1]*IS[5]+yy1[6];
        yy1[6]=yy1[6]*Rt[i];
        yy1[7]=yy1[6]*IS[1]+yy1[5]*IS[2]+yy1[4]*IS[3]+yy1[3]*IS[4]+yy1[2]*IS[5]+yy1[1]*IS[6]+yy1[7];
        yy1[7]=yy1[7]*Rt[i];
        Ct[i]=sum(yy1);
        Ctacum[i]=Ct[i]+Ctacum[i-1];


        yy2[2]=yy1[1]*IM[1]+yy2[2];
        yy2[3]=yy1[2]*IM[1]+yy1[1]*IM[2]+yy2[3];
        yy2[4]=yy1[3]*IM[1]+yy1[2]*IM[2]+yy1[1]*IM[3]+yy2[4];
        yy2[5]=yy1[4]*IM[1]+yy1[3]*IM[2]+yy1[2]*IM[3]+yy1[1]*IM[4]+yy2[5];
        yy2[6]=yy1[5]*IM[1]+yy1[4]*IM[2]+yy1[3]*IM[3]+yy1[2]*IM[4]+yy1[1]*IM[5]+yy2[6];
        yy2[7]=yy1[6]*IM[1]+yy1[5]*IM[2]+yy1[4]*IM[3]+yy1[3]*IM[4]+yy1[2]*IM[5]+yy1[1]*IM[6]+yy2[7];

      Rec[i]=dot_product(Ct[1:(i-1)],IR[(1500-i+2):1500]);
      if(momento1[i]==1){
        IFR[i]=IFR0*2*inv_logit(-weekly_effect_IFR[1]-alpa_IFR_CFR*CFR[i]-alpa_IFR_s*(Total_V[i])/N-alpa_IFR_c*log(Ctacum[i]-Rec[i]+1)-
      momento_IFRt[1]*momento1[i]-momento_IFRt[2]*momento2[i]-momento_IFRt[3]*momento3[i]-
      momento_IFRt[4]*momento4[i]+dot_product(-alpa_IFR_j_M,Variantes[i,]));
      }
      if(momento2[i]==1){
        IFR[i]=IFR0*2*inv_logit(-weekly_effect_IFR[2]-alpa_IFR_CFR*CFR[i]-alpa_IFR_s*(Total_V[i])/N-alpa_IFR_c*log(Ctacum[i]-Rec[i]+1)-
      momento_IFRt[1]*momento1[i]-momento_IFRt[2]*momento2[i]-momento_IFRt[3]*momento3[i]-
      momento_IFRt[4]*momento4[i]+dot_product(-alpa_IFR_j_M,Variantes[i,]));
      }
      if(momento3[i]==1){
        IFR[i]=IFR0*2*inv_logit(-weekly_effect_IFR[3]-alpa_IFR_CFR*CFR[i]-alpa_IFR_s*(Total_V[i])/N-alpa_IFR_c*log(Ctacum[i]-Rec[i]+1)-
      momento_IFRt[1]*momento1[i]-momento_IFRt[2]*momento2[i]-momento_IFRt[3]*momento3[i]-
      momento_IFRt[4]*momento4[i]+dot_product(-alpa_IFR_j_M,Variantes[i,]));
      }
      if(momento4[i]==1){
        IFR[i]=IFR0*2*inv_logit(-weekly_effect_IFR[4]-alpa_IFR_CFR*CFR[i]-alpa_IFR_s*(Total_V[i])/N-alpa_IFR_c*log(Ctacum[i]-Rec[i]+1)-
      momento_IFRt[1]*momento1[i]-momento_IFRt[2]*momento2[i]-momento_IFRt[3]*momento3[i]-
      momento_IFRt[4]*momento4[i]+dot_product(-alpa_IFR_j_M,Variantes[i,]));
      }
      
      
      dt[i]=IFR[i]*sum(yy2);
      if(dt[i]<0.0001){dt[i]=0.0001;}
      dtacum[i]=dt[i]+dtacum[i-1];
    
  }
  borde1=Ctacum[16]-2330;
  borde2=Ctacum[37]-7798;
  borde3=Ctacum[93]-964035;
  borde4=Ctacum[117]-1105443;
}



model {
  
  kapa~normal(0,0.1);
  R0~normal(1.40,kapa);

  kapa_IFR0~normal(0,0.0001);
  IFR0~normal(0.01,kapa_IFR0);

  // jota_s~normal(0,0.1);
  // alpa_s~normal(0,jota_s);
  // 
  // jota_c~normal(0,0.1);
  // alpa_c~normal(0,jota_c);
  
  for(k in 1:nv){
    alpa_j_M[k]~normal(0,0.5);
    alpa_IFR_j_M~normal(0,0.5);
  }
  
  for (k in 1:4){
    j_Rt[k]~normal(0,0.5);
    momento_Rt[k]~normal(0,j_Rt[k]);
    j_IFRt[k]~normal(0,0.5);
    momento_IFRt[k]~normal(0,j_IFRt[k]);
    
  }
  
  j_vac_susep~normal(0,0.09);;
  vac_susep~normal(0.9,j_vac_susep);;
  
  
  // 
  // j_salva~normal(0,0.1);
  // salva_Rt~ normal(0,j_salva);
  

  jota_IFR_s~normal(0,0.1);
  alpa_IFR_s~normal(0,jota_IFR_s);

  jota_IFR_c~normal(0,0.1);
  alpa_IFR_c~normal(0,jota_IFR_c);
  jota_IFR_CFR~normal(0,0.5);
  alpa_IFR_CFR~normal(0,jota_IFR_CFR);
  
  jota_alpa_error~normal(0,0.1);
  alpa_error~normal(1,jota_alpa_error);
  jota_alpa_error_Ct~normal(0,0.1);
  alpa_error_Ct~normal(1,jota_alpa_error_Ct);
  
  
  tau ~ exponential(0.03);
  fi~normal(0,5);
  c_inicial~exponential(1/tau);
  
  weekly_effect_R0[1] ~ normal(0, 0.01);
  weekly_effect_R0[2] ~ normal(0,weekly_sd_R0 *sqrt(1-pow(weekly_rho_R0,2)-pow(weekly_rho1_R0,2) - 2 * pow(weekly_rho_R0,2) * weekly_rho1_R0/(1-weekly_rho1_R0)));
  weekly_effect_R0[3:(4+1)] ~ normal( weekly_effect_R0[2:4]* weekly_rho_R0 + weekly_effect_R0[1:(4-1)]* weekly_rho1_R0, 
                                            weekly_sd_R0 *sqrt(1-pow(weekly_rho_R0,2)-pow(weekly_rho1_R0,2) - 2 * pow(weekly_rho_R0,2) * weekly_rho1_R0/(1-weekly_rho1_R0)));
  
  
  weekly_effect_IFR[1] ~ normal(0, 0.01);
  weekly_effect_IFR[2] ~ normal(0,weekly_sd_IFR *sqrt(1-pow(weekly_rho_IFR,2)-pow(weekly_rho1_IFR,2) - 2 * pow(weekly_rho_IFR,2) * weekly_rho1_IFR/(1-weekly_rho1_IFR)));
  weekly_effect_R0[3:(4+1)] ~ normal( weekly_effect_IFR[2:4]* weekly_rho_IFR + weekly_effect_IFR[1:(4-1)]* weekly_rho1_IFR, 
                                            weekly_sd_IFR *sqrt(1-pow(weekly_rho_IFR,2)-pow(weekly_rho1_IFR,2) - 2 * pow(weekly_rho_IFR,2) * weekly_rho1_IFR/(1-weekly_rho1_IFR)));
  

  target+=normal_lpdf(IFR[1]|IFR0,0.0005);
  target+=normal_lpdf(Rt[1]|R0,0.0001);

  target+=normal_lpdf(IFR[16]-IFR[17]|0,0.0001);
  target+=normal_lpdf(IFR[37]-IFR[38]|0,0.0001);
  target+=normal_lpdf(IFR[93]-IFR[94]|0,0.0001);

  target+=normal_lpdf(Rt[16]-Rt[17]|0,0.0001);
  target+=normal_lpdf(Rt[37]-Rt[38]|0,0.0001);
  target+=normal_lpdf(Rt[93]-Rt[94]|0,0.0001);

  // Muertes_datos ~ neg_binomial_2(dt,fi);
   Muertes_datos ~ poisson(dt);
}

