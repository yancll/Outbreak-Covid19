data {
  int<lower=0> n;
  int<lower=0> m;
  int<lower=0> N;
  int<lower=0> nv;
  vector [n]Casos_datos;
  int Casos_d_acum[n];
  int Muertes_datos[n];
  // vector [n]Muertes_datos;
  vector [n]Nuevos_V;
  vector [1500]IS; // SI in reverse order
  vector [1500]IM; // SI in reverse order
  vector [1500]IR; // SI in reverse order
  vector [1500]IV; // SI in reverse order
  vector [1500]IC;
  vector [4] lar;
  vector [n]Total_V;
  matrix [n,nv] Variantes;
  vector [n]CFR;
  vector [n] momento1;
  vector [n] momento2;
  vector [n] momento3;
  vector [n] momento4;
}


parameters {
  real<lower=0> tau;
  real<lower=0,upper=1> lamda;
  real<lower=0> R0;
  real<lower=0.007,upper=0.016> IFR0;
  real<lower=0> kapa;
  real<lower=0> kapa_IFR0;
  real<lower=0.5,upper=1> alpa_error;
  real<lower=0> jota_alpa_error;
  real salva_Rt;
  real<lower=0.5,upper=1> alpa_error_Ct;
  real<lower=0> jota_alpa_error_Ct;
  // real<lower=0,upper=1> alpa_techo_v;
  // real<lower=0,upper=1> alpa_techo_c;
  // real<lower=0> b_techo_v;
  // real<lower=0> b_techo_c;
  real alpa_error_combinado;
  real<lower=0>jota_alpa_error_combinado;
  real<lower=0> fi;
  real<lower=0> f11;
  real<lower=0>j_salva;
  real<lower=0,upper=500> c_inicial;
  real alpa_s;
  real alpa_c;
  real<lower=0> jota_s;
  real<lower=0> jota_c;
  real alpa_IFR_s;
  real alpa_IFR_c;
  real <lower=0>alpa_IFR_CFR;
  real <lower=0>jota_IFR_CFR;
  real<lower=0> kapa_lamda;
  real<lower=0> jota_IFR_s;
  real<lower=0> jota_IFR_c;
  // real<lower=0> jota_IFR_Act;
  vector [nv] alpa_j_M;
  vector  [nv]alpa_IFR_j_M;
  
  vector  [4] momento_Rt;
  vector<lower=0>  [4]j_Rt;
  vector  [4] momento_IFRt;
  vector<lower=0>  [4]j_IFRt;
}

// momento_Rt[1]*momento1[i]+momento_Rt[2]*momento1[i]+momento_Rt[3]*momento1[i]
// momento_IFRt[1]*momento1[i]+momento_IFRt[3]*momento1[i]+momento_IFRt[3]*momento1[i]


transformed parameters {
  vector<lower=0,upper=1> [n] St;
  vector<lower=0,upper=120100>  [n]Ct;
  vector<lower=0,upper=650>  [n]dt;
  vector<lower=0,upper=2000000>  [n]Ctacum;
  vector<lower=0,upper=9000>  [n]dtacum;
  // vector<upper=N>  [n]Vt_ac_2;
  // vector<upper=N>  [n]Vt_2;
  real prueba1;
  real prueba2;
  vector  [n]Ct_R_ac;
  vector  [n]Ct_R_ac_VC;
  // 
  // vector  [n]Ct_RS_ac;
  // vector  [n]Ct_RS;
  // 
  vector  [n]Vt_R_ac;
  // vector  [n]Vt_R;
  // 
  // vector  [n]Vt_RS_ac;
  // vector  [n]Vt_RS;

  // vector  [n]Vtecho_t;
  // vector  [n]Ctecho_t;
  // vector [n] Vtecho_aster_t;
  // vector [n]Ctecho_aster_t;
  vector<lower=0,upper=4> [n]Rt;
  vector<lower=0,upper=0.02> [n]IFR;
  vector<lower=0> [n]Rec;
  real y;
  // vector [5]y1;
  vector<lower=0> [7]yy1;
  vector<lower=0> [7]yy2;
  vector<lower=0> [7]yy3;
  // vector<lower=0> [7]yy4;
  // vector<lower=0> [7]yy5;
 
  St[1]=1;
  // Vtecho_t[1]=0;
  // Ctecho_t[1]=0;
  y=0;
  for(j in 1:nv){
    y=y+alpa_j_M[1]*Variantes[1,j];
  }
   // Rt[1]=R0*2*St[1]*inv_logit(+dot_product(alpa_j_M,Variantes[1,])+momento_Rt[1]*momento1[1]+momento_Rt[2]*momento2[1]+momento_Rt[3]*momento3[1]+momento_Rt[4]*momento4[1]);
  Rt[1]=R0*2*St[1]*inv_logit(+dot_product(alpa_j_M,Variantes[1,])+salva_Rt*R0+momento_Rt[1]*momento1[1]+momento_Rt[2]*momento2[1]+momento_Rt[3]*momento3[1]+momento_Rt[4]*momento4[1]);
  // Rt[1]=R0*2*inv_logit(alpa_s*St[1]+dot_product(alpa_j_M,Variantes[1,])+salva_Rt);
  // y=0;
  // for(j in 1:nv){
  //   y=y+alpa_IFR_j_M[j]*Variantes[1,j];
  // }
  Ct[1]=c_inicial;
  IFR[1]=IFR0*2*inv_logit(alpa_IFR_s*Total_V[1]+momento_IFRt[1]*momento1[1]+momento_IFRt[2]*momento2[1]+momento_IFRt[3]*momento3[1]+momento_IFRt[4]*momento4[1]+alpa_IFR_c*(Ct[1])+dot_product(alpa_IFR_j_M,Variantes[1,]))*alpa_IFR_CFR*CFR[1];
  
  
  dt[1]=0.0001;
  prueba1=inv_logit(alpa_s*St[1]+dot_product(alpa_j_M,Variantes[1,]));
  prueba2=inv_logit(alpa_IFR_s*Total_V[1]+alpa_IFR_c*(Ct[1])+dot_product(alpa_IFR_j_M,Variantes[1,]));
  Ctacum[1]=Ct[1];
  dtacum[1]=0;
  // Vt_ac_2[1]=0;
  // Vt_2[1]=0;
  // Vt_2[2]=0;
  // Vt_2[3]=0;
  // Vt_2[4]=0;
  // 
  Ct_R_ac[1]=0;
  // Ct_R[1]=0;
  // 
  // Ct_RS_ac[1]=0;
  // Ct_RS[1]=0;
  // 
  Vt_R_ac[1]=0;
  // Vt_R[1]=0;
  // 
  // Vt_RS_ac[1]=0;
  // Vt_RS[1]=0;
  dtacum[1]=0;
  Ct_R_ac_VC[1]=0;
  Ct_R_ac_VC[1]=0;
  Rec[1]=0;
  
  
  for (i in 2:n){
    // Vt_ac_2[i]=Total_V[i-1]*(Ctacum[i-1]-Casos_d_acum[i-1])/N;
    // Vt_2[i]=Vt_ac_2[i]-Vt_ac_2[i-1];
    
      // yy4=rep_vector(0,7);
      // yy5=rep_vector(0,7);
      // for (j in 1:(i-1)) {
      //   if((7*(j-1)+7)<(1500-lar[4]+7)){
      //     yy4[1]+=sum(IV[(7*(j-1)+1):(7*(j-1)+7)])*Ct[i-j]/7;
      //     // yy4[2]+=sum(ICV[(7*(j-1)+2):(7*(j-1)+8)])*Ct[i-j]/7;
      //     // yy4[3]+=sum(ICV[(7*(j-1)+3):(7*(j-1)+9)])*Ct[i-j]/7;
      //     // yy4[4]+=sum(ICV[(7*(j-1)+4):(7*(j-1)+10)])*Ct[i-j]/7;
      //     // yy4[5]+=sum(ICV[(7*(j-1)+5):(7*(j-1)+11)])*Ct[i-j]/7;
      //     // yy4[6]+=sum(ICV[(7*(j-1)+6):(7*(j-1)+12)])*Ct[i-j]/7;
      //     // yy4[7]+=sum(ICV[(7*(j-1)+7):(7*(j-1)+13)])*Ct[i-j]/7;
      // 
      //     yy5[1]+=sum(IV[(7*(j-1)+1):(7*(j-1)+7)])*Nuevos_V[i-j]/7;
      //     // yy5[2]+=sum(ICV[(7*(j-1)+2):(7*(j-1)+8)])*Nuevos_V[i-j]/7;
      //     // yy5[3]+=sum(ICV[(7*(j-1)+3):(7*(j-1)+9)])*Nuevos_V[i-j]/7;
      //     // yy5[4]+=sum(ICV[(7*(j-1)+4):(7*(j-1)+10)])*Nuevos_V[i-j]/7;
      //     // yy5[5]+=sum(ICV[(7*(j-1)+5):(7*(j-1)+11)])*Nuevos_V[i-j]/7;
      //     // yy5[6]+=sum(ICV[(7*(j-1)+6):(7*(j-1)+12)])*Nuevos_V[i-j]/7;
      //     // yy5[7]+=sum(ICV[(7*(j-1)+7):(7*(j-1)+13)])*Nuevos_V[i-j]/7;
      //   }
      // }

      
      
      Vt_R_ac[i]=dot_product(Nuevos_V[1:(i-1)],IC[(1500-i+2):1500]);
      // Vt_R[i]=Vt_R_ac[i]-Vt_R_ac[i-1];
      Ct_R_ac[i]=dot_product(Ct[1:(i-1)],IC[(1500-i+2):1500]);
      Ct_R_ac_VC[i]=dot_product(Ct[1:(i-1)]-Ct[1:(i-1)].*Nuevos_V[1:(i-1)]/N,IC[(1500-i+2):1500]);
      // 
      // Ct_RS_ac[i]=Ct_R_ac[i]-Vt_ac_2[i]-dtacum[i-1];
      // Ct_RS[i]=Ct_RS_ac[i]-Ct_RS_ac[i-1];
      // Vt_RS_ac[i]=Vt_R_ac[i]-Vt_ac_2[i-1];
      
      St[i]=(N+alpa_error_Ct*Ct_R_ac_VC[i]+alpa_error*Vt_R_ac[i]-Total_V[i-1]-Ctacum[i-1]+Total_V[i-1]*(Ctacum[i-1])/N)/N;
      // St[i]=(N-Ctacum[i-1]-Total_V[i-1]+alpa_error*Vt_RS_ac[i-1]+Ct_RS_ac[i-1])/N;
      // St[i]=(N-Ctacum[i-1]-Total_V[i-1]+Vt_ac_2[i-1]+Ct_RS_ac[i-1]+Vt_RS_ac[i-1])/N;
      // St[i]=(N)/N;
      // Vtecho_aster_t[i-1]=Vtecho_t[i-1]/(St[i-1]*N)*Ct[i-1];
      // Ctecho_aster_t[i-1]=Ctecho_t[i-1]/(N*St[i-1])*Ct[i-1];
      // St[i]=1-sum(Ct[1:(i-1)])/N-Total_V[i-1]/N-sum(dt[1:(i-1)])/N +sum(Vtecho_t[1:i])/N+sum(Ctecho_t[1:i])/N-sum(Vtecho_aster_t[1:(i-1)])/N-sum(Ctecho_aster_t[1:(i-1)])/N;
      // St[i]=St[i-1]-Vtecho_t[i-1]/N-Ctecho_t[i-1]/N-Ct[i-1]/N+Vtecho_aster_t[i-1]/N+Ctecho_aster_t[i-1]/N;
      
      // y=0;
      // for(j in 1:nv){
      //   y=y+alpa_j_M[j]*Variantes[i,j];
      // }
      
      
      // if(i<n){Rt[i]=R0*2*St[i]*inv_logit(dot_product(alpa_j_M,Variantes[i,])+momento_Rt[1]*momento1[i]+momento_Rt[2]*momento2[i]+momento_Rt[3]*momento3[i]+momento_Rt[4]*momento4[i]);}
      // if(i==n){Rt[i]=R0*2*St[i]*inv_logit(dot_product(alpa_j_M,Variantes[i,])+momento_Rt[1]*momento1[i]+momento_Rt[2]*momento2[i]+momento_Rt[3]*momento3[i]+momento_Rt[4]*momento4[i]);}
       if(i<n){Rt[i]=R0*2*St[i]*inv_logit(dot_product(alpa_j_M,Variantes[i,])+momento_Rt[1]*momento1[i]+momento_Rt[2]*momento2[i]+momento_Rt[3]*momento3[i]+momento_Rt[4]*momento4[i]+salva_Rt*((Casos_datos[i]+Casos_datos[1+i])/(Casos_datos[i]+Casos_datos[i-1])));}
      if(i==n){Rt[i]=R0*2*St[i]*inv_logit(dot_product(alpa_j_M,Variantes[i,])+momento_Rt[1]*momento1[i]+momento_Rt[2]*momento2[i]+momento_Rt[3]*momento3[i]+momento_Rt[4]*momento4[i]+salva_Rt*((Casos_datos[i])/(Casos_datos[i-1])));}
      
       // Rt[i]=R0*2*inv_logit(alpa_s*St[i]+dot_product(alpa_j_M,Variantes[i,])+salva_Rt*(Casos_datos[i]/Casos_datos[i-1]));
      yy1=rep_vector(0,7);
      yy2=rep_vector(0,7);
      yy3=rep_vector(0,7);
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
        // if((7*(j-1)+7)<(1500-lar[3]+7)){
        //    yy3[1]+=sum(IR[(7*(j-1)+1):(7*(j-1)+7)])*Ct[i-j]/7;
        //    // yy3[2]+=sum(IR[(7*(j-1)+2):(7*(j-1)+8)])*Ct[i-j]/7;
        //    // yy3[3]+=sum(IR[(7*(j-1)+3):(7*(j-1)+9)])*Ct[i-j]/7;
        //    // yy3[4]+=sum(IR[(7*(j-1)+4):(7*(j-1)+10)])*Ct[i-j]/7;
        //    // yy3[5]+=sum(IR[(7*(j-1)+5):(7*(j-1)+11)])*Ct[i-j]/7;
        //    // yy3[6]+=sum(IR[(7*(j-1)+6):(7*(j-1)+12)])*Ct[i-j]/7;
        //    // yy3[7]+=sum(IR[(7*(j-1)+7):(7*(j-1)+13)])*Ct[i-j]/7;
        // }
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

        // yy3[2]=yy1[1]*IR[1]+yy3[2];
        // yy3[3]=yy1[2]*IR[1]+yy1[1]*IR[2]+yy3[3];
        // yy3[4]=yy1[3]*IR[1]+yy1[2]*IR[2]+yy1[1]*IR[3]+yy3[4];
        // yy3[5]=yy1[4]*IR[1]+yy1[3]*IR[2]+yy1[2]*IR[3]+yy1[1]*IR[4]+yy3[5];
        // yy3[6]=yy1[5]*IR[1]+yy1[4]*IR[2]+yy1[3]*IR[3]+yy1[2]*IR[4]+yy1[1]*IR[5]+yy3[6];
        // yy3[7]=yy1[6]*IR[1]+yy1[5]*IR[2]+yy1[4]*IR[3]+yy1[3]*IR[4]+yy1[2]*IR[5]+yy1[1]*IR[6]+yy3[7];

        Rec[i]=dot_product(Ct[1:(i-1)],IR[(1500-i+2):1500]);
      // At[i]=sum(Ct[1:i]-Rec[1:i]);
      IFR[i]=IFR0*2*inv_logit(alpa_IFR_s*(Total_V[i])/N+alpa_IFR_c*log(Ctacum[i]-Rec[i])+momento_IFRt[1]*momento1[i]+momento_IFRt[2]*momento2[i]+momento_IFRt[3]*momento3[i]+momento_IFRt[4]*momento4[i]+dot_product(alpa_IFR_j_M,Variantes[i,]))*alpa_IFR_CFR*CFR[i];
      // IFR[i]=IFR0*inv_logit(alpa_IFR_s*Total_V[i] +alpa_IFR_c*(Casos_datos[i]-Casos_datos[i-1])+sum(Variantes*alpa_IFR_j_M));
      // Recacum[i]+=Recacum[i-1];
      // Ctacum[i]+=Ct[i];
      // At[i]=Ctacum[i]-Rec[i];
      dt[i]=IFR[i]*sum(yy2);
      if(dt[i]<0.0001){dt[i]=0.0001;}
      dtacum[i]=dt[i]+dtacum[i-1];
    
  }
}



model {
  // b_techo_v~normal(0,0.5);
  // alpa_techo_v~normal(0,b_techo_v);
  // b_techo_c~normal(0,0.5);
  // alpa_techo_c~normal(0,b_techo_c);
  kapa~normal(0,0.1);
  R0~normal(1.40,kapa);

  kapa_IFR0~normal(0,0.01);
  IFR0~normal(0.01,kapa_IFR0);

  // target+=normal_lpdf(prueba2|0.5,0.005);
  target+=normal_lpdf(IFR[1]|IFR0,0.0005);
  target+=normal_lpdf(Rt[1]|R0,0.0001);
  // if(sum(momento1)>0){target+=normal_lpdf(dot_product(momento1,dt)/dot_product(momento1,Ct)|dot_product(momento1,IFR),0.001);}
  // if(sum(momento2)>0){target+=normal_lpdf(dot_product(momento2,dt)/dot_product(momento2,Ct)|dot_product(momento2,IFR),0.001);}
  // if(sum(momento3)>0){target+=normal_lpdf(dot_product(momento3,dt)/dot_product(momento3,Ct)|dot_product(momento3,IFR),0.001);}
  // if(sum(momento4)>0){target+=normal_lpdf(dot_product(momento4,dt)/dot_product(momento4,Ct)|dot_product(momento4,IFR),0.001);}
  
  
  
  
  
  
  
  jota_s~normal(0,0.1);
  alpa_s~normal(0,jota_s);

  jota_c~normal(0,0.1);
  alpa_c~normal(0,jota_c);
  
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
  

  kapa_lamda~normal(0,0.1);
  // lamda~normal(1,kapa_lamda);
  // target+=normal_lpdf(lamda|1,kapa_lamda);
  j_salva~normal(0,0.1);
  salva_Rt~ normal(0,j_salva);
  

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
  
  jota_alpa_error_combinado~normal(0,0.1);
  alpa_error_combinado~normal((alpa_error+alpa_error_Ct)/2,jota_alpa_error_combinado);
  
  // jota_IFR_Act~normal(0,0.5);
  // alpa_IFR_Act~normal(0,jota_IFR_Act);
  // 
  
  tau ~ exponential(0.03);
  fi~normal(0,5);
  f11~normal(0,0.00005);
  c_inicial~exponential(1/tau);
  // Muertes_datos ~ normal(dt,fi);
  // target+=neg_binomial_2_lpmf(Muertes_datos|dt,fi);
  Muertes_datos ~ neg_binomial_2(dt,fi);
}

















