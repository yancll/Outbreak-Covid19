
data {
  int<lower=0> n;
  int<lower=0> m;
  // int<lower=0> N;
  int<lower=0> nv;
  // int Casos_datos[n];
  int Muertes_datos[n];
  // real Nuevos_V[n];
  // real L_inf_dia[1500];
  // real L_sup_dia[1500];
  // real Total_V[n];
  // matrix [n,nv] Variantes;
}

transformed data {
  real IS[1500]; // SI in reverse order
  real IM[1500]; // SI in reverse order
  real IR[1500]; // SI in reverse order
  // for (i in 1:1500){
  //    IS[i]=gamma_cdf(L_sup_dia[i], 5.3*(5.3/0.26)/0.26,(5.3/0.26)/0.26)-gamma_cdf(L_inf_dia[i], 5.3*(5.3/0.26)/0.26,(5.3/0.26)/0.26);
  //    IM[i]=gamma_cdf(L_sup_dia[i], 21.33660,1/1.03174)-gamma_cdf(L_inf_dia[i], 21.33660,1/1.03174);
  //    IR[i]=gamma_cdf(L_sup_dia[i], 18.62067,1/1.18892)-gamma_cdf(L_inf_dia[i], 18.62067,1/1.18892);
  //
  // }

}


parameters {
  real<lower=0> tau;
  real lamda;
  real R0;
  real IFR0;
  real<lower=0> kapa;
  real<lower=0> kapa_IFR0;
   real<lower=0> kapa_lamda;
  real<lower=0> a_v;
  real<lower=0> b_v;
  real<lower=0> alpa_techo_v;
  real<lower=0> a_c;
  real<lower=0> b_c;
  real<lower=0> alpa_techo_c;
  real<lower=0> b_techo_v;
  real<lower=0> b_techo_c;
  real<lower=0> fi;
  real f11;
  real<lower=0> c_inicial;
  
  real alpa_s;
  real alpa_c;
  real<lower=0> jota_s;
  real<lower=0> jota_c;
  
  real alpa_IFR_s;
  real alpa_IFR_c;
  real<lower=0> jota_IFR_s;
  real<lower=0> jota_IFR_c;
  
  real  alpa_j_M[nv];
  real  alpa_IFR_j_M[nv];
  
}

transformed parameters {
  real  St[n];
  real  At[n];
  real  Ct[n];
  vector  [n]dt;
  vector  [n]dt_cuad;
  real  Vtecho_t[n];
  real  Ctecho_t[n];
  real  Vtecho_aster_t[n];
  real  Ctecho_aster_t[n];
  real  Rt[n];
  real  IFR[n];
  real y;
  real y1[5];
  real yy1[7];
  real yy2[7];
  real yy3[7];
  real yy4[7];
  real yy5[7];
  real IV[1500];
  real IC[1500];
  vector [n]salida;
   
  for (i in 1:n){
    if(i==1){ 
      // St[i]=1;Vtecho_t[i]=0;Ctecho_t[i]=0;
      // y=0;
      // for(j in 1:nv){
      //   y=y+alpa_j_M[j]*Variantes[i,j];
      // }
      // Rt[i]=R0;
      // y=0;
      // for(j in 1:nv){
      //   y=y+alpa_IFR_j_M[j]*Variantes[i,j];
      // }
      // IFR[i]=IFR0*inv_logit(alpa_IFR_s*Total_V[i]+alpa_IFR_c*(Casos_datos[i])+y);
      // Ct[i]=c_inicial;
      // At[i]=Ct[i];
      // dt[i]=0;
      // salida[i]=(dt[i]*(IFR[i]*At[i])/N)/(dt[i]-(IFR[i]*At[i])/N);
    }
   if(i>1){
     
     //  for (j in 1:1500){
     //    IV[i]=gamma_cdf(L_sup_dia[i], a_v,b_v)-gamma_cdf(L_inf_dia[i], a_v,b_v);
     //    IC[i]=gamma_cdf(L_sup_dia[i], a_c,b_c)-gamma_cdf(L_inf_dia[i], a_c,b_c);
     // }
     //  yy4=rep_array(0,7);
     //  yy5=rep_array(0,7);
     //  for (j in 1:(i-1)) {
     //    yy4[1]=sum(IV[(7*(j-1)+1):(7*(i-1)+7)])*Ct[i-j]/7+yy4[1];
     //    yy4[2]=sum(IV[(7*(j-1)+2):(7*(i-1)+8)])*Ct[i-j]/7+yy4[2];
     //    yy4[3]=sum(IV[(7*(j-1)+3):(7*(i-1)+9)])*Ct[i-j]/7+yy4[3];
     //    yy4[4]=sum(IV[(7*(j-1)+4):(7*(i-1)+10)])*Ct[i-j]/7+yy4[4];
     //    yy4[5]=sum(IV[(7*(j-1)+5):(7*(i-1)+11)])*Ct[i-j]/7+yy4[5];
     //    yy4[6]=sum(IV[(7*(j-1)+6):(7*(i-1)+12)])*Ct[i-j]/7+yy4[6];
     //    yy4[7]=sum(IV[(7*(j-1)+7):(7*(i-1)+13)])*Ct[i-j]/7+yy4[7];
     //    
     //    yy5[1]=sum(IC[(7*(j-1)+1):(7*(i-1)+7)])*Ct[i-j]/7+yy5[1];
     //    yy5[2]=sum(IC[(7*(j-1)+2):(7*(i-1)+8)])*Ct[i-j]/7+yy5[2];
     //    yy5[3]=sum(IC[(7*(j-1)+3):(7*(i-1)+9)])*Ct[i-j]/7+yy5[3];
     //    yy5[4]=sum(IC[(7*(j-1)+4):(7*(i-1)+10)])*Ct[i-j]/7+yy5[4];
     //    yy5[5]=sum(IC[(7*(j-1)+5):(7*(i-1)+11)])*Ct[i-j]/7+yy5[5];
     //    yy5[6]=sum(IC[(7*(j-1)+6):(7*(i-1)+12)])*Ct[i-j]/7+yy5[6];
     //    yy5[7]=sum(IC[(7*(j-1)+7):(7*(i-1)+13)])*Ct[i-j]/7+yy5[7];
     //    
     //  }
     //    
     //  
     //  Vtecho_t[i]=alpa_techo_v*sum(yy4);
     //  Ctecho_t[i]=alpa_techo_c*sum(yy5);
     //  Vtecho_aster_t[i-1]=Vtecho_t[i-1]/(N*St[i-1])*Ct[i-1];
     //  Ctecho_aster_t[i-1]=Ctecho_t[i-1]/(N*St[i-1])*Ct[i-1];
     //  St[i]=1-sum(Ct[1:(i-1)])/N-Total_V[i-1]/N-sum(dt[1:(i-1)])/N -sum(Vtecho_t[1:i])/N-sum(Ctecho_t[1:i])/N+sum(Vtecho_aster_t[1:(i-1)])/N+sum(Ctecho_aster_t[1:(i-1)])/N;
     //  y=0;
     //  for(j in 1:nv){
     //    y=y+alpa_j_M[j]*Variantes[i,j];
     //  }
     //  Rt[i]=R0*inv_logit(alpa_s*St[i]+alpa_c*Casos_datos[i]/Casos_datos[i-1]+y);
     //  yy1=rep_array(0,7);
     //  yy2=rep_array(0,7);
     //  yy3=rep_array(0,7);
     //  for (j in 1:(i-1)) {
     //    yy1[1]=sum(IS[(7*(j-1)+1):(7*(i-1)+7)])*Ct[i-j]/7+yy1[1];
     //    yy1[2]=sum(IS[(7*(j-1)+2):(7*(i-1)+8)])*Ct[i-j]/7+yy1[2];
     //    yy1[3]=sum(IS[(7*(j-1)+3):(7*(i-1)+9)])*Ct[i-j]/7+yy1[3];
     //    yy1[4]=sum(IS[(7*(j-1)+4):(7*(i-1)+10)])*Ct[i-j]/7+yy1[4];
     //    yy1[5]=sum(IS[(7*(j-1)+5):(7*(i-1)+11)])*Ct[i-j]/7+yy1[5];
     //    yy1[6]=sum(IS[(7*(j-1)+6):(7*(i-1)+12)])*Ct[i-j]/7+yy1[6];
     //    yy1[7]=sum(IS[(7*(j-1)+7):(7*(i-1)+13)])*Ct[i-j]/7+yy1[7];
     //    
     //    yy2[1]=sum(IM[(7*(j-1)+1):(7*(i-1)+7)])*Ct[i-j]/7+yy2[1];
     //    yy2[2]=sum(IM[(7*(j-1)+2):(7*(i-1)+8)])*Ct[i-j]/7+yy2[2];
     //    yy2[3]=sum(IM[(7*(j-1)+3):(7*(i-1)+9)])*Ct[i-j]/7+yy2[3];
     //    yy2[4]=sum(IM[(7*(j-1)+4):(7*(i-1)+10)])*Ct[i-j]/7+yy2[4];
     //    yy2[5]=sum(IM[(7*(j-1)+5):(7*(i-1)+11)])*Ct[i-j]/7+yy2[5];
     //    yy2[6]=sum(IM[(7*(j-1)+6):(7*(i-1)+12)])*Ct[i-j]/7+yy2[6];
     //    yy2[7]=sum(IM[(7*(j-1)+7):(7*(i-1)+13)])*Ct[i-j]/7+yy2[7];
     //    
     //    yy3[1]=sum(IR[(7*(j-1)+1):(7*(i-1)+7)])*Ct[i-j]/7+yy3[1];
     //    yy3[2]=sum(IR[(7*(j-1)+2):(7*(i-1)+8)])*Ct[i-j]/7+yy3[2];
     //    yy3[3]=sum(IR[(7*(j-1)+3):(7*(i-1)+9)])*Ct[i-j]/7+yy3[3];
     //    yy3[4]=sum(IR[(7*(j-1)+4):(7*(i-1)+10)])*Ct[i-j]/7+yy3[4];
     //    yy3[5]=sum(IR[(7*(j-1)+5):(7*(i-1)+11)])*Ct[i-j]/7+yy3[5];
     //    yy3[6]=sum(IR[(7*(j-1)+6):(7*(i-1)+12)])*Ct[i-j]/7+yy3[6];
     //    yy3[7]=sum(IR[(7*(j-1)+7):(7*(i-1)+13)])*Ct[i-j]/7+yy3[7];
     //  }
     //    yy1[1]=yy1[1]*Rt[i];
     //    yy1[2]=yy1[1]*IS[1]+yy1[2];
     //    yy1[2]=yy1[2]*Rt[i];
     //    yy1[3]=yy1[2]*IS[1]+yy1[1]*IS[2]+yy1[3];
     //    yy1[3]=yy1[3]*Rt[i];
     //    yy1[4]=yy1[3]*IS[1]+yy1[2]*IS[2]+yy1[1]*IS[3]+yy1[4];
     //    yy1[4]=yy1[4]*Rt[i];
     //    yy1[5]=yy1[4]*IS[1]+yy1[3]*IS[2]+yy1[2]*IS[3]+yy1[1]*IS[4]+yy1[5];
     //    yy1[5]=yy1[5]*Rt[i];
     //    yy1[6]=yy1[5]*IS[1]+yy1[4]*IS[2]+yy1[3]*IS[3]+yy1[2]*IS[4]+yy1[1]*IS[5]+yy1[6];
     //    yy1[6]=yy1[6]*Rt[i];
     //    yy1[7]=yy1[6]*IS[1]+yy1[5]*IS[2]+yy1[4]*IS[3]+yy1[3]*IS[4]+yy1[2]*IS[5]+yy1[1]*IS[6]+yy1[7];
     //    yy1[7]=yy1[7]*Rt[i];
     //    Ct[i]=sum(yy1)*St[i];
     //    
     //  
     //    yy2[2]=yy1[1]*IM[1]+yy2[2];
     //    yy2[3]=yy1[2]*IM[1]+yy1[1]*IM[2]+yy2[3];
     //    yy2[4]=yy1[3]*IM[1]+yy1[2]*IM[2]+yy1[1]*IM[3]+yy2[4];
     //    yy2[5]=yy1[4]*IM[1]+yy1[3]*IM[2]+yy1[2]*IM[3]+yy1[1]*IM[4]+yy2[5];
     //    yy2[6]=yy1[5]*IM[1]+yy1[4]*IM[2]+yy1[3]*IM[3]+yy1[2]*IM[4]+yy1[1]*IM[5]+yy2[6];
     //    yy2[7]=yy1[6]*IM[1]+yy1[5]*IM[2]+yy1[4]*IM[3]+yy1[3]*IM[4]+yy1[2]*IM[5]+yy1[1]*IM[6]+yy2[7];
     //    
     //    yy3[2]=yy1[1]*IR[1]+yy3[2];
     //    yy3[3]=yy1[2]*IR[1]+yy1[1]*IR[2]+yy3[3];
     //    yy3[4]=yy1[3]*IR[1]+yy1[2]*IR[2]+yy1[1]*IR[3]+yy3[4];
     //    yy3[5]=yy1[4]*IR[1]+yy1[3]*IR[2]+yy1[2]*IR[3]+yy1[1]*IR[4]+yy3[5];
     //    yy3[6]=yy1[5]*IR[1]+yy1[4]*IR[2]+yy1[3]*IR[3]+yy1[2]*IR[4]+yy1[1]*IR[5]+yy3[6];
     //    yy3[7]=yy1[6]*IR[1]+yy1[5]*IR[2]+yy1[4]*IR[3]+yy1[3]*IR[4]+yy1[2]*IR[5]+yy1[1]*IR[6]+yy3[7];
     //    
     //  y=0;
     //  for(j in 1:nv){
     //    y=y+alpa_IFR_j_M[j]*Variantes[i,j];
     //  }
     //  IFR[i]=IFR0*inv_logit(alpa_IFR_s*Total_V[i]+alpa_IFR_c*(Casos_datos[i]-Casos_datos[i-1])+y);
     //  dt[i]=IFR[i]*sum(yy2) ;
         dt_cuad[i]=Muertes_datos[i];
     //  
     //  At[i]=sum(Ct[1:(i-1)])-(lamda)*sum(yy3);
    }
  }
}


model {
  // b_techo_v~normal(0,0.5);
  // alpa_techo_v~normal(0,b_techo_v);
  // 
  // a_v~normal(0,10);
  // b_v~normal(0,10);
  // 
  // alpa_techo_c~normal(0,b_techo_c);
  // b_techo_c~normal(0,0.5);
  a_c~normal(0,10);
  // b_c~normal(0,10);
  // 
  // kapa~normal(0,0.5);
  // R0~normal(1.40,kapa);
  // 
  // kapa_IFR0~normal(0,0.5);
  // IFR0~normal(0.01,kapa_IFR0);
  // 
  // jota_s~normal(0,0.5);
  // alpa_s~normal(0,jota_s);
  // 
  // jota_c~normal(0,0.5);
  // alpa_c~normal(0,jota_c);
  // 
  // alpa_j_M~normal(0,0.5);
  // 
  // kapa_lamda~normal(0,0.5);
  // lamda~normal(1/5.3,kapa_lamda);
  // 
  // alpa_IFR_j_M~normal(0,0.5);
  // 
  // jota_IFR_s~normal(0,0.5);
  // alpa_IFR_s~normal(0,jota_IFR_s);
  // 
  // jota_IFR_c~normal(0,0.5);
  // alpa_IFR_c~normal(0,jota_IFR_c);
  // tau ~ exponential(0.03);
  fi~normal(0,10);
  f11~normal(0,5);
  // c_inicial~exponential(1/tau);
  Muertes_datos ~ neg_binomial_2(dt_cuad+f11,fi);
}

