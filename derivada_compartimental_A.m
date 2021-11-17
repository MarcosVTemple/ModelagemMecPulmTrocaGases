function dx = derivada_compartimental_A(x, u)
  # x = [ n_A_O2; n_A_CO2; n_A_N; n_cap_O2; n_cap_CO2; n_t_O2; n_t_CO2]
  
  R = 8.314;                      # Cte universal gases ideais [Pa.m3/mol.K]
  T = 273+36.5;                   # Temperatura corporal [K]  
  VA_t = 2.2/1000;                # volume alveolar antes da inspiração 2.2 [L] - 0.0022 [m3]
  Patm = 100000;                 # 100000 [Pa] - 1 [atm] - 760 [mmHg]
 # Patm =  10000;                 # Teste pressão de referencia menor
  f_O2 = 0.2094;                  # fração do gas na atm
  f_CO2 = 0.0038;                 # fração do gas na atm
  f_N = 0.7868;                   # fração do gas na atm
  #D_O2 = 0.00002*133.322/1000;                 # capacidade de difusao do gas O2 atraves da memb. respiratoria [m3/s/mmHg]
  #D_CO2 = 0.003*133.322/1000;                # capacidade de difusao do gas CO2 atraves da memb. respiratoria [m3/s/mmHg]
  D_O2_Alb = 32.253e-10;                   # 26 [ml/min.mmHg] - 0.00043 [L/s.mmHg] - 32.253e-10 [m³/s.Pa]
  D_O2 = (D_O2_Alb)*((Patm)/(R*T))*1000;        # convertido em mmols
  D_CO2_Alb = 22.502e-09;                  # 180 [ml/min.mmHg] - 0.003 [L/s.mmHg] - 22.502e-09 [m³/s.Pa]
  D_CO2 = D_CO2_Alb*((Patm)/(R*T))*1000;        # convertido em mmols
  
  #derivada capilar
    Q_b = (5.6/60)/1000;         #5.6 L/min - 5.6/1000 m3
    #Q_b = Q_b*60 ;         #AUMENTADO APNEIA #5.6 L/min - 5.6/1000 m3
    #Q_b = (25.6/60)/1000;         #AUMENTADO SENOIDAL 25.6 L/min - 5.6/1000 m3 (Corrida)
    sigma = 0.98;
    Vt = 37.35/1000;       # Volume tecidos musculares/nao musc/capilares tec
    Vcap = 0.1/1000;       # m3
  
  #derivada tecidual
    n_t_O2_fis = 208;                 # número de mols do gás O2 nos tecidos (fisiológico) [mmol]
    n_t_CO2_fis = 76;                 # número de mols do gás CO2 nos tecidos (fisiológico) [mmol]
    c_t_O2_fis = 19607;               # n_t_O2/Vt   [mmol/m3]
    c_t_CO2_fis = 4610.1;              # n_t_CO2/Vt [mmol/m3]
    
    #Repouso
    Q_O2_Alb = ((0.245/1000)/60);          # proporção do consumo do gás O2 [m3/s] 200-300 ml/min
    Q_O2 = (Q_O2_Alb*(Patm/(R*T)))*1000;          # proporção do consumo do gás O2 [mmol/s] 
    #Aumentado
    #Q_O2_Alb_H = ((2.88/1000)/60);    # AUMENTADO: proporção do consumo do gás O2 [m3/s] 2880 ml/min
    #Q_O2 = (Q_O2_Alb_H*(Patm/(R*T)))*1000;          # AUMENTADO: proporção do consumo do gás O2 [mmol/s] 
    #Q_O2 = Q_O2*30000;          # AUMENTADO APNEIA: proporção do consumo do gás O2 [mmol/s] 
    
    Q_CO2 = Q_O2*0.85;                # proporção da produção do gás CO2 [mmol/s] - para QR = 0.85
 
  n_A_O2 = x(1);
  n_A_CO2 = x(2);
  n_A_N = x(3);
  n_cap_O2 = x(4);
  n_cap_CO2 = x(5);
  n_t_O2 = x(6);
  n_t_CO2 = x(7);
  
  nT = n_A_O2 + n_A_CO2 + n_A_N; # numero de mmols no ALVEOLO
  nT_cap = ((-n_cap_O2-n_cap_CO2)*(nT/n_A_N))/(1-(nT/n_A_N));
  
  u=0;        # Apneia
  # derivada alveolar
  if u >= 0
    dn_A_O2 = (D_O2*Patm)*(-(n_A_O2/nT)+(n_cap_O2/nT_cap)) + ((Patm*f_O2)/(R*T))*u;
    dn_A_CO2 = (D_CO2*Patm)*(-(n_A_CO2/nT)+(n_cap_CO2/nT_cap)) + ((Patm*f_CO2)/(R*T))*u;
    dn_A_N = 0 + ((Patm*f_N)/(R*T))*u;
    PA = Patm*(f_CO2+f_O2+f_N); 
  else      
    PA_O2 = n_A_O2*(Patm)/nT; 
    PA_CO2 = n_A_CO2*(Patm)/nT; 
    PA_N = n_A_N*(Patm)/nT;         
    PA = PA_O2+PA_CO2+PA_N;
    
    dn_A_O2 = (D_O2*Patm)*(-(n_A_O2/nT)+(n_cap_O2/nT_cap)) + ((PA_O2)/(R*T))*u;
    dn_A_CO2 = (D_CO2*Patm)*(-(n_A_CO2/nT)+(n_cap_CO2/nT_cap)) + ((PA_CO2)/(R*T))*u;
    dn_A_N = ((Patm*f_N)/(R*T))*u; 
  endif
  
  # derivada capilar   
    dn_cap_O2 =  (D_O2*Patm)*((n_A_O2/nT)-(n_cap_O2/nT_cap)) + (Q_b*sigma)*((n_t_O2/Vt)-(n_cap_O2/Vcap));
    dn_cap_CO2 =  (D_CO2*Patm)*((n_A_CO2/nT)-(n_cap_CO2/nT_cap)) + (Q_b*sigma)*((n_t_CO2/Vt)-(n_cap_CO2/Vcap));    
  
  # derivada tecidual
    K_O2 = Q_O2/c_t_O2_fis;
    K_CO2 = 0.85*K_O2;
    
    dn_t_O2 = Q_b*sigma*((-n_t_O2/Vt)+(n_cap_O2/Vcap)) -  K_O2*(n_t_O2/Vt);
    dn_t_CO2 = Q_b*sigma*((-n_t_CO2/Vt)+(n_cap_CO2/Vcap)) + K_CO2*(n_t_O2/Vt);
  
  dx = [dn_A_O2; dn_A_CO2; dn_A_N; dn_cap_O2; dn_cap_CO2; dn_t_O2; dn_t_CO2];
endfunction
