# Marcos Temple
# 25/05/2021
close all
clear all
clc
# Troca de gases

# 1- DEFINIÇÕES INICIAIS

  # Parâmetros
  f = 0.25;                #[rad/s]
  #f = 0.83;                #[AUMENTADO] rad/s] - 50 inc/min
  N = 60000;
  dt = 0.001;                        #[s]
  t = 0:dt:((N-1)*dt);              #[s]
  Pfis = 1.5;                       #Amplitude de Pw
  phi = 0;                          #Fase
  
  # Unidades STPD
  f_O2 = 0.2094;                  # fração do gas na atm
  f_CO2 = 0.0038;                 # fração do gas na atm
  f_N = 0.7868;                   # fração do gas na atm
  VA_t = 2.2/1000;                # volume alveolar antes da inspiração[L] 2.2 l - 0.0022 m3
  Patm = 100000;                 # 100000 Pa - 1 atm - 760 [mmHg]
  R = 8.314;                      # Cte universal gases ideais Pa.m3/mol.K
  T = 273+36.5;                   # Temperatura corporal [K]  
  nT = ((Patm*VA_t)/(R*T))*1000; # numero de mols total N2, O2, CO2 alveolo[mmol]  
  P_inO2 = Patm*f_O2;
  P_inCO2 = Patm*f_CO2;
  P_inN2 = Patm*f_N;          fprintf('P_inO2: %.3f Pa, P_inCO2: %.3f Pa, P_inN2: %.3f Pa, nT: %.4f mmol \n\n',P_inO2, P_inCO2, P_inN2, nT)
  x_O2 = P_inO2/Patm;        # fração molar do gas
  x_CO2 = P_inCO2/Patm;      # fração molar do gas
  x_N2 = P_inN2/Patm;        # fração molar do gas
  
  #x0 = [nT*f_O2;nT*f_CO2;nT*f_N;      # x = [ n_A_O2; n_A_CO2; n_A_N; 
  x0 = [nT*f_O2;nT*f_CO2;nT*f_N;      # x = [ n_A_O2; n_A_CO2; n_A_N; 
        0.52;0.014;                         #       n_cap_O2; n_cap_CO2;  
        191;76];                         #       n_t_O2; n_t_CO2 ]
  x = zeros(7,N);         
  x(:,1) = x0;
  y = zeros(5,N);
  u_plot = zeros(1,N);
  nT_plot = zeros(1,N);
  nT_plot(:,1) = x0(1)+x0(2)+x0(3);
  PA = zeros(3,N);
  PA(:,1) = [nT*f_O2*((Patm)/nT);
             nT*f_CO2*((Patm)/nT);
             nT*f_N*((Patm)/nT)];

# 2- INTEGRAÇÃO NUMÉRICA - Runge-Kutta 4ªOrdem
  for i=1:(N-1)
    if t(i) == 3.689
      x(:,i);
    end
    phi = 2*pi*f*t(i);
    u = entrada_compartimental(phi, Pfis);  #QA*(Patm*x/R*T)
    u_plot(i) = u; 
    
    phi_umeio = phi+2*pi*f*(t(i)+dt/2);
    u_meio = entrada_compartimental(phi_umeio, Pfis);
    
    phi_udt = phi+2*pi*f*(t(i+1));
    u_dt = entrada_compartimental(phi_udt, Pfis);
    
    % calc 
    k1 = derivada_compartimental_A(x(:,i), u);
    k2 = derivada_compartimental_A(x(:,i)+k1*dt/2, u_meio);
    k3 = derivada_compartimental_A(x(:,i)+k2*dt/2, u_meio);
    k4 = derivada_compartimental_A(x(:,i)+k3*dt, u_dt);
    
    % Calcula estado em t+dt
    x(:,i+1) = x(:,i) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;  
    
    phi = phi_udt;
    nT_plot(i+1) = x(1,i+1)+x(2,i+1)+x(3,i+1);    # Calc nT
    
    PA(:, i+1) = [x(1,i+1)*((Patm)/nT_plot(i+1));
                  x(2,i+1)*((Patm)/nT_plot(i+1));
                  x(3,i+1)*((Patm)/nT_plot(i+1))]; 
    
    
    
  end #for rk

# 3- PLOTAGEM DOS RESULTADOS
  n_A_T = nT_plot;
  n_A_O2 = x(1,:);
  n_A_CO2 = x(2,:);
  n_A_N = x(3,:);
  n_cap_O2 = x(4,:);
  n_cap_CO2 = x(5,:);
  n_t_O2 = x(6,:);
  n_t_CO2 = x(7,:);
  
  alfa_O2 = 9.83e-06*1000; # mmol/m3.Pa
  alfa_CO2 = (0.0307 + 0.00057*(37-(T-273))+0.00002*(37-(T-273))^2)*1333.22; #mmol/(m3.mmHg)
  
  Vcap = 0.1/1000;       # m3
  n_cap_T = (n_A_T.*(-n_cap_O2-n_cap_CO2))/(n_A_N-n_A_T);
  P_cap_O2 = Patm*(n_cap_O2/n_cap_T)*0.00750062;
  #P_cap_O2 = (n_cap_O2/Vcap)/alfa_O2;
  P_cap_CO2 = Patm*(n_cap_CO2/n_cap_T)*0.00750062;
  #P_cap_CO2 = (n_cap_CO2/Vcap)/alfa_CO2;
  
  # Alvéolo
  figure(1)
  subplot(2,1,1)
  plot(t, n_A_O2,'r', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 12) 
  xlabel("Tempo [s]");
  ylabel("nº mols [mmol]");
  legend("n O2");
  title("nº de Mols dos Gases no Compartimento Alveolar - Função Senoidal com Corrida");
  
  figure(1)
  subplot(2,1,2)
  plot(t, n_A_CO2,'b', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 12) 
  xlabel("Tempo [s]");
  ylabel("nº mols [mmol]");
  legend("n CO2");
  print -djpg _mol_alv_fsen_corr.jpg;

#  figure(1)
#  subplot(3,1,3)
 # plot(t, n_A_N, 'k', 'linewidth', 3);
  #set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 
  #xlabel("Tempo [s]");
  #ylabel("nº mols [mmol]");
 # legend("n N2");
 # print -djpg mol_cap_apn_corr.jpg;
  
  # Capilares
  figure(2)
  subplot(2,1,1)
  plot(t, n_cap_O2,'r', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 12) 
  xlabel("Tempo [s]");
  ylabel("nº mols [mmol]");
  legend("n O2");
  title("nº de Mols dos Gases no Compartimento dos Capilares - Função Senoidal com Corrida");
  
  figure(2)
  subplot(2,1,2)
  plot(t, n_cap_CO2,'b', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 12) 
  xlabel("Tempo [s]");
  ylabel("nº mols [mmol]");
  legend("n CO2");
  print -djpg _mol_cap_fsen_corr.jpg;
  
  # Tecidos
  figure(3)
  subplot(2,1,1)
  plot(t, n_t_O2,'r', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize',12) 
  xlabel("Tempo [s]");
  ylabel("nº mols [mmol]");
  legend("n O2");
  title("nº de Mols dos Gases no Compartimento Tecidual - Função Senoidal com Corrida");
  
  figure(3)
  subplot(2,1,2)
  plot(t, n_t_CO2,'b', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 12) 
  xlabel("Tempo [s]");
  ylabel("nº mols [mmol]");
  legend("n CO2");
  print -djpg _mol_tec_fsen_corr.jpg;
  
  # Número de mols total e ventilação alveolar
  figure(4)
  plot(t, nT_plot, 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 12) 
  xlabel("Tempo [s]");
  ylabel("nº mols [mmol]");
  legend("nT");
  #title("Número de mols total - Alvéolo");

  #figure(6)
  #plot(t, u_plot,'r');
  #xlabel("Tempo [s]");
  #ylabel("Ventilação [m^3/s]");
  #legend("QA");
  #title("Ventilação Alveolar");
  
  # Pressões parciais
  figure(5)
  subplot(2,1,1)
  plot(t, PA(1,:)*0.00750062, 'r', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize',12) 
  xlabel("Tempo [s]");
  ylabel("P parcial [mmHg]");
  legend("P O2");
  title("Pressões Parciais do Ar nos Alvéolos - Função Senoidal com Corrida");
  
  figure(5)
  subplot(2,1,2)
  plot(t, PA(2,:)*0.00750062, 'b', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize',12) 
  xlabel("Tempo [s]");
  ylabel("P parcial [mmHg]");
  legend("P CO2");
  print -djpg _pressao_alv_fsen_corr.jpg;
  
  #figure(5)
  #subplot(3,1,3)
  #plot(t, PA(3,:)*0.00750062, 'k');
  #xlabel("Tempo [s]");
  #ylabel("P parcial [mmHg]");
  #legend("P N2");
  
  # Pressões parciais
  figure(6)
  subplot(2,1,1)
  plot(t, P_cap_O2, 'r', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 12) 
  xlabel("Tempo [s]");
  ylabel("P parcial [mmHg]");
  legend("P O2");
  title("Pressões Parciais do Ar nos Capilares - Função Senoidal com Corrida");
  
  figure(6)
  subplot(2,1,2)
  plot(t, P_cap_CO2, 'b', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 12) 
  xlabel("Tempo [s]");
  ylabel("P parcial [mmHg]");
  legend("P CO2");
  print -djpg _pressao_cap_fsen_corr.jpg;