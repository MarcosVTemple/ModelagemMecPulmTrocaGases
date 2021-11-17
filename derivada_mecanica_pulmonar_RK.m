function dx = derivada_mecanica_pulmonar_RK(x,t,u)
# x - condicoes iniciais nulas Q e Q'
# t - vetor temporal [s]
 
 # Constantes

  Cs = 0.005;        #[L/cmH2O]
  Rc = 1;           #[cmH2O.s/L]
  Rp = 0.5;          #[cmH2O.s/L]
  Rp = 5;          #[cmH2O.s/L]
  Cl = 0.2;          #[L/cmH2O]
  
 # Variaveis de interesse
  x1 = x(1);                             #Q
  x2 = x(2);                             #dQ 
  dx1 = x(2);                            #dQ
  dx2 = (-(1/Cl)*x1 - ((Cs*Rc/Cl)+Rp+Rc)*x2 + u)/(Rp*Rc*Cs);   
  
  dx = [dx1;dx2];  
endfunction