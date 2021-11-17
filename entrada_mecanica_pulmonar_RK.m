function u = entrada_mecanica_pulmonar_RK(f,phi, Pfis)
  u_aux = Pfis*sin(phi);           #Pw
  u_der = 2*pi*f*Pfis*cos(phi);    #dPw
  u = [u_aux;u_der];           #[mmHg]
endfunction
