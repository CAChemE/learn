function z=CSTR_PID_1tank_ODEsystem(t,y)

global k0 Ea alfa nu Cp Tref DH_Tref UA V0 Tr tau Rg R S Qv0 CA_set CA_peturbacion Talim tau_controlador Kc C_entrada

C1=y(1:S); % column vector !!!
T1 =y(S+1);
C1 = C1'; % row vector

CA_manipulada=y(1*S+2);

C1alim=C_entrada;
C1alim(1)=CA_peturbacion+CA_manipulada;


% Reactor 1

T1alim=Talim;
k1=k0.*exp(-Ea./(Rg*T1)); % min-1

for i=1:R
    r1(i)=k1(i)*prod(C1.^nu(i,:)); % mol/(m³·min-1)  (1×2)
end
BM1 =1/tau *(C1alim-C1 )+r1*alfa;% (1×2)*(2×3) mol/(m³·h)

DCp=Cp*alfa';    % J/(mol·K)  (1×3)*(3×2)
DH =DH_Tref +DCp*(T1 -Tref);    % J/mol  (1×2)
%nalim=C1alim*Qv0;
%BE1=(UA *(Tr -T1 )-nalim*Cp'*(T1 -T1alim )-DH *r1'*V0 )/(V0*C1*Cp' ); % ºC/min
BE1=0;

% Corriente manipulada
error=CA_set-C1(1);
dCA_manipulada_dt=Kc*(-BM1(1)+1/tau_controlador*error);


z=[BM1'; BE1; dCA_manipulada_dt];



