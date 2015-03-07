function z=CSTR_PID_3tanks_ODEsystem(t,y)

global k0 Ea alfa nu Cp Tref DH_Tref UA V0 Tr tau Rg R S Qv_Feed CA_set CA_peturbacion Talim tau_controlador Kc C_entrada

C1=y(1:S);
T1 =y(S+1);

C2=y(S+2:2*S+1);
T2 =y(2*S+2);

C3=y(2*S+3:3*S+2);
T3 =y(3*S+3);

CA_manipulada=y(3*S+4);

C1alim=C_entrada;
C1alim(1)=CA_peturbacion+CA_manipulada;

% transform into row vectors
C1 = C1';
C2 = C2';
C3 = C3';


%

% Reactor 1

T1alim=Talim;
k1=k0.*exp(-Ea./(Rg*T1)); % min-1
for i=1:R
    r1(i)=k1(i)*prod(C1.^nu(i,:)); % mol/(m³·min-1)  (1×2)
end
BM1 =1/tau *(C1alim-C1 )+r1*alfa;% (1×2)*(2×3) mol/(m³·h)

DCp=Cp*alfa';    % J/(mol·K)  (1×3)*(3×2)
DH =DH_Tref +DCp*(T1 -Tref);    % J/mol  (1×2)
%nalim=C1alim*Qv_Feed;
%BE1=(UA *(Tr -T1 )-nalim*Cp'*(T1 -T1alim )-DH *r1'*V0 )/(V0*C1*Cp' ); % ºC/min
BE1=0;


% Reactor 2
C2alim=C1;
T2alim=Talim;
k2=k0.*exp(-Ea./(Rg*T2)); % min-1
for i=1:R
    r2(i)=k2(i)*prod(C2.^nu(i,:)); % mol/(m³·min-1)  (1×2)
end
BM2 =1/tau *(C2alim -C2 )+r2*alfa; % (1×2)*(2×3) mol/(m³·h)
DCp=Cp*alfa';    % J/(mol·K)  (1×3)*(3×2)
DH =DH_Tref +DCp*(T2 -Tref);    % J/mol  (1×2)
%nalim=C2alim*Qv_Feed;
%BE2=(UA *(Tr -T2 )-nalim*Cp'*(T2 -T2alim )-DH *r2'*V0 )/(V0*C2*Cp' ); % ºC/min
BE2=0;


% Reactor 3
C3alim=C2;
T3alim=Talim;
k3=k0.*exp(-Ea./(Rg*T3)); % min-1
for i=1:R
    r3(i)=k3(i)*prod(C3.^nu(i,:)); % mol/(m³·min-1)  (1×2)
end
BM3 =1/tau *(C2 -C3 )+r3*alfa;% (1×2)*(2×3) mol/(m³·h)
DCp=Cp*alfa';    % J/(mol·K)  (1×3)*(3×2)
DH =DH_Tref +DCp*(T3 -Tref);    % J/mol  (1×2)
%nalim=C3alim*Qv_Feed;
%BE3=(UA *(Tr -T3 )-nalim*Cp'*(T3 -T3alim )-DH *r2'*V0 )/(V0*C2*Cp' ); % ºC/min
BE3=0;

% Corriente de alimentacion al reactor 1
%dC1alim_dt=zeros(1,S);
error=CA_set-C3(1);
%dCA_manipulada_dt=Kc*(-BM3(1)+(1/tau_controlador)*error);
dCA_manipulada_dt=Kc*(-BM3(1)+1/tau_controlador*error);
%dC1alim_dt(2)=0;

z=[BM1'; BE1; BM2'; BE2; BM3'; BE3; dCA_manipulada_dt];



