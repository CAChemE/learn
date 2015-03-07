clear all
close all
clc
% Problema 7, curso 2007-08, RCTA en régimen no estacionario
% 


%global k0 Ea alfa nu Cp Tref DH_Tref Talim UA V0 Tr tau nalim Calim Rg R S
global k0 Ea alfa nu Cp Tref DH_Tref UA V0 Tr tau Rg R S Qv_Feed CA_set CA_peturbacion Talim tau_controlador Kc C_entrada
%%
%{       
Se quiere estudiar el funcionamiento en régimen no estacionario de un RCTA
que opera a volumen constante y con un tiempo de residencia de 2 s. La
 corriente de alimentación al reactor contiene únicamente reactivo A con
 una concentración de 0.8 de kmol/m³. El sistema reaccionante es el siguiente:
A-->B  k1=40 min^{-1 }, orden de reacción =1.2
B-->C  k2=20 min^{-1 } orden de reacción = 0.5
En un determinado instante se produce una peturbación en la corriente de
alimentación, incrementándose la concentración de A en dicha corriente
en 0.2 de kmol/m³. Se pide:
a) (30%) Dibujar la evolución de la concentración de los componentes con
 el tiempo a partir del estado incial estacionario donde únicamente 
existe componente A en el reactor a una concentración de 0.4 kmol/m³.

b) (45 %) Se  desea mantener la concentración del reactivo A en el reactor
 en su valor estacionario, 0.4 kmol/m³ (concentración de A de consigna). Con
tal fin se manipula la concentración de A en la corriente de alimentación.
De esta forma, la concentración del reactivo A en la corriente de alimentación deja
de ser constante y se convierte en una función del tiempo determinada 
por la siguiente expresión:
error=CA_consigna-CA
CA_manipulad(t)=0.8+Kc*[(error)+1/tau_controlador(int(error*dt)].
donde Kc=0.5 y  tau_controlador=0.1
Dibujar la evoulción de la concentración de los componentes con el tiempo
si se activa el sistema de control. 

c) (25 %) Repetir el apartado interior sustituyendo el reactor anterior,
por un sistema compuesto por 3 reactores idénticos conectados en serie. Los
parámetros del sistema de control para este sistema son Kc=30 y
 1/tau_controlador=15.



NOTA: obsérvese que dCA_manipulad(t)_dt=+Kc*[-dCA/dt+1/tau_controlador*(error)].

%
---Reactor (operating conditions) 
Type of reactor = ideal CSTR at unsteady state
T ~= f(t)
V ~= f(t)
Size of the reactor: residence time := tau = 2 min

---Reaction system
alpha_{ij}A_j(liquid) = 0;
number of components = 3
number of reactions = 2
partial reaction order matrix is known
k_{i}, Ea_{i}

---Feed stream data
Concentrations for all the components in the feed stream are known


---Initial values of the state variables
Concentration for all the components at t=0 := C0

%}
%%
% DATOS DEL REACTOR
% RCTA con rebosamiento contínuo (V=cte)
%UA=1.3e6; %J/(min·K)
%Tr=10+273; % K
%V0=40;  % m³

Qv_Feed=10;     %m³·min^(-1)  
tau=2;  % min
%tau=V0/Qv0;  %

T0=273.15; % Temperatura de todo el sistema [K]
Rg=8.31441; %kJ·(kmol·k)^(-1)


%--- PARAMETROS DEL SIST. REACCIONANTE--------------------------------------------------------
% Reacción en fase líquida
% Matriz con los coeficientes estequiométricos

%        ];
alfa= [-1 1 0   
       0 -1 1];
nu=   [ 1.2  0 0  
        0 0.5 0];

    [R S]=size(alfa);


% Datos cinéticos
%k0=39.1*ones(1,R);
k0= [40 20] ;  % min^{-1 }<
Ea=9900*ones(1,R);

%k=k0*exp(-Ea/(Rg*T0));

% Datos termodinámicos

% Pesos molec ulares
PM=45*ones(1,S);    % kg/kmol



% Datos termodinámicos
Cp=1e3*ones(1,S); % kJ/(kmol·K)


Tref=293; % K
DH_Tref=-1e5*ones(1,R); % kJ/kmol




%--- FIN PARAMETROS DEL SIST. REACCIONANTE--------------------------------------------------------



% Intervalo simulación
t0=0;
tf=20; % min

%% 1 único reactor
%
C_entrada=zeros(1,S);
C_entrada(1)=[0.8 ]; % vector concentraciones corriente de entrada al sistema [kmol·m^(-3)]

C01=zeros(1,S);
C01(1)=0.4; % concentracion de A a t=0 en el reactor 1 [kmol·m^3]
T01=T0;

CA_peturbacion=0.2;
CA0_manipulada=C_entrada(1);

y0=[C01 T01 CA0_manipulada];

Talim=T0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CA_set=0.4;    %Concentración de consigna del componente A en el reactor 3
Kc=5;  %ganacia controlador
%Kc=0;  %ganacia controlador
tau_controlador=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t y]=ode45(@CSTR_PID_1tank_ODEsystem,[t0 tf], y0,4000);

C1=y(:,1:S);
T1=y(:,S+1);

CA_manipulada=y(:,1*S+2);

C1Aalim=CA_peturbacion+CA_manipulada;


figure(1)
subplot(211)
plot (t,[C1 ] )
title(['Estado estacionario inicial  CA0 = ' num2str(y0(1)) ' mol/m³, CB_0 = ' num2str(y0(2)) ' mol/m³' ] )
legend('C1_A', 'C1_B', 'C1_C')
xlabel('t (min)');
ylabel('C (kmol/m³)');
grid

figure(1)
subplot(212)
plot (t,[ CA_manipulada.*Qv_Feed ] )
title(['control action' ] )
legend('F_{A,manipulated}')
xlabel('t (min)');
ylabel('F_{A}^{Feed} (kmol/min)');
grid


%% 3 Reactores en serie
C_entrada=zeros(1,S);
C_entrada(1)=[0.8 ]; % vector concentraciones corriente de entrada al sistema [kmol·m^(-3)]

C01=zeros(1,S);
C01(1)=0.4; % concentracion de A a t=0 en el reactor 1 [kmol·m^3]
T01=T0;

C02=zeros(1,S);
C02(1)=[0.2 ]; % concentracion a t=0 en el reactor 2 [kmol·m^3]
T02=T0;

C03=zeros(1,S);
C03(1)=[0.1 ]; % concentracion a t=0 en el reactor 2 [kmol·m^3]
T03=T0;


%C1alim0=Cent+Cent_peturbacion;
CA_peturbacion=0.2;
CA0_manipulada=C_entrada(1);

y0=[C01 T01 C02 T02 C03 T03 CA0_manipulada];


Talim=T0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CA_set=0.1;    %Concentración de consigna del componente A en el reactor 3
Kc=30;  %ganacia controlador
%Kc=0;
tau_controlador=15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t y]=ode45(@CSTR_PID_3tanks_ODEsystem,[t0 tf], y0,1000);

C1=y(:,1:S);
T1=y(:,S+1);

C2=y(:,S+2:2*S+1);
T2=y(:,2*S+2);

C3=y(:,2*S+3:3*S+2);
T3=y(:,3*S+3);

CA_manipulada=y(:,3*S+4);

C1Aalim=CA_peturbacion+CA_manipulada;


figure(2)
subplot(211)
plot (t,[C1(:,1) C2(:,1) C3(:,1)  ] )
title(['3 RCTA en serie:  C1A_0 = ' num2str(C1(1)) ' mol/m³, C2A_0 = ' num2str(C2(1)) ' mol/m³ , C3A_0 = ' num2str(C3(1)) ' mol/m³' ] )
legend('C1_A', 'C2_A', 'C3_A')
xlabel('t (min)');
ylabel('C (kmol/m³)');
grid


figure(2)
subplot(212)
plot (t,[ CA_manipulada.*Qv_Feed ] )
title(['control action' ] )
legend('F_{A,manipulated}')
xlabel('t (min)');
ylabel('F_{A}^{Feed} (kmol/min)');
grid

