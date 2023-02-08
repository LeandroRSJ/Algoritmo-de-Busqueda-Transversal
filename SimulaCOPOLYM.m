% Program used for CSTR simulation 
%% compute right-hand sides of balance equations

ts=20;
paso=300;
%XINIT=[0.01 0.01 0.01 0.01 0.01 0 100 0.01 0.01 0.01 10 100 0.01 0.01 0.01 0.01 0.01 0 0.01 0.01 0 0.01 0 0];
XINIT=[0.2535972 5.8386115 0.0020120 2.7575710 0.3663659 0 353.4259628 0.8223480 0.4179869 0.0054007 119.7059113 4177280.5766291 0.2535972 5.8386115 0.0020120 2.7575710 0.3663659 0 0.3156237 7.2666578 0 3.4320359 0 0 336.15 0.0055];
%XINIT=1.2*XINIT;
%% Simulation
%% PARAMETROS
M1=100.12;    % peso molecular del metil metacrilato
M2=86.09;     % peso molecular del vinil acrílico
%%
n=19;
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
[timer,Xreal]=ode15s('ModeloMMA',(0:paso:ts*paso),XINIT,options);
Xo=XINIT;

% COVARIANCE Q R

% Q covarianza del ruido del proceso
% R covarianza del ruido de las mediciones

% se asume un mínimo porcentaje de los valores de las variables de estado
% estacionario que sea igual al máximo desvío(tres desvíos estándar):

% Porcentaje*Xee1 = 3*(desvío estándar)

Porcentje_q = 8e-3;
Porcentje_q1 = 5e-4;

% Desvío estándar del ruido del proceso

q = (((Porcentje_q*1e-2).*XINIT(1,:))./3); %resto de las variables de estado
q(1,7)= (((Porcentje_q1*1e-2).*XINIT(1,7))./3); %Tr
q(1,25)= (((Porcentje_q1*1e-2).*XINIT(1,25))./3); %Tj

% Covarianza del ruido del proceso

Q=diag((q).^2);
Q=1.2*Q;
% se asume un porcentaje de los valores de las variables medidas de acuerdo
% con la presición de los dispositivos de medición, que sea igual al máximo 
% desvío(tres desvíos estándar):

% Porcentaje*Xee1 = 3*(desvío estándar)

Porcentje_C = 5;     %Conv. Total y Concentración
Porcentje_T = 1;   %Temperatura
Porcentje_MW = 8;     %Mw (Peso molecular del polímero en peso)PROBAR ESTE!!!

% Desvío estándar del ruido de la medición

r1 =((Porcentje_T*1e-2)*(XINIT(1,7))/3);
r2 =((Porcentje_T*1e-2)*(XINIT(1,25))/3);
r3 =((Porcentje_T*1e-2)*(XINIT(1,26))/3);
r4 =((Porcentje_MW*1e-2)*(XINIT(1,12)/XINIT(1,11))/3);
r5 =((Porcentje_C*1e-2)*(XINIT(1,1))/3);
r6 =((Porcentje_C*1e-2)*(XINIT(1,2))/3);
r7 =((Porcentje_C*1e-2)*(XINIT(1,13))/3);
r8 =((Porcentje_C*1e-2)*(XINIT(1,14))/3);
r9 =((Porcentje_C*1e-2)*(XINIT(1,19))/3);
r10 =((Porcentje_C*1e-2)*(XINIT(1,20))/3);
r11 =((Porcentje_C*1e-2)*(XINIT(1,11)/(XINIT(1,11)+XINIT(1,1)*M1+XINIT(1,2)*M2))/3);


% Covarianza del ruido de la medición

R= diag([(r1)^2;
         (r2)^2; 
         (r3)^2;
         (r4)^2;
         (r5)^2;
         (r6)^2; 
         (r7)^2;
         (r8)^2;
         (r9)^2;
         (r10)^2;
         (r11)^2]);
  
% Po Covarianza del error de estimación del punto incial
% se asume un mínimo porcentaje de los valores de las variables de estado
% estacionario que sea igual al máximo desvío(tres desvíos estándar):

% Este porcentaje debe ser chico porque el punto inical es asumida 
% una buena estimación.

% Porcentaje*Xee1 = 3*(desvío estándar)

Porcentje_po = 1e-3;

% Desvío estándar del error de estimación inicial

po=((Porcentje_po*1e-2)*XINIT(1,:)/3);

% Covarianza del error de estimación inicial

Po= diag((po).^2);    

for i=1:length(XINIT(1,:))
    if Q(i,i)==0;
       Q(i,i)=1e-15;   
    end
    
    if Po(i,i)==0;
       Po(i,i)=1e-15;   
    end
end

%% Trayectoria real

nx=26; %número de variables del vector estado
ny=11;  %número de variables medidas
nw=length(Q);
nv=length(R);

aleat_w(:,:)=randn(nw,ts); %ruidos aleatorios para el proceso
aleat_v(:,:)=randn(nv,ts); %ruidos aleatorios para las mediciones

Xreal(1,:)=XINIT(1,:);

for k=2:ts

% Se evalúa el modelo de las mediciones para cada Xr(k,:)

hxr(:,k)=[Xreal(k,7); %Tr
          Xreal(k,25); %Tj
          Xreal(k,26); %Gpi 
          Xreal(k,12)/(Xreal(k,11));%Mw (peso molecular del pol en peso)
          Xreal(k,1); %CMMA
          Xreal(k,2); %CVA
          Xreal(k,13); %CMMA
          Xreal(k,14); %CVA
          Xreal(k,19); %CMMA
          Xreal(k,20); %CVA
          Xreal(k,11)/(Xreal(k,11)+Xreal(k,1)*M1+Xreal(k,2)*M2)];  %Conv Total.
              
% Se genera las mediciones con ruidos aleatorios

measur_noise(:,k)=(aleat_v(:,k)).*sqrt(diag(R));

Y(:,k) = hxr(:,k) + measur_noise(:,k);

end


minS7=1;
mQ=[];
lenghtsimulation=ts;












