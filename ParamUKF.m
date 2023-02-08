%% Parámetros de UKF

ALPHA=1; BETA=2; KAPPA=0; 
LV=26; %variables de estado
LAMBDA=((ALPHA^2)*(LV+KAPPA))-LV;

nm(:,1)=LAMBDA/(LV+LAMBDA);
nc(:,1)=(LAMBDA/(LV+LAMBDA))+1-(ALPHA^2)+BETA;

for i=2:2*LV+1
    nm(:,i)=1/(2*(LV+LAMBDA));
    nc(:,i)=1/(2*(LV+LAMBDA));
end

load MW; %para Observabilidad
load Cota; %para UKF
Cota=S7;
