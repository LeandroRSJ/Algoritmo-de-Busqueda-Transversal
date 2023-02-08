function dx=ModeloMMA(t,x)
%1-Colocar todos los parametros del sistema. 
%2-Colocar las ecuaciones algebraicas y diferenciales que describen el
%modelo
%3-Las ecuaciones diferenciales deben colocarce de la forma
%dx(i,:)=f(x(i)).
%4-Recordar colocar las ecuaciones en su correspondiente secuencia de
%calculo.
%Ejemplo
%k=0.16;
%A=k*0.5+x(1);
%dx(1,:)=x(1)^2+2;
%dx(2,:)=x(2)*x(1)/A;
%--------------------------------------------------------------------------
% global ti

NoComp=6;
tau=1;
Ai=4.5E14;
Acaa=4.209E11;
Acbb=1.6E9;
Adaa=0;
Adbb=0;
Apaa=3.207E6;
Apab=1.233E5; 
Apba=2.10E8; 
Apbb=6.308E6;
Axaa=32.08;
Axab=1.234;
Axas=86.6;
Axat=2085;
Axba=5.257E4;
Axbb=1577;
Axbs=1514;
Axbt=4.163E5;
Aza=2.2;
Azb=1.13E5;
Ei=1.25E5;
Ecaa=2.69E4;
Ecbb=4E3;
Edaa=0;
Edbb=0;
Epaa=2.42E4;
Epab=2.42E4;
Epba=1.8E4;
Epbb=2.42E4;%1.8E4
Exaa=2.42E4;
Exab=2.42E4;
Exas=2.42E4;
Exat=2.42E4;
Exba=1.8E4;
Exbb=1.8E4;
Exbs=1.8E4;
Exbt=2.42E4;%1.8E4;
Eza=0;
Ezb=0;
Hpaa=-54E3;
Hpba=-54E3;
Hpab=-86E3;
Hpbb=-86E3;
ror=8.79E2;
Cr=2.01;
Ur=6E-2;
Vr=1;
Sr=4.6;
M1=100.12;
M2=86.09;
M3=164.21;
M4=78.11;
M5=44.05;
M6=168.11;
Rg=8.3144; 
Tj=336.15;
Cf7=353;

Ff1=4.994E-5;
Ff2=2.904E-4;
Ff3=3.045E-7;
Ff4=1.28E-4;
Ff5=1.703E-5;
Ff6=0;
landaf1=0;
landaf2=0;
psip0f=0;
psip1f=0;
psip2f=0;

RP=0.95;

%Ecuaciones cineticas.
ki= Ai*exp(-Ei/(Rg*x(7)));
kcaa=Acaa*exp(-Ecaa/(Rg*x(7)));
kcbb=Acbb*exp(-Ecbb/(Rg*x(7)));
kdaa=Adaa*exp(-Edaa/(Rg*x(7)));
kdbb=Adbb*exp(-Edbb/(Rg*x(7)));
kpaa=Apaa*exp(-Epaa/(Rg*x(7)));
kpab=Apab*exp(-Epab/(Rg*x(7)));
kpba=Apba*exp(-Epba/(Rg*x(7)));
kpbb=Apbb*exp(-Epbb/(Rg*x(7)));
kxaa=Axaa*exp(-Exaa/(Rg*x(7)));
kxab=Axab*exp(-Exab/(Rg*x(7)));
kxas=Axas*exp(-Exas/(Rg*x(7)));
kxat=Axat*exp(-Exat/(Rg*x(7)));
kxba=Axba*exp(-Exba/(Rg*x(7)));
kxbb=Axbb*exp(-Exbb/(Rg*x(7)));
kxbs=Axbs*exp(-Exbs/(Rg*x(7)));
kxbt=Axbt*exp(-Exbt/(Rg*x(7)));
kza=Aza*exp(-Eza/(Rg*x(7)));
kzb=Azb*exp(-Ezb/(Rg*x(7)));
kcab=((kcaa*kcbb+1E-20)^(1/2));
kdab=((kdaa*kdbb+1E-20)^(1/2));    

% Balance de Masa

% Velocidad instanteanea de polimerizacion

% Calculo concentracion de Radicales
beta=((kpab+kxab)*x(2))/((kpba+kxba)*x(1));
I1=kcaa+kdaa+2*beta*(kcab+kdab)+beta^2*(kcbb+kdbb);
I2=x(6)*(kza+beta*kzb);
I3=-2*ki*x(3)*tau;
Carad=(-I2+sqrt(I2^2-4*I1*I3))/(2*I1);
Cbrad=beta*Carad;
R1=((kpaa+kxaa)*Carad+(kpba+kxba)*Cbrad)*x(1);
R2=((kpbb+kxbb)*Cbrad+(kpab+kxab)*Carad)*x(2);
R3=ki*x(3);
R4=(kxas*Carad+kxbs*Cbrad)*x(4);
R5=(kxat*Carad+kxbt*Cbrad)*x(5);
R6=(kza*Carad+kzb*Cbrad)*x(6);
Gpi2=(R1*M1+R2*M2)*Vr;

Qf=(Ff1*M1+Ff2*M2+M3*Ff3+M4*Ff4+M5*Ff5+M6*Ff6)/ror;
titar=Vr/Qf;

Cf1=Ff1/Qf; 
Cf2=Ff2/Qf; 
Cf3=Ff3/Qf; 
Cf4=Ff4/Qf; 
Cf5=Ff5/Qf; 
Cf6=Ff6/Qf; 

% Cf1=5.15; 
% Cf2=2.97; 
% Cf3=4.6E-3; 
% Cf4=1.31; 
% Cf5=1.17; 
% Cf6=0; 

% Conv1=(Cf1-x(1))/(Cf1+1E-20); 
% Conv2=(Cf2-x(2))/(Cf2+1E-20); 
% Conv3=(Cf3-x(3))/(Cf3+1E-20); 
% Conv4=(Cf4-x(4))/(Cf4+1E-20); 
% Conv5=(Cf5-x(5))/(Cf5+1E-20); 
% Conv6=(Cf6-x(6))/(Cf6+1E-20); 

Qfint1=Ff1*M1/ror;
Qfint2=Ff2*M2/ror;
Qfint3=Ff3*M3/ror;
Qfint4=Ff4*M4/ror;
Qfint5=Ff5*M5/ror;
Qfint6=Ff6*M6/ror;

dx(1,:)=((Cf1-x(1))/titar)-R1;
dx(2,:)=((Cf2-x(2))/titar)-R2;
dx(3,:)=((Cf3-x(3))/titar)-R3;
dx(4,:)=((Cf4-x(4))/titar)-R4;
dx(5,:)=((Cf5-x(5))/titar)-R5;
dx(6,:)=((Cf6-x(6))/titar)-R6;
    
% Ecuaciones Adicionales de los balances de Masa y energia

% Termino de generacion

% Balace de Energia
Tj=x(25);
dx(7,:)=((Cf7-x(7))/titar)+(((-Hpaa)*kpaa*x(1)*Carad+(-Hpba)*kpba*x(1)*Cbrad)/(ror*Cr))+(((-Hpab)*kpab*x(2)*Carad+(-Hpbb)*kpbb*x(2)*Cbrad)/(ror*Cr))-(Ur*Sr*(x(7)-Tj))/(Vr*ror*Cr);

%Composicion del copolimero muerto.
dx(8,:)=((landaf1-x(8))/titar)+R1;
dx(9,:)=((landaf2-x(9))/titar)+R2;
  
Yap=x(8)/(x(8)+x(9)+1E-10);

% Momentos.

% % Ecuaciones adicionales para los momentos

L1=kxas*x(4)+kxaa*x(1)+kxab*x(2)+kxat*x(5)+kza*x(6)+kdaa*Carad+kdab*Carad;
L2=kxbs*x(4)+kxbb*x(2)+kxba*x(1)+kxbt*x(5)+kzb*x(6)+kdbb*Cbrad+kdab*Carad;

alpha1=kpaa*x(1)/((kcaa+kdaa)*Carad+(kcab+kdab)*Cbrad+(kpaa+kxaa)*x(1)+(kpab+kxab)*x(2)+kxat*x(5)+kxas*x(4)+kza*x(6));
alpha2=kpbb*x(2)/((kcbb+kdbb)*Cbrad+(kcab+kdab)*Carad+(kpbb+kxbb)*x(2)+(kpba+kxba)*x(1)+kxbt*x(5)+kxbs*x(4)+kzb*x(6));

c1=((2*ki*x(3)*tau+x(4)*(kxas*Carad+kxbs*Cbrad))/(kpaa*(x(1)+x(2))))+((x(5)*(kxat*Carad+kxbt*Cbrad))/(kpaa*(x(1)+x(2))))+((kxaa*Carad+kxba*Cbrad)/kpaa);
c4=((2*ki*x(3)*tau+x(4)*(kxas*Carad+kxbs*Cbrad))/(kpbb*(x(1)+x(2))))+((x(5)*(kxat*Carad+kxbt*Cbrad))/(kpbb*(x(1)+x(2))))+((kxbb*Cbrad+kxab*Carad)/kpbb);

r1=kpaa/kpab;
r2=kpbb/kpba;

gamma=kpba/kpab;

c2=c4*r2*gamma;
c3=c1*r2/gamma;

g=1/(r1*r2);

V1=c2*g-c1;
V2=c3*g-c4;

B1=1-(alpha1+alpha2)+alpha1*alpha2*(1-g);
B2=(M1+M2)*(1-g)*alpha1*alpha2-alpha1*M1-alpha2*M2;
B3=alpha1*c1+alpha1*alpha2*V1;
B4=alpha2*c4+alpha1*alpha2*V2;

%Momento de orden 0 de A
psiarad0=B3/B1;

%Momento de orden 1 de A
psiarad1=((alpha1*alpha2*V1*(M1+M2)+alpha1*c1*M1)/B1)-(B2*B3/(B1^2));

%Momento de orden 2 de A
psiarad2=((-B3*((alpha1+alpha2)*M1*M2+B2*(M1+M2-1)))/(B1^2))+((alpha1*alpha2*V1*(M1+M2-1)*(M1+M2))/(B1))+((alpha1*c1*M1*(M1-1))/B1)+(2*B2^2*B3/(B1^3))+psiarad1-(2*B2*(alpha1*alpha2*V1*(M1+M2)+alpha1*c1*M1))/(B1^2);

%Momento de orden 0 de B
psibrad0=B4/B1;

%Momento de orden 1 de B
psibrad1=((alpha1*alpha2*V2*(M2+M1)+alpha2*c4*M2)/B1)-(B2*B4/(B1^2));

%Momento de orden 2 de B
psibrad2=((-B4*((alpha1+alpha2)*M1*M2+B2*(M1+M2-1)))/(B1^2))+((alpha1*alpha2*V2*(M1+M2-1)*(M1+M2))/(B1))+((alpha2*c4*M2*(M2-1))/B1)+(2*B2^2*B4/(B1^3))+psibrad1-(2*B2*(alpha1*alpha2*V2*(M1+M2)+alpha2*c4*M2))/(B1^2);

%Momento de orden 0 de P
dx(10,:)=(((psip0f-x(10))/titar)+1/2*kcaa*psiarad0^2+kcab*psiarad0*psibrad0+1/2*kcbb*psibrad0^2+L1*psiarad0+L2*psibrad0);

%Momento de orden 1 de P
dx(11,:)=(((psip1f-x(11))/titar)+kcaa*psiarad0*psiarad1+kcab*(psiarad0*psibrad1+psibrad0*psiarad1)+kcbb*psibrad0*psibrad1+L1*psiarad1+L2*psibrad1);

%Momento de orden 2 de P
dx(12,:)=(((psip2f-x(12)/titar))+kcaa*(psiarad1^2+psiarad0*psiarad2)+kcab*(2*psiarad1*psibrad1+psibrad2*psiarad0+psiarad2*psibrad0)+kcbb*(psibrad1^2+psibrad0*psibrad2)+L1*psiarad2+L2*psibrad2);

% SEPARADOR

Frout1=x(1)*Qf;
Frout2=x(2)*Qf;
Frout3=x(3)*Qf;
Frout4=x(4)*Qf;
Frout5=x(5)*Qf;
Frout6=x(6)*Qf;

Qsint1=Frout1*M1/ror;
Qsint2=Frout2*M2/ror;
Qsint3=Frout3*M3/ror;
Qsint4=Frout4*M4/ror;
Qsint5=Frout5*M5/ror;
Qsint6=Frout6*M6/ror;

Qp=Qf-(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6);

Csf1=Frout1/(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);
Csf2=Frout2/(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);
Csf3=Frout3/(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);
Csf4=Frout4/(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);
Csf5=Frout5/(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);
Csf6=Frout6/(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);

dx(13,:)=((Csf1-x(13))/titar);
dx(14,:)=((Csf2-x(14))/titar);
dx(15,:)=((Csf3-x(15))/titar);
dx(16,:)=((Csf4-x(16))/titar);
dx(17,:)=((Csf5-x(17))/titar);
dx(18,:)=((Csf6-x(18))/titar);

% PURGA

Fsout1=x(13)*(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);
Fsout2=x(14)*(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);
Fsout3=0;
Fsout4=x(16)*(Qsint1+Qsint2+Qsint3+Qsint4+Qsint5+Qsint6+Qp);
Fsout5=0;
Fsout6=0;

Qsout1=RP*Qsint1;
Qsout2=RP*Qsint2;
Qsout3=RP*Qsint3;
Qsout4=RP*Qsint4;
Qsout5=RP*Qsint5;
Qsout6=RP*Qsint6;

%Hold Tank

Chf1=Fsout1/(Qsout1+Qsout2+Qsout4);
Chf2=Fsout2/(Qsout1+Qsout2+Qsout4);
Chf3=Fsout3/(Qsout1+Qsout2+Qsout4);
Chf4=Fsout4/(Qsout1+Qsout2+Qsout4);
Chf5=Fsout5/(Qsout1+Qsout2+Qsout4);
Chf6=Fsout6/(Qsout1+Qsout2+Qsout4);

Fhout1=x(19)*(Fsout1/(x(13)+1E-10));
Fhout2=x(20)*(Fsout2/(x(14)+1E-10));
Fhout3=x(21)*(Fsout3/(x(15)+1E-10));
Fhout4=x(22)*(Fsout4/(x(16)+1E-10));
Fhout5=x(23)*(Fsout5/(x(17)+1E-10));
Fhout6=x(24)*(Fsout6/(x(18)+1E-10));

Qh1=Fhout1*M1/ror;
Qh2=Fhout2*M2/ror;
Qh3=Fhout3*M3/ror;
Qh4=Fhout4*M4/ror;
Qh5=Fhout5*M5/ror;
Qh6=Fhout6*M6/ror;

F2h1=x(19)*(Qh1+Qh2+Qh4);
F2h2=x(20)*(Qh1+Qh2+Qh4);
F2h3=x(21)*(Qh1+Qh2+Qh4);
F2h4=x(22)*(Qh1+Qh2+Qh4);
F2h5=x(23)*(Qh1+Qh2+Qh4);
F2h6=x(24)*(Qh1+Qh2+Qh4);

dx(19,:)=((Chf1-x(19))/titar);
dx(20,:)=((Chf2-x(20))/titar);
dx(21,:)=((Chf3-x(21))/titar);
dx(22,:)=((Chf4-x(22))/titar);
dx(23,:)=((Chf5-x(23))/titar);
dx(24,:)=((Chf6-x(24))/titar);

Trs=353.1;
dx(25,:)=0.0001*(Trs-x(7)); %Tj
dx(26,:)=Gpi2-x(26); %Gpi-x(26)

% MEZCLADOR

Q6m1=(F2h1*M1)/ror;
Q6m2=(F2h2*M2)/ror;
Q6m3=(F2h3*M3)/ror;
Q6m4=(F2h4*M4)/ror;
Q6m5=(F2h5*M5)/ror;
Q6m6=(F2h6*M6)/ror;

Q8m1=Qfint1+Q6m1;
Q8m2=Qfint2+Q6m2;
Q8m3=Qfint3+Q6m3;
Q8m4=Qfint4+Q6m4;
Q8m5=Qfint5+Q6m5;
Q8m6=Qfint6+Q6m6;
end