function [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
%% UFK Function
%% PARAMETROS
M1=100.12;    % peso molecular del metil metacrilato
M2=86.09;     % peso molecular del vinil acrílico
%--------------------------------------------------------------------------
Xest_UKF(:,1)=Xo;
Pest_UKF(:,:,1)=Po;

%- compute Y and R using the sensor network proposed

Y=Y(find(q),:);
R=R(find(q),find(q));
EP1=[];

for k=2:lenghtsimulation

    for i=1:LV;
    if Pest_UKF(i,i,k-1)==0;
       Pest_UKF(i,i,k-1)=1E-15;
    end
    if Xest_UKF(i,k-1)<1E-10;
       Xest_UKF(i,k-1)=0;
    end
    end
   
    % Calculate error covariance matrix square root and the sigma-points
    sqrt_Pest(:,:,k)=chol(Pest_UKF(:,:,k-1),'lower');  
    sigma_point(:,:,k)=[Xest_UKF(:,k-1),(Xest_UKF(:,k-1)*ones(1,LV)+(sqrt(LV+LAMBDA)*sqrt_Pest(:,:,k))),(Xest_UKF(:,k-1)*ones(1,LV)-(sqrt(LV+LAMBDA)*sqrt_Pest(:,:,k)))];
                       
    for i=1:2*LV+1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %options=odeset('RelTol',1e-7,'AbsTol',1e-10);
    [timei,Xi]=ode15s('ModeloMMA',((((k-2)*300)):30:(((k-1)*300))),sigma_point(:,i,k)');%,options
    if size (Xi,1)==11
    CHI_pred(:,i,k)=Xi(11,:)';
    nm_CHI_pred(:,i,k)=nm(:,i)*CHI_pred(:,i,k);
     
     x1=CHI_pred(7,i,k);x2=CHI_pred(25,i,k);x3=CHI_pred(26,i,k); 
     x4=CHI_pred(12,i,k)/CHI_pred(11,i,k);
     x5=CHI_pred(1,i,k);x6=CHI_pred(2,i,k);
     x7=CHI_pred(13,i,k);x8=CHI_pred(14,i,k);
     x9=CHI_pred(19,i,k);x10=CHI_pred(20,i,k);
     x11=CHI_pred(11,i,k)/(CHI_pred(11,i,k)+M1*CHI_pred(1,i,k)+M2*CHI_pred(2,i,k));
     YY=[x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11];
     
     PSI_pred(:,i,k)=YY(find(q),:); 
     nm_PSI_pred(:,i,k)=nm(:,i)*PSI_pred(:,i,k);
    else
     CHI_pred(:,i,k)=Xi(1,:)';
     nm_CHI_pred(:,i,k)=nm(:,i)*CHI_pred(:,i,k);
     
     PSI_pred(:,i,k)=YY(find(q),:);
     nm_PSI_pred(:,i,k)=nm(:,i)*PSI_pred(:,i,k);
    end
    end

    Xpred_UKF(:,k)=sum(nm_CHI_pred(:,:,k),2);
    Ypred_UKF(:,k)=sum(nm_PSI_pred(:,:,k),2);

    for i=1:2*LV+1
     nc_Xpred(:,:,i)=nc(:,i)*((CHI_pred(:,i,k)-Xpred_UKF(:,k))*((CHI_pred(:,i,k)-Xpred_UKF(:,k)))');
     nc_Ypred(:,:,i)=nc(:,i)*((PSI_pred(:,i,k)-Ypred_UKF(:,k))*((PSI_pred(:,i,k)-Ypred_UKF(:,k)))');
     nc_XYpred(:,:,i)=nc(:,i)*((CHI_pred(:,i,k)-Xpred_UKF(:,k))*((PSI_pred(:,i,k)-Ypred_UKF(:,k)))');
    end
    
Pxpred_UKF(:,:,k)=Q+(sum(nc_Xpred,3));
Py_UKF(:,:,k)=R+(sum(nc_Ypred,3));
Pxy_UKF(:,:,k)=(sum(nc_XYpred,3));

% 3. ACTUALIZACIÓN

% Ganancia del filtro K 

K_UKF(:,:,k)=Pxy_UKF(:,:,k)*(Py_UKF(:,:,k)^-1);

% El valor de las variables y su respectiva covarianza se actualizan con la
% ganancia del filtro K y la medición obtenida para cada instante k.

Xest_UKF(:,k)=Xpred_UKF(:,k)+ K_UKF(:,:,k)*(Y(:,k)-Ypred_UKF(:,k));

Pest_UKF(:,:,k)=Pxpred_UKF(:,:,k)-(K_UKF(:,:,k)*Py_UKF(:,:,k)*K_UKF(:,:,k)'); 

EP1=[EP1 (diag(Q+(sum(nc_Xpred,3))))];

end
Variables_interes=[8 9];
for i=1:(lenghtsimulation-1)
    EP1(Variables_interes,i)=EP1(Variables_interes,i)./EP1(Variables_interes,1);
end
EP1=sum(EP1,2);EP1=sum(EP1)-EP1(12);EP1=EP1/(length(Variables_interes)*(lenghtsimulation-1));S7=EP1;
end