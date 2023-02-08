%% PROGRAM for Sensor Network Design with Obsv and UKF constraints
%% METODO para Obsv y UKF 
   clear;clc;
%% define cost vector
tic;
   CotaCosto=1500;
   VectorCosto;
%% compute for OCM and UKF constraints
  SimulaCOPOLYM; %para UKF
  ParamUKF; %para UKF y Obs
  %--------------------------------------------------------------------------
% PROCEDIMIENTO DE BUSQUEDA TRANSVERSAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enumeración sin criterio de parada %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;
%--------------------------------------------------------------------------
% PROCEDIMIENTO DE BUSQUEDA EN ARBOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enumeración en arbol con criterio de finalizacion %%%%%%%%%%%%%%%%%%%%%%%
tic;
for k=1:10;
for i1=1:nv;  
           q=zeros(nv,1); q(i1)=1;
           if Co*q<=CotaCosto
                [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Nivel 2
for i1=1:nv; for i2=i1+1:nv;   
           q=zeros(nv,1); q(i1)=1;q(i2)=1;
           if Co*q<=CotaCosto
               [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end;end;  


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nivel 3
for i1=1:nv; for i2=i1+1:nv; for i3=i2+1:nv; 
            q=zeros(nv,1); q(i1)=1;q(i2)=1;q(i3)=1;
           if Co*q<=CotaCosto
               [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end;end;end;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4
for i1=1:nv; for i2=i1+1:nv; for i3=i2+1:nv; for i4=i3+1:nv;   
           q=zeros(nv,1); q(i1)=1;q(i2)=1;q(i3)=1; q(i4)=1;
           if Co*q<=CotaCosto
               [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end;end;end;end;  


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5
for i1=1:nv; for i2=i1+1:nv; for i3=i2+1:nv; for i4=i3+1:nv;for i5=i4+1:nv;   
           q=zeros(nv,1); q(i1)=1;q(i2)=1;q(i3)=1; q(i4)=1; q(i5)=1; 
           if Co*q<=CotaCosto
               [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end;end;end;end;end;  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6
for i1=1:nv; for i2=i1+1:nv; for i3=i2+1:nv; for i4=i3+1:nv;for i5=i4+1:nv; for i6=i5+1:nv;   
           q=zeros(nv,1); q(i1)=1;q(i2)=1;q(i3)=1; q(i4)=1; q(i5)=1; q(i6)=1;
           if Co*q<=CotaCosto
               [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end;end;end;end;end;end;  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%7
for i1=1:nv; for i2=i1+1:nv; for i3=i2+1:nv; for i4=i3+1:nv;for i5=i4+1:nv; for i6=i5+1:nv; for i7=i6+1:nv;    
           q=zeros(nv,1); q(i1)=1;q(i2)=1;q(i3)=1; q(i4)=1; q(i5)=1; q(i6)=1; q(i7)=1;
           if Co*q<=CotaCosto
               [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end;end;end;end;end;end;end;  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%8
for i1=1:nv; for i2=i1+1:nv; for i3=i2+1:nv; for i4=i3+1:nv;for i5=i4+1:nv; for i6=i5+1:nv; for i7=i6+1:nv; for i8=i7+1:nv;    
           q=zeros(nv,1); q(i1)=1;q(i2)=1;q(i3)=1; q(i4)=1; q(i5)=1; q(i6)=1; q(i7)=1; q(i8)=1;
           if Co*q<=CotaCosto
               [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end;end;end;end;end;end;end;end;  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%9
for i1=1:nv; for i2=i1+1:nv; for i3=i2+1:nv; for i4=i3+1:nv;for i5=i4+1:nv; for i6=i5+1:nv; for i7=i6+1:nv; for i8=i7+1:nv; for i9=i8+1:nv;    
           q=zeros(nv,1); q(i1)=1;q(i2)=1;q(i3)=1; q(i4)=1; q(i5)=1; q(i6)=1; q(i7)=1; q(i8)=1; q(i9)=1;
           if Co*q<=CotaCosto
               [S1] = Observabilidad(q,W,nx,ny);
               if S1 == nx;
                  [S7]=PrecisionUKF(nm,nc,Q,R,Xo,Po,Xreal,Y,q,lenghtsimulation,LV,LAMBDA);
                      if S7<=Cota 
                         mQ=q;
                         Cota=S7;
                      end
               end
           end
end;end;end;end;end;end;end;end;end;  
toc;
end;
       