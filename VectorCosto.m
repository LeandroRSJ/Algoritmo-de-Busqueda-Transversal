%Vector de Costo
 Co=[2 2 20 80 100 100 100 100 100 100 100]*5;
 %Co=sort(Co);
 ub=sum(Co);lb=0;
%Inicializo vac�os las matrices para guardar las redes de sensores que son
%soluci�n y el costo de las mismas
 CF=[];Costo=[];