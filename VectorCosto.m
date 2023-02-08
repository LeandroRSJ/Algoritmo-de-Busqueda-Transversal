%Vector de Costo
 Co=[2 2 20 80 100 100 100 100 100 100 100]*5;
 %Co=sort(Co);
 ub=sum(Co);lb=0;
%Inicializo vacíos las matrices para guardar las redes de sensores que son
%solución y el costo de las mismas
 CF=[];Costo=[];