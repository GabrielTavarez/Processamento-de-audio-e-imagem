function [w, erro, y_chapeu] = lms(x, d, M, mu)

% w = matriz dos coeficientes do filtro transversal
%   cada coluna k representa os valores de w na k-�sima intera��o
% erro = sinal de erro (d-y)
%
% x = sequ�ncia de observa��o
% d = sequ�ncia desejada
% M = n�mero de coeficientes do filtro adaptativo
% mu = passo de adapta��o

N = length(x);
w = zeros(M,N);
erro = zeros(1,N);
phi = zeros(M,1);
y_chapeu = zeros(1,N);

for n=1:N-1
    phi = [x(n); phi(1:M-1)];
    y_chapeu(n) = w(:,n)'*phi; 
    erro(n) = d(n) - y_chapeu(n);
    w(:,n+1) = w(:,n) + mu*erro(n)*phi;
end
end