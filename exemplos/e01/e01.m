% *************************************************************
%  e01/e01.m   script MATLAB / GNU Octave
% 
%  autor: P. R. Zingano <zingano@gmail.com>
% 
% 
%  Programa de Pós-Graduação em Matemática Aplicada
%  Instituto de Matemática / UFRGS
% 
%  MAP0202 - Métodos Numéricos para Eq. Diferenciais
% 
%  Exemplo 01: Solução da eq. diferencial
% 
%    du/dt = f(t,u(t)) ,
% 
%    f(t,u) = (1 + exp(-t)) * u*u / (1 + u*u*u*u)
% 
%    u(t0) = u0
% 
%    pelo método de Euler com passo uniforme. 
%
% *************************************************************

f = @(t,z) (1+exp(-t))*z^2/(1+z^4);

t0 = 0; tF = 10;
h = 1e-3;
t = t0:h:tF;
nt = length(t);
u = zeros(1,nt);

u(1) = 1;

for n = 1:nt-1
  tn = t(n);
  un = u(n);
  kn = h * f(tn,un);
  u(n+1) = un + kn;
end

%u % imprime a solução no console

figure(10)
plot(t,u,'-r')
xlabel('t')
title('Plot of solution u(t) on [0,10]')