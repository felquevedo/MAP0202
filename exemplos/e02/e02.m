% *************************************************************
%  e02/e02.m   script MATLAB / GNU Octave
% 
%  autor: P. R. Zingano <zingano@gmail.com>
% 
% 
%  Programa de Pós-Graduação em Matemática Aplicada
%  Instituto de Matemática / UFRGS
% 
%  MAP0202 - Métodos Numéricos para Eq. Diferenciais
% 
%  Exemplo 02: Solução do P.V.I.
% 
%    du/dt = f(t,u(t),v(t)) ,
% 
%    dv/dt = g(t,u(t),v(t)) ,
% 
%    f(t,u,v) = (1 + exp(-t)) * u * v / (1 + u*v*v) ,
% 
%    g(t,u,v) =  v / u,
% 
%    u(t0) = 1 , v(t0) = 1/2
% 
%    pelo método de Euler com passo uniforme.
%
% *************************************************************

f = @(t,z) [ (1+exp(-t))*z(1)*z(2)/(1+z(1)*z(2)^2); z(2)/z(1) ];

t0 = 0; tF = 10;
h = 1e-3;
t = t0:h:tF;
nt = length(t);
w = zeros(2,nt);

w(:,1) = [1;0.5];

for n = 1:nt-1
  tn = t(n);
  wn = w(:,n);
  kn = h * f(tn,wn);
  w(:,n+1) = wn + kn;
end

figure(10), plot(t,w(1,:),'-r'), xlabel('t')
title('Plot of u(t) on [0,10]')

figure(20), plot(t,w(2,:),'-r'), xlabel('t')
title('Plot of v(t) on [0,10]')
