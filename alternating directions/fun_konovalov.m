function [y] = fun_konovalov(i,j,k)
global hx;
global hy;
global ht;
y=(j*hy)*(ht*(k+(1/2)))^2;
end

