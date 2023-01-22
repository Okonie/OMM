  % задача №2 ОММ вариант 4
 %схема переменных направлений
 
 % ut = (uxx + uyy)+yt^2
 % x = [0; 1], y = [0;2]
 % t = [0, 10] (это время амплитуда упадет в ~e раз)

 %граничные условия 
 %ux(0)=u(1)=0
 %uy(0)=uy(2)=0
 %начальное условие u(t=0) = (x^2-1)*cos(pi*y)
 %аналитическое решение u(x, y, t)
 
 a=1;
global hx;
global hy;
global ht;

%Число шагов по x;
Nx = 200;
%Число шагов по y;
Ny = 200;
%Число шагов по t;
Nt = 200;
%Отрезок по x;
x0=0;
x1=1;
%Отрезок по y;
y0=0;
y1=2;
%Отрезок по t;
t0=0;
t1=0.01;
%Шаг h по x,а q по t;
hx=(x1-x0)/(Nx);
hy=(y1-y0)/(Ny);
ht=(t1-t0)/(Nt);
 x_grid = linspace(x0, x1, Nx); %для построения сетки
 y_grid = linspace(y0, y1, Ny);
 t_grid = linspace(t0, t1, Nt);
 x_mesh = zeros(Ny, Nx);
 for p = 1:Ny
 x_mesh(p, :) = x_grid;
 end

 y_mesh = zeros(Ny, Nx);
 for p = 1:Nx
 y_mesh(:, p) = y_grid'; % ' - транспонирование
 end

  %сетки для вывода u(y, t) при x = pi/4 = const
 y_mesh_x = zeros(Ny, Nt);
 for p = 1:Nt
 y_mesh_x(:, p) = y_grid';
 end

 t_mesh_x = zeros(Ny, Nt);
 for p = 1:Ny
 t_mesh_x(p, :) = t_grid;
 end
 
 x_mesh_y = zeros(Nx, Nt);
 for p = 1:Nt
 x_mesh_y(:, p) = x_grid';
 end

 t_mesh_y = zeros(Nx, Nt);
 for p = 1:Nx
 t_mesh_y(p, :) = t_grid;
 end
 
u_x = zeros(Ny, Nt);
u_y = zeros(Nx, Nt);

u_an=zeros(Ny,Nx);
u_an=exp(-5*pi^2*t1).*sin(2*pi*x_mesh).*sin(pi*y_mesh);


f2 = figure(2);
 surf(x_mesh, y_mesh, u_an, 'LineStyle', 'none');
 title('Аналитическое решение');
 xlabel('x');
 ylabel('y');
 zlabel('u');
 colorbar;
 
pointx=zeros(Nx,Ny); %матрица для численного решения
for i=1:1:Nx
    for j=1:1:Ny
      pointx(i,j)=(sin(2*pi*i*hx))*sin(pi*hy*j);  
    end
end

for k=0:1:Nt-1
     %сохраняем решение
 u_x(:, k+1) = pointx(Nx/4, :)';
  u_y(:, k+1) = pointx(:, Ny/4);
f1 = figure(1);
mesh(x_mesh, y_mesh, pointx'); 
xlabel('x');
ylabel('y');
zlabel('u');
ax=ht/(2*hx^2);
bx=ht/(2*hx^2);
cx=1+ht/(hx^2);

%граничные условия для x 
%for j=1:1:Ny 
%pointx(1,j)=0; 
%pointx(Nx,j)=0; 
%end 

%%%%%%%%%%%%%%%%%%%%%%%
%   Прямая  прогонка  %
%%%%%%%%%%%%%%%%%%%%%%%

%Формируем правый столбец
for i=2:1:Nx-1
    for j=2:1:Ny-1
fx(i,j)= pointx(i,j)+ht*(pointx(i,j+1)-2*pointx(i,j)+pointx(i,j-1))/(2*hy^2)*a;
    end
end
%Формируем левую часть
for j=2:1:Ny-1
    alphax(2,j)=0;
    betax(2,j)=0;
    for i=2:1:Nx-1
        alphax(i+1,j)=bx/(cx-alphax(i,j)*ax);
        betax(i+1,j)=(ax*betax(i,j)+fx(i,j))/(cx-alphax(i,j)*ax);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%  Обратная прогонка  %
%%%%%%%%%%%%%%%%%%%%%%%
for j=2:1:Ny-1
    pointx(Nx,j)=0; % граничное условие дирихле
    for i=Nx-1:-1:1
        pointx(i,j)=alphax(i+1,j)*pointx(i+1,j)+betax(i+1,j);
    end
end


%-----------------------------------------------------------

ay=ht/(2*hy^2);
by=ht/(2*hy^2);    
cy=1+ht/(hy^2);

%%%%%%%%%%%%%%%%%%%%%%%
%   Прямая прогонка   %
%%%%%%%%%%%%%%%%%%%%%%%

%Формируем правый столбец
for j=2:1:Ny-1
    for i=2:1:Nx-1
        fy(i,j)=pointx(i,j)+ht*(pointx(i+1,j)-2*pointx(i,j)+pointx(i-1,j))/(2*hx^2)*a;
    end
end

%Формируем левую часть
for i=2:1:Nx-1
    alphay(i,2)=0; 
    betay(i,2)=0;
    for j=2:1:Ny-1
        alphay(i,j+1)=by/(cy-alphay(i,j)*ay);
        betay(i,j+1)=(ay*betay(i,j)+fy(i,j))/(cy-alphay(i,j)*ay);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%  Обратная прогонка  %
%%%%%%%%%%%%%%%%%%%%%%%

for i=2:1:Nx-1
    pointx(i,Ny)=0;%betay(i,Ny)/(1-alphay(i,Ny));
    for j=Ny-1:-1:1
        pointx(i,j)=alphay(i,j+1)*pointx(i,j+1)+betay(i,j+1);
    end
end
for i=1:1:Nx-1  %из граничных условий неймана
pointx(i,1)=0; 
pointx(i,Ny)=0; 
end

%-----------------------------------------------------------



% 
%    filename = 'test.gif';
%     %surf(pointy);
%     %mesh(pointx);
%     colorbar;
%    % axis([0 100 0 100 -1 2]);
%     f = getframe;
%     [im,map] = rgb2ind(f.cdata,256);
%     im(:,:,1,j) = rgb2ind(f.cdata,map);
%     imwrite(im,map,filename,'DelayTime',0,'Loopcount',0);
%     
%     % СНАЧАЛА ПРОГНАТЬ fun_karpov.m, ПОТОМ ПРОГНАТЬ main_karpov.m ииии,
%     % наконец, заккоментировать 120-ую строку, и расскомментировать 30.
%     %Та - даааам, все работает.  
end

 f5 = figure(5);
 surf(x_mesh, y_mesh, abs((pointx' - u_an)), 'LineStyle', 'none');
 title('Ошибка');
 xlabel('x');
 ylabel('y');
 zlabel('u');
 colorbar;

 f6 = figure(6);
 surf(t_mesh_x, y_mesh_x, u_x, 'LineStyle', 'none');
 title('Зависимость u(y,t) при x = 1/2');
 xlabel('t');
 ylabel('y');
 zlabel('u');
 colorbar;
% axis([0, tN, 0, 2, -1, 1]);

 f7 = figure(7);
 surf(t_mesh_y, x_mesh_y, u_y, 'LineStyle', 'none');
 title('Зависимость u(y,t) при y = 1');
 xlabel('t');
 ylabel('y');
 zlabel('u');
 colorbar;
% axis([0, tN, 0, 2, -1, 1]);
















