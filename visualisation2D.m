function visualisation2D()
close all;
figure(2);
% here we plot the 'coarse signal'
data =load('coarseEvolution.txt');
T = data(:,1)/400; % convert to physical time
signal = data(:,2);
plot(T, signal);
xlabel('time');
ylabel('signal');
title('coarse signal');

figure(1);
map = load('location.txt');
config = load('detailedEvolution.txt');

T = config(:,1);
config = config(:,2:end);

Lx =map(1,:);
Ly =map(2,:);
Lz =map(3,:);

[phi,el,rho] = cart2sph(Lx, Ly, Lz);
sc = scatter(phi, el, 10, config(1,:), 'filled');
hold on;

% this bit construct the consere-area 2D mapping (Mollweide). 
% Only compute this the first time, the structure will be saved into a temp file named
% tempTheta.mat, any subsequent call of this function will read the saved
% file. If the geometry is changed, the temp file should be removed.


options = optimset(optimset('fsolve'), 'TolFun', 1.0e-8, 'TolX',1.0e-8);
theta = zeros(size(phi));

kk = el(2);
display(testFunc(0,kk));

fn = fopen('tempTheta.mat');

if (fn==-1)
    for i=1:size(phi,2)
        display(i);

        if el>0
            theta(i) = fsolve(@(x) testFunc(x,el(i)), 1.0, options);
        else
            theta(i) = fsolve(@(x) testFunc(x,el(i)), -1.0, options);
        end
    end
    save('tempTheta', 'theta');
else
    fclose(fn);
    load('tempTheta');
end



phi0 = max(phi)*0.5+min(phi)*0.5;

x = 2*sqrt(2)/pi*(phi-phi0).*cos(theta);
y = sqrt(2)*sin(theta);

colormap parula;

xlabel('longitude');
ylabel('latitude');
set(gca, 'xlim', max(x)*[-1.01 1.01], 'ylim', max(y)*[-0.7 1], 'clim', [-10 110]);
set(gca, 'clim', [0 90]);
set(gca, 'xtick', [], 'ytick', []);
axis off;

set(gca, 'color', [0.8 0.8 0.8]);
set_paper_size(0.75*[32, 18]);


for t= 1:4:size(T,1)

  cla;
  delete(sc);
  
  sc = scatter(x, y, 10, config(t,:), 'filled');
  tl = title(strcat(num2str(1/400*T(t)), ' sec')); %added by eugene to display time
  set(tl, 'fontname', 'times new roman', 'fontsize', 30, 'color', 'w');
  set(gcf, 'color', 'k');
  set(gca, 'color', 'k');
  drawnow;
  pause(0.01);
  
end

  

end
function y=testFunc(x, psi)
    y = 2*x+sin(2*x)-pi*sin(psi);
end
function set_paper_size(siz)

    
    
    set(gcf,'paperunit', 'centimeters');
    set(gcf,'papersize', siz);
    set(gcf,'paperposition', [0 0 siz]);
    set(gcf, 'paperpositionmode', 'manual');
    set(gcf, 'color', [0.8 0.8 0.8]);
    set(gcf, 'InvertHardCopy', 'off');
    
    

end
