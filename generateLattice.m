function generate_lattice()

 close all;
  
  sn_location = [pi/2-pi/12 pi/2]; 		% sinus node location
  sn_radius = 0.08;						% sinus node radius
  
  mv_location = [pi 0];			   		% mitral valve location
  mv_radius = 85/2/pi/21.2;				% mitral valve radius
  final_rotation = pi - mv_location(1);  
  
  pv_location = [pi/2, pi/3 ;
                 pi/2-pi/5, pi/3 ;
                 pi/2, -pi/3;
                 pi/2-pi/5, -pi/3;];	% PV location
             
  pv_radius = 5/21.2*[1;1;1;1];			% PV radius
  pv_band = 0.1;						% The band around PVs which can burst
  pv_burst_radius = 0.1;				% The radius of the group of cells bursting at the same itme

  initial_defect_std = 0.4;				% initial defect stdev
  initial_defect_center = mv_location*0.65;	% initial defect center


    [L, T] = generate_icosahedron(5);
    
    Lx = L(:,1);
    Ly = L(:,2);
    Lz = L(:,3);
    
  number_neighbors = zeros(size(Lx));
    
  threshold = 0.12; 						% distance within this threshold can be excited
  
  for i=1:size(Lx,1)

      target = [Lx(i), Ly(i), Lz(i)];
      
      Lx_test = Lx;
      Ly_test = Ly;
      Lz_test = Lz;
      
      Lx_test(i) = [];
      Ly_test(i) = [];
      Lz_test(i) = [];
      
      theta = acos(target(1)*Lx_test + target(2)*Ly_test + target(3)*Lz_test );
      number_neighbors(i) = sum(theta <= threshold);

  end

  
          % take out mitral valve
          mv_temp = [sin(mv_location(1))*cos(mv_location(2)), sin(mv_location(1))*sin(mv_location(2)), cos(mv_location(1))];
          theta = acos([Lx Ly Lz]*mv_temp');
          N_mitral_valve = sum( theta < mv_radius );

          Lx = Lx(theta >= mv_radius );
          Ly = Ly(theta >= mv_radius );
          Lz = Lz(theta >= mv_radius );
          number_neighbors = number_neighbors(theta >= mv_radius);


          % take out pulmonary veins
          for i=1:4

              pv_temp = [sin(pv_location(i,1))*cos(pv_location(i,2)), sin(pv_location(i,1))*sin(pv_location(i,2)), cos(pv_location(i,1))];
              theta = acos([Lx Ly Lz]*pv_temp');
              N_pv = sum( theta < pv_radius(i) );
        %       display(N_pv / N);
              Lx = Lx(theta >= pv_radius(i) );
              Ly = Ly(theta >= pv_radius(i) );
              Lz = Lz(theta >= pv_radius(i) );

              number_neighbors = number_neighbors(theta >= pv_radius(i));


          end

  
  
  % generating a neighboring map
  nbhd_map = zeros(size(Lx,1), max(number_neighbors)+1);
  
  for i=1:size(Lx,1)
      
      
      target = [Lx(i), Ly(i), Lz(i)];
      theta = acos( target(1)*Lx + target(2)*Ly + target(3)*Lz );
      ind = find(theta <= threshold);
      
      nbhd_map(i,1) = size(ind,1);
      for j=1:size(ind,1)
          nbhd_map(i,j+1) = ind(j);
      end
      
      
  end
  
  % generating activation points near the pulmonary veins
  pv_activation = [];
  for i=1:4
      
      pv_temp = [sin(pv_location(i,1))*cos(pv_location(i,2)), sin(pv_location(i,1))*sin(pv_location(i,2)), cos(pv_location(i,1))];
      theta = acos([Lx Ly Lz]*pv_temp');
     
      ind = find(theta <= pv_radius(i) + pv_band );
      
      pv_activation = [pv_activation; ind];   
  end
   
  % generating the cellular maps for each activation
  
  
  max_n = 0;
  
  for i=1:size(pv_activation,1)
    
      ind = pv_activation(i);
      target = [Lx(ind), Ly(ind), Lz(ind)];
      theta = acos( target(1)*Lx + target(2)*Ly + target(3)*Lz );
      
      ind = find(theta <= pv_burst_radius);
      
      max_n = max(max_n, size(ind,1));
      
      
  end
  
  activation_map = zeros(size(pv_activation, 1), max_n+1);
  
  for i=1:size(pv_activation,1)
    
      ind = pv_activation(i);
      target = [Lx(ind), Ly(ind), Lz(ind)];
      theta = acos( target(1)*Lx + target(2)*Ly + target(3)*Lz );
      
      ind = find(theta <= pv_burst_radius);
      
      activation_map(i, 1) = size(ind,1);
      activation_map(i, 2:size(ind,1)+1) = ind';

  end
  
   
    
  % determine the sinus location
  [snx, sny, snz] = sph2cart(sn_location(2),pi/2 - sn_location(1), 1.0);
  theta = acos(snx*Lx + sny*Ly + snz*Lz);
  ind = find(theta <= sn_radius);
  sn_list = ind;
    
  
  % now set up an initial defect distribution
  [idcx, idcy, idcz] = sph2cart(initial_defect_center(2),pi/2 - initial_defect_center(1), 1.0);
  theta = acos(idcx*Lx + idcy*Ly + idcz*Lz);
  initial_defect_probability = exp(-theta.^2 / 2 ./ initial_defect_std^2);
  
  
  % SN and PV bursts location cannot have defect
  for i=1:size(initial_defect_probability)
      
      if sum(sn_list == i)>0
          initial_defect_probability(i) = 0;
      end
      
      if sum(pv_activation == i)>0
          initial_defect_probability(i) = 0;
      end
      
  end
  
  
    initial_defect_probability = initial_defect_probability / sum(initial_defect_probability);
   
  
  % final rotation to orient the mitral valve at the south pole
  
  R = [cos(final_rotation) 0 sin(final_rotation);0 1 0;-sin(final_rotation) 0 cos(final_rotation)];
  
  L_temp = [Lx ,Ly, Lz];
  L_temp = L_temp';
  
  L_temp = R*L_temp;
  
  L_temp = L_temp';
  
  Lx = L_temp(:,1);
  Ly = L_temp(:,2);
  Lz = L_temp(:,3);
  
  L = L_temp;
  
  
  % preparing for output;
  fn = fopen('./geometry.txt', 'w');
  % this file saves the dimensionalities / matrices for the simulator
  
  % Number of particles
  fprintf(fn, '%d \n', size(Lx,1));
  
  % Initial defect probabilities
  for i=1:size(Lx,1)
    fprintf(fn, '%.12f \t', initial_defect_probability(i));
  end
  fprintf(fn, '\n');
  
  % Activation points
  fprintf(fn, '%d %d \n', size(activation_map,1), size(activation_map, 2));
  for i=1:size(activation_map,1)
      for j=1:size(activation_map,2)
          fprintf(fn, '%d \t', activation_map(i,j)-1);
      end
      fprintf(fn, '\n');
  end
  
  % Neighborhood maps
  fprintf(fn, '%d \n', size(nbhd_map, 2));
  for i=1:size(nbhd_map,1)
      j = 1;
      fprintf(fn, '%d \t', nbhd_map(i,j));
      
      for j=2:size(nbhd_map,2)
          fprintf(fn, '%d \t', nbhd_map(i,j)-1);
      end
      fprintf(fn, '\n');
  end
  
  %  sinus cells
  fprintf(fn, '%d \n', size(sn_list, 1));
  for i=1:size(sn_list,1)
      fprintf(fn, '%d \t', sn_list(i)-1);
  end
  fprintf(fn, '\n');
  fclose(fn);
  
  fn = fopen('./location.txt', 'w');
  
  % Particle location
  
  for i=1:size(Lx,1)
    fprintf(fn, '%.12f \t', Lx(i));
  end
  fprintf(fn, '\n');
  for i=1:size(Ly,1)
    fprintf(fn, '%.12f \t', Ly(i));
  end
  fprintf(fn, '\n');
  
  for i=1:size(Lz,1)
    fprintf(fn, '%.12f \t', Lz(i));
  end
  fprintf(fn, '\n');
  
  fclose(fn);
  
  % visualisation
  figure(1);
  
  [x,y,z]  = sphere(400);
  surf(0.99*x,0.99*y,0.99*z, 'edgecolor', 'none', 'facecolor', [0 0 0]);
  hold on;
  scatter3(Lx, Ly, Lz, 15, number_neighbors, 'filled');
  hold off;
  view(30, 35);
  colorbar;
  
   
  axis equal; box on;
  xlabel('x');
  ylabel('y');
  zlabel('z');
  
  
  limit = [-1.05 1.05];
  
  set(gca, 'xlim', limit, 'ylim', limit, 'zlim', limit, 'clim', [min(number_neighbors), max(number_neighbors)]);
  title('density fluctuation');
  

  figure(3);
  
  [x,y,z]  = sphere(400);
  surf(0.99*x,0.99*y,0.99*z, 'edgecolor', 'none', 'facecolor', [0 0 0]);
  hold on;
  scatter3(Lx, Ly, Lz, 15, initial_defect_probability, 'filled');
  scatter3(Lx(sn_list), Ly(sn_list), Lz(sn_list), 15, 'filled', 'markerfacecolor', 'r');
  
  hold off;
  view(30, 35);
  colorbar;
  
   
  axis equal; box on;
  xlabel('x');
  ylabel('y');
  zlabel('z');
  
  
  limit = [-1.05 1.05];
  
  set(gca, 'xlim', limit, 'ylim', limit, 'zlim', limit, 'clim', [min(initial_defect_probability), max(initial_defect_probability)]);
  title('initial defect distribution, and SN (red)');
  
    
  figure(4);
  
  
  [x,y,z]  = sphere(400);
  surf(0.99*x,0.99*y,0.99*z, 'edgecolor', 'none', 'facecolor', [0 0 0]);
  hold on;
  
  color_l = zeros(size(Lx,1), 1);
  color_l(sn_list) = 1.0;
  
  scatter3(Lx, Ly, Lz, 15, color_l, 'filled');
  hold off;
  view(30, 35);
  colorbar;
  
   
  axis equal; box on;
  xlabel('x');
  ylabel('y');
  zlabel('z');
  
  
  limit = [-1.05 1.05];
  
  set(gca, 'xlim', limit, 'ylim', limit, 'zlim', limit, 'clim', [min(color_l), max(color_l)]);
  title('Sinus node');
  
  
  
  
 figure(2);
  
  for i=1:1
      cla;
      
      [x,y,z]  = sphere(400);
      surf(0.99*x,0.99*y,0.99*z, 'edgecolor', 'none', 'facecolor', [0 0 0]);
      hold on;
      color_l = zeros(size(Lx,1), 1);
      color_l(pv_activation) = 1.0;
      
      % randomly pick one reaction center
      DD = size(activation_map,1);
      ind = floor(1+DD*rand(1,1));
      numbers = activation_map(ind,1);
      color_l(activation_map(ind, 2:numbers+1)) = 2.0;
      
      scatter3(Lx, Ly, Lz, 15, color_l, 'filled');
      hold off;
      view(30, 35);
      colorbar;


      axis equal; box on;
      xlabel('x');
      ylabel('y');
      zlabel('z');


      limit = [-1.05 1.05];

      set(gca, 'xlim', limit, 'ylim', limit, 'zlim', limit, 'clim', [min(color_l), max(color_l)]);
      title('PV activation candidates');

      view(87,31);
      pause(0.1);
      
  end

  figure(3);
  
    
end

function [vertex, radius] = generate_icosahedron(n)

close all;

    if nargin==0
        n = 0;
    end
        
    % find out the length of the edge;
    x = linspace(0, 2, 10000);
    
    theta1 = acos( (2 - x.^2)/2 ) ;
    theta2 = acos( (2 - 3*x.^2/4)/2);
   
    
    ER = fzero(@edgeLength, 1.097);
  
    radius = ER/2^n;
    

    % now generate a list of points
    L = zeros(0,2);
    
    % first point at north pole
    L(1,:) = [0, 0];
    
    % Next 5 points 
    for i=1:5
        L(1+i,:) = [ER, 2*pi*(i-1)/5];
    end
    
    % Next 5 points 
    for i=1:5
        L(6+i,:) = [pi-ER, 2*pi*(i-1)/5+2*pi/10];
    end
    
    % final point at south pole

    L(12,:) = [pi, 0];
    
    
    [Lx, Ly, Lz] = sph2cart(L(:,2), pi/2 - L(:,1), 1.0);
    
    
    faceMap = [1 2 3;
               1 3 4;
               1 4 5;
               1 5 6;
               1 6 2;
               2 7 3;
               3 8 4;
               4 9 5;
               5 10 6;
               6 11 2;
               2 7 11;
               3 8 7;
               4 8 9;
               5 9 10;
               6 10 11;
               12 7 8;
               12 8 9;
               12 9 10;
               12 10 11;
               12 11 7;];
           
  

   triangle_list = zeros(size(faceMap,1), 9);
   
   for i=1:size(faceMap,1)
       triangle_list(i,1:3) = [Lx(faceMap(i,1)), Ly(faceMap(i,1)), Lz(faceMap(i,1))];
       triangle_list(i,4:6) = [Lx(faceMap(i,2)), Ly(faceMap(i,2)), Lz(faceMap(i,2))];
       triangle_list(i,7:9) = [Lx(faceMap(i,3)), Ly(faceMap(i,3)), Lz(faceMap(i,3))];
   end
   
  
   for iteration = 1:n

       next_generation_triangle_list = zeros(size(triangle_list,1)*4, 9);
       
       for i=1:size(triangle_list,1)
           
            vertex = reshape(triangle_list(i,:), [3 3]);
            vertex = vertex';
            
            vertex1 = vertex(1,:);
            vertex2 = vertex(2,:);
            vertex3 = vertex(3,:);
            
            mid12 = 0.5*(vertex1+vertex2);
            mid23 = 0.5*(vertex2+vertex3);
            mid31 = 0.5*(vertex3+vertex1);
            
            mid12 = mid12/  sqrt(mid12(1)^2 + mid12(2)^2 + mid12(3)^2);
            mid23 = mid23/  sqrt(mid23(1)^2 + mid23(2)^2 + mid23(3)^2);
            mid31 = mid31/  sqrt(mid31(1)^2 + mid31(2)^2 + mid31(3)^2);
            
            
            next_generation_triangle_list( (i-1)*4 + 1, : ) = [vertex1 mid12 mid31];
            next_generation_triangle_list( (i-1)*4 + 2, : ) = [vertex2 mid12 mid23];
            next_generation_triangle_list( (i-1)*4 + 3, : ) = [vertex3 mid23 mid31];
            next_generation_triangle_list( (i-1)*4 + 4, : ) = [mid12 mid23 mid31];
            
       
       end
       
   
       triangle_list = next_generation_triangle_list;
  
   end
   

  vertex = [triangle_list(:,1:3); triangle_list(:,4:6); triangle_list(:,7:9);];
  vertex = unique(vertex, 'rows');
   

    [x,y,z]  = sphere(400);
    surf(0.99*x,0.99*y,0.99*z, 'edgecolor', 'none', 'facecolor', [1 0 0 ]);
    hold on;
    scatter3(vertex(:,1), vertex(:,2), vertex(:,3), 5, 'filled', 'markerfacecolor', 'k');
    hold off;
    display(strcat('total number of particle:', num2str(size(vertex, 1))));
    
  axis equal; box on;
  xlabel('x');
  ylabel('y');
  zlabel('z');
  
  limit = [-1.05 1.05];
  
  set(gca, 'xlim', limit, 'ylim', limit, 'zlim', limit, 'clim', [-10 30]);
  
    
end

function out = edgeLength(theta)


    point1 = [0 , pi/2 - theta, 1];
    point2 = [2*pi/10 , - pi/2 + theta, 1];
    
    [x1,y1,z1] = sph2cart(point1(1), point1(2), point1(3));
    [x2,y2,z2] = sph2cart(point2(1), point2(2), point2(3));
    
    
    
    out = acos(x1*x2 + y1*y2 + z1*z2) - theta;
    
end
