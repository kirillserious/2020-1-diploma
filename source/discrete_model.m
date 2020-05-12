classdef discrete_model < handle
   properties (Access = private)
      mA
      mB
      mM
      mN
      mT
      
      tStart
      tFinish
      startCoord
      
      interval
      delay
      
      mPhi
      mGamma1
      mGamma2
      
      mHPhi
      mHGamma
      mHM
      mHT
      
      startU
      startZ
      n
      
      mPs
      all_us
      all_xs
      all_times
      all_coords
   end
   
   methods
        function obj = discrete_model(A, B, M, N, T, tspan, coord_start, delay, interval)
            obj.mA = A;
            obj.mB = B;
            obj.mM = M;
            obj.mN = N;
            obj.mT = T;
            obj.tStart  = tspan(1);
            obj.tFinish = tspan(2);
            obj.startCoord = coord_start;
            obj.delay = delay;
            obj.interval = interval;
            
            obj.setFirstMatrices();
            obj.setSecondMatrices();
            obj.setStartU();
            
            obj.n = ceil((obj.tFinish - obj.tStart)/obj.interval);
        end
        
        function [Phi, Gamma1, Gamma2] = getFirstMatrices(o)
           Phi = o.mPhi;
           Gamma1 = o.mGamma1;
           Gamma2 = o.mGamma2;
        end
        
        function [HPhi, HGamma, HM, HT] = getSecondMatrices(o)
            HPhi = o.mHPhi;
            HGamma = o.mHGamma;
            HM = o.mHM;
            HT = o.mHT;
        end
        
        function startU = getStartU(o)
            startU = o.startU;
        end
        
        function fit(o)
            o.all_times = [];
            o.all_coords = [];
            o.all_us = [transpose(o.startU)];
            o.all_xs = [transpose(o.startCoord)];
            o.setPs();
            
            prev_u = o.startU;
            curr_u = 0;
            coord  = o.startCoord;
            time = o.tStart;
            for i = 1:o.n
                curr_u = -(o.mN + transpose(o.mHGamma)*o.mPs(:,:,i)*o.mHGamma)\(transpose(o.mHGamma)*o.mPs(:,:,i)*o.mHPhi*[coord; prev_u]);
                o.all_us = [o.all_us; transpose(curr_u)];
                
                
                [times, coords] = o.get_real_system(coord, [time, time + o.delay], prev_u);
                o.all_times = [o.all_times; times];
                o.all_coords = [o.all_coords; coords];
                [times, coords] = o.get_real_system(transpose(coords(end,:)), [time + o.delay, time + o.interval], curr_u);
                o.all_times = [o.all_times; times];
                o.all_coords = [o.all_coords; coords];
                
                coord = o.get_next_discrete_system(coord, curr_u, prev_u);
                o.all_xs = [o.all_xs; transpose(coord)];
       
                
                prev_u = curr_u;
                time = time + (o.tFinish - o.tStart)/o.n;
            end
            
        end
        
        function draw_coords(o)
            myfigure(16), hold on, grid on;
            
            plot(o.all_times, o.all_coords(:,1), 'linewidth', 2);
            plot(o.all_times, o.all_coords(:,2), 'linewidth', 2);
            xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14);
            legend('$x_1$', '$x_2$', 'interpreter', 'latex', 'fontsize', 14);
            xlim([o.tStart, o.tFinish]);
        end
        
        function draw_control(o)
            myfigure(16), hold on, grid on;
            
            times = [o.tStart, o.tStart + o.delay : o.interval : o.tFinish];
            stairs([times, o.tFinish], [o.all_us; o.all_us(end)], 'linewidth', 2);
            xlabel('$t$', 'interpreter', 'latex', 'fontsize', 14);
            ylabel('$u$', 'interpreter', 'latex', 'fontsize', 14);
            xlim([o.tStart, o.tFinish]);
        end
        
        function draw_strange(o)
            myfigure(16), hold on, grid on;
            plot(0:o.n, o.all_us, '-*', 'linewidth', 2);
            xlabel('$k$', 'interpreter', 'latex', 'fontsize', 14);
            ylabel('$u^k$', 'interpreter', 'latex', 'fontsize', 14);
        end
        
        function [times, coords] = get_real_system(o, start_coord, tspan, control)
           fcnSys = @(t, x) o.mA * x + o.mB * control;
           opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
           [times, coords] = ode45(fcnSys, tspan, start_coord, opts);

        end
        
        function coord = get_next_discrete_system(o, start_coord, control, prev_control)
            coord = o.mPhi * start_coord + o.mGamma1 * control + o.mGamma2 * prev_control;
        end
   end
   
   methods (Access = private)
       % Creates Phi, Gamma_1, Gamma_2
       function setFirstMatrices(o)
           o.mPhi = expm(o.mA .* o.interval);
           
           fcnIntegral = @(s) expm(o.mA .* s); 
           o.mGamma2 = integral(fcnIntegral, 0, o.delay, 'ArrayValued', true);
           o.mGamma2 = o.mGamma2 * o.mB;
           o.mGamma1 = integral(fcnIntegral, o.delay, o.interval, 'ArrayValued', true);
           o.mGamma1 = o.mGamma1 * o.mB;
       end
       
       % Creates all tilda matrices
       function setSecondMatrices(o)
          % Create mHPhi
          sz = size(o.mPhi, 2) + size(o.mGamma2, 2);
          o.mHPhi = zeros(sz);
          for i = 1:size(o.mPhi, 1)
              for j = 1:size(o.mPhi, 2)
                  o.mHPhi(i, j) = o.mPhi(i, j);
              end
          end
          for i = 1:size(o.mGamma2, 1)
              for j = 1:size(o.mGamma2, 2)
                  o.mHPhi(i, size(o.mPhi, 2) + j) = o.mGamma2(i,j);
              end
          end
          % Create mHGamma
          o.mHGamma = [o.mGamma1; eye(size(o.mGamma1, 2))];
          % Create mHM
          o.mHM = zeros(sz);
          for i = 1:size(o.mM, 1)
              for j = 1:size(o.mM, 2)
                  o.mHM(i,j) = o.mM(i,j);
              end
          end
          % Create mHT
          o.mHT = zeros(sz);
          for i = 1:size(o.mT, 1)
              for j = 1:size(o.mT, 2)
                  o.mHT(i,j) = o.mT(i,j);
              end
          end
       end
       
       % Creates u_0
       function setStartU(o)
          mP = care(o.mA, o.mB, o.mM, o.mN, [], []);
          o.startU = -o.mN \ transpose(o.mB) * mP * o.startCoord;
          o.startZ = [o.startCoord; o.startU];
       end
       
       % Creates all P^k matrices
       function setPs(o)
            o.mPs = zeros(size(o.mHPhi, 1), size(o.mHPhi, 2), o.n + 1);
            o.mPs(:, :, o.n+1) = o.mHT;
            
            for k = o.n : -1 : 1
                o.mPs(:, :, k) = o.mHM+ ...
                    transpose(o.mHPhi) * o.mPs(:, :, k+1) * o.mHPhi - ...
                    transpose(o.mHPhi) * o.mPs(:, :, k+1) * o.mHGamma * ...
                    inv(o.mN + transpose(o.mHGamma) * o.mPs(:, :, k+1) * o.mHGamma) * ...
                    transpose(o.mHGamma) * o.mPs(:, :, k+1) * o.mHPhi;
            end
       end
   end
end