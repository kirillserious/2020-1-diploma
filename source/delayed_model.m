classdef delayed_model < handle
    % DELAYED_MODEL Модель системы с запаздыванием
    
    properties (Access = private)
        matrixA
        matrixB
        matrixM
        matrixN
        matrixT
        
        time_start
        time_finish
        coord_start
        delay
        interval
        
        matrixP
        
        coords
        coord_times
        
        controls
        known_coords
        times
        
        value
    end
    
    methods
        function obj = delayed_model(A, B, M, N, T, tspan, coord_start, delay, interval)
            obj.matrixA = A;
            obj.matrixB = B;
            obj.matrixM = M;
            obj.matrixN = N;
            obj.matrixT = T;
            
            obj.time_start  = tspan(1);
            obj.time_finish = tspan(2);
            obj.coord_start = coord_start;
            obj.delay = delay;
            obj.interval = interval;
            
            obj.matrixP = care(obj.matrixA, ...
                               obj.matrixB, ...
                               obj.matrixM, ...
                               obj.matrixN, ...
                               [], []);
            
            obj.coords        = [];
            obj.coord_times   = [];
            
            obj.times = obj.time_start : obj.interval : obj.time_finish;
            obj.controls     = [];
            obj.known_coords = [];
            
            discrete_delay = ceil(delay/interval);
            current_control = - obj.matrixN \                     ...
                              transpose(obj.matrixB) *            ...
                              obj.matrixP *                       ...
                              obj.coord_start;
            current_coord   = obj.coord_start;
            obj.controls = [obj.controls; transpose(current_control)];
            obj.known_coords = [obj.known_coords; transpose(current_coord)];
            for i = 2:numel(obj.times)    
                [coord_times, coords, current_coord] = obj.get_system(current_coord,  ...
                                               obj.times(i-1),   ...
                                               obj.times(i), ...
                                               current_control);
                
                obj.coord_times = [obj.coord_times; coord_times];
                obj.coords = [obj.coords; coords];
                obj.known_coords = [obj.known_coords; transpose(current_coord)];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                current_control = expm(obj.matrixA .* (obj.times(i) - obj.times(i-1))) * transpose(obj.known_coords(i - 1, :));
                integral_fcn = @(s) expm(A .*(obj.times(i) - s)) * obj.matrixB * obj.controls(i-1);
                current_control = expm(obj.matrixA .* (obj.times(i) - obj.times(i-1))) * current_control + integral(integral_fcn,   ...
                                                                     obj.times(i-1),   ...
                                                                     obj.times(i), ...
                                                                     'ArrayValued', true);
                %if i - discrete_delay < 1
                %    current_control = obj.coord_start;
                %    for j = 1 : i-1
                %        integral_fcn = @(s) expm(A .*(obj.times(j+1) - s)) * obj.matrixB * obj.controls(j);
                %        current_control = expm(obj.matrixA .* (obj.times(j+1) - obj.times(j))) * current_control + integral(integral_fcn,   ...
                %                                                     obj.times(j),   ...
                %                                                     obj.times(j+1), ...
                %                                                     'ArrayValued', true);
                %    end
                %    
                %else 
                %    current_control = expm(obj.matrixA .* (obj.times(i-discrete_delay+1) - obj.times(i - discrete_delay))) * transpose(obj.known_coords(i - discrete_delay, :));
                %    for j = i-discrete_delay : i-1
                %        integral_fcn = @(s) expm(A .*(obj.times(j+1) - s)) * obj.matrixB * obj.controls(j);
                %        current_control = expm(obj.matrixA .* (obj.times(j+1) - obj.times(j))) * current_control + integral(integral_fcn,   ...
                %                                                     obj.times(j),   ...
                %                                                     obj.times(j+1), ...
                %                                                     'ArrayValued', true);
                %    end
                %end
                current_control = - obj.matrixN \ transpose(obj.matrixB) * obj.matrixP * current_control;
                obj.controls = [obj.controls; transpose(current_control)];
            end
            
            first_part = trapz(obj.coord_times,                        ...
                diag(obj.coords * obj.matrixM * transpose(obj.coords)) ...
                );
            second_part = sum(diag(obj.controls * obj.matrixN * transpose(obj.controls))) * obj.interval;
            third_part = obj.coords(end, :) * obj.matrixT * transpose(obj.coords(end, :));
            obj.value = first_part + second_part + third_part;
            
        end
        
        function value = get_value(obj)
            value = obj.value;
        end
        
        function draw_coords(obj)
            myfigure(16), hold on, grid on;
            
            plot(transpose(obj.coord_times), obj.coords, 'linewidth', 2);
            
            xlabel('$t$');
            legend('$x_1$', '$x_2$', 'interpreter', 'latex', 'fontsize', 14);
            
        end
        
        function draw_control(obj)
            myfigure(16), hold on, grid on;
            
            stairs(transpose(obj.times), obj.controls, 'linewidth', 2);
            
            xlabel('$t$');
            ylabel('$u$');
        end
    end
    
    methods (Access = private)
        function [coord_times, coords, coord_to] = get_system(obj, coord_from, time_from, time_to, control)
           system_fcn = @(t, x) obj.matrixA * x + obj.matrixB * control;
           
           [coord_times, coords] = ode45(system_fcn, [time_from time_to], coord_from);
           
           coord_to = transpose(coords(end, :));
        end
        
        
         
        %function control = get_control(obj, time)
        %    if time - obj.delay < time_start + interval
        %        prev_time  = time_start
        %        prev_coord = obj.coord_start
        %        pred_coord = expm(A*(time))
        %    else
        %        
        %    [~, prev_index] = min(abs(obj.times - (time - )))
        %    prev_coord
        %    current_coord = 
        %    control = - obj.matrixN \ transpose(obj.matrixB) * obj.matrixP
        %end
    end
end

