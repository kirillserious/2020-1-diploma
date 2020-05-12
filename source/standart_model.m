classdef standart_model < handle
    % STANDART_MODEL Стратегия без запаздывания с постоянными матрицами
    
    properties (Access = private)
        matrixA
        matrixB
        matrixM
        matrixN
        matrixT
        
        time_start
        time_finish
        coord_start
        
        matrixP
        
        coords
        controls
        coord_times
        control_times
        
        value
    end
    
    methods
        function obj = standart_model(A, B, M, N, T, tspan, coord_start)
            obj.matrixA = A;
            obj.matrixB = B;
            obj.matrixM = M;
            obj.matrixN = N;
            obj.matrixT = T;
            
            obj.time_start  = tspan(1);
            obj.time_finish = tspan(2);
            obj.coord_start = coord_start;
            
            obj.matrixP = care(obj.matrixA, ...
                               obj.matrixB, ...
                               obj.matrixM, ...
                               obj.matrixN, ...
                               [], []);
            
            obj.coords        = [];
            obj.controls      = [];
            obj.coord_times   = [];
            obj.control_times = [];
            
            options = odeset('RelTol',1e-9);
            [obj.coord_times, obj.coords] = ode45( ...
                @(time, coord) obj.get_system(time, coord, obj.get_control(time, coord)), ...
                [obj.time_start, obj.time_finish], ...
                obj.coord_start,                   ...
                options                            ...
            );
        
            first_part = trapz(obj.coord_times,                    ...
                diag(obj.coords * obj.matrixM * transpose(obj.coords)) ...
                );
            second_part = trapz(obj.control_times,                     ...
                diag(obj.controls * obj.matrixN * transpose(obj.controls)) ...
                );
            third_part = obj.coords(end, :) * obj.matrixT * transpose(obj.coords(end, :));
            obj.value = first_part + second_part + third_part;
        end
        
        function value = get_value(obj)
            value = obj.value;
        end
        
        function draw_control(obj)
            set(0,'DefaultTextInterpreter', 'latex');
            %  Размер шрифта такой же как в ЛаТеХ-документе
            set(0, 'DefaultAxesFontSize', 14);
            set(0, 'DefaultTextFontSize', 14);
            myfigure(16), hold on, grid on;
            plot(obj.control_times, obj.controls, 'linewidth', 2);
            xlabel('$t$');
            ylabel('$u$');
        end
        
        function draw_coords(obj)
            myfigure(16), hold on, grid on;
            plot(obj.coord_times, obj.coords, 'linewidth', 2);
            xlabel('$t$');
            legend('$x_1$', '$x_2$', 'interpreter', 'latex', 'fontsize', 14);
        end
        
        function draw_trajectory(obj)
            figure, hold on, grid on;
            plot(obj.coords(:, 1), obj.coords(:, 2));
        end
    end
    
    methods (Access = private)
        
        function dxdt = get_system(obj, time, coord, control)
            dxdt = obj.matrixA * coord + obj.matrixB * control;
        end
        
        function control = get_control(obj, time, coord)
            control = -obj.matrixN \ transpose(obj.matrixB) * obj.matrixP * coord;
            
            obj.control_times = [obj.control_times; time];
            obj.controls      = [obj.controls; transpose(control)];
        end
    end
end

