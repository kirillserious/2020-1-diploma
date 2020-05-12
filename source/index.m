%
%    Алгоритмы представлены для постоянных матриц
%
%
%
%% Данные

t_start = 1;
t_finish = 3;
x_start = [4; 100];


A = [-2, -0.02; -1, -10];
B = [2; 0];
M = [1 0; 0 10];
N = [1];
T = [1 0; 0 1];
%%
sm = standart_model(A, B, M, N, T, [t_start, t_finish], x_start);

%%
%sm.get_value()
sm.draw_control()
sm.draw_coords()
%sm.draw_trajectory()
%%
sm.get_value()
%%
dm = delayed_model(A, B, M, N, T, [t_start, t_finish], x_start, 0.5, 0.01);
%%
sm.get_value()
dm.get_value()
%%
dm.draw_control()
dm.draw_coords()
%%
sm.draw_control()
sm.draw_coords()
%% Управление с запаздыванием
eps = 0.1;

f = @(t, x, u) A*x + B*u;
P = care(A, B, M, N, [], []);


for current_t = t_start:eps:t_finish
    
end
%%
times = [0.01 0.025 0.05 0.1 0.25 0.5];
values = [];
for i = 1:numel(times)
  dm = delayed_model(A, B, M, N, T, [t_start, t_finish], x_start, 0.5, times(i));
  values = [values, dm.get_value()];
  'here'
end
%%
myfigure(16), hold on, grid on;
plot(times, ((values - values(1)).^2), '-*', 'linewidth', 2);
xlabel('$\varepsilon$');
ylabel('$(J_\varepsilon - J)^2$');
    
 
%%
dm = delayed_model(A, B, M, N, T, [t_start, t_finish], x_start, 0.5, 0.25);
dm.draw_control();
%%
myfigure(16), hold on, grid on;
plot(0.1:0.1:0.5, [18.5675   34.7678   52.4638   66.0417   80.3556], '-*', 'linewidth', 2);
plot(0.1:0.1:0.5, times, '-*', 'linewidth', 2);
legend('Без оптимизации', 'С оптимизацией');
xlabel('$h$');
ylabel('cpu time, sec.');

%% 
P = care(A(1), B(1), M(1), N(1), [], []);
%%
eps = 0.1;
times = t_start : eps : t_finish;

figure, hold on, grid on;

x_current = x_start;
for i = t_start : eps : t_finish
    i
    %[times, xs] = ode45(@(t) control(t, x_current), [i, i + eps], x_current);
    
    %x_current = xs(end);
    %plot(times, xs, 'r');
end
%%
[times, xs] = ode45(@(t, x) ode_function(t, x, control(t, x)), [0, 7], transpose(x_start));

plot(times, xs);
%plot(xs(:, 1), xs(:, 2));
grid on;
%%
myfigure(9);
%% Dicrete model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval = 0.2;
delay = 0.05;
modelDiscr = discrete_model(A, B, M, N, T, [t_start, t_finish], x_start, delay, interval);
%%
[times, coords] = modelDiscr.get_real_system(x_start, [t_start, t_finish], 5);
myfigure(16), hold on, grid on;
plot(coords(:, 1), coords(:, 2), '-', 'linewidth', 2);

d_coords = transpose(x_start);
coord = x_start;
for i = 1 : ceil((t_finish - t_start)/interval)
    coord = modelDiscr.get_next_discrete_system(coord, 5, 5);
    d_coords = [d_coords; transpose(coord)];
end
plot(d_coords(:, 1), d_coords(:, 2), '*', 'linewidth', 2);
xlabel('$x_1$', 'interpreter', 'latex', 'fontsize', 14);
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize', 14);
legend('Траектория непрерывной системы', 'Орбита дискретной системы', 'interpreter', 'latex', 'fontsize', 14);

%% 
'here'
modelDiscr.fit()


modelDiscr.draw_coords();
%%
modelDiscr.draw_control();
modelDiscr.draw_strange();




function dxdy = ode_function(t, x, control)
    dxdy = A(t) * x + B(t) * control;
end

function result = control(t, x)
    global B N P
    result = - N(t) \ (transpose(B(t)) * P * x);
end

