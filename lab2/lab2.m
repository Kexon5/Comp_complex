path = 'C:\Users\Андрей\Documents\MATLAB';
X = intval([infsup(-5, 5), infsup(-5, 5)]);

% rasstrigin function
[Z, WorkList, diams] = globopt0(X);
iter = 1:1:length(diams);
plot(iter, diams);
hold on;
xlim([0, length(diams)]);
xlabel('Iterations');
ylabel('Diameter of area');
title('Constriction of area');

saveas(gcf, fullfile(path, 'Rastrigin function'), 'png'); 


% Three-hump camel function
[Z, WorkList, diams] = globopt0(X);
solution = 0;

answer = [];
diff = []
for i = 1 : length(WorkList)
    answer(i) = WorkList(i).Estim;
    diff(i) = abs(answer(i) - solution);
end

iter = 1:1:length(answer);
plot(iter, diff);

hold on;
xlim([0, length(answer)]);
xlabel('Iterations');
ylabel('Absolute difference');
title('Convergence of method');
saveas(gcf, fullfile(path, 'Three-hump camel function'), 'png');

min_value = diff(1)
index_min = 0;
for i = 1 : length(WorkList)
   if diff(i) - min_value <= 1e-2
       index_min = i;
       min_value = diff(i)
   end
end
key

x_center = [];
y_center = [];
for i = 101 : length(WorkList)
    x_center(i) = WorkList(i).Box(1).mid;
    y_center(i) = WorkList(i).Box(2).mid;
end

% Trajectory of bar center
x = linspace(-5,5);
y = linspace(-5,5);
[X,Y] = meshgrid(x,y);
Z = 2 .* X .^2 - 1.05 .* X .^ 4 + X .^ 6 / 6 + X .* Y + Y .^ 2;
contour(X,Y,Z, 20)
hold on
plot(x_center, y_center); 
xlabel('Center X');
ylabel('Center Y');
title('Trajectory of bar center');
saveas(gcf, fullfile(path, 'Trajectory center'), 'png'); 
