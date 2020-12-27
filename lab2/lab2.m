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

for i = 1:30
    disp(WorkList(i).Box);
    s = ['f(y) = ', num2str(WorkList(i).Estim)];
    disp(s);
end

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
