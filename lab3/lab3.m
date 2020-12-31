addpath('IntLinIncR2_UTF8');
addpath('IntLinIncR3_UTF8');

Aconst = [19 11;10 15;14 12]
A=[infsup(17,21) infsup(9,13);infsup(8,12) infsup(13,17);infsup(10,18) infsup(10,14)]
x=[0.3;0.2]
bconst=Aconst*x
b=[infsup(6.4,9.4);infsup(4.5,7.5);infsup(4.6,8.6)]

[V, P1, P2, P3, P4]=EqnTolR2(inf(A),sup(A),inf(b),sup(b))
Cminim=cond(inf(A));
MatrixM=[0 0;0 0;0 0]
b1=inf(A);

for i1=inf(A(1,1)):(sup(A(1,1))-inf(A(1,1))):sup(A(1,1))
for i2=inf(A(1,2)):(sup(A(1,2))-inf(A(1,2))):sup(A(1,2))
for j1=inf(A(2,1)):(sup(A(2,1))-inf(A(2,1))):sup(A(2,1))
for j2=inf(A(2,2)):(sup(A(2,2))-inf(A(2,2))):sup(A(2,2))
for k1=inf(A(3,1)):(sup(A(3,1))-inf(A(3,1))):sup(A(3,1))
for k2=inf(A(3,2)):(sup(A(3,2))-inf(A(3,2))):sup(A(3,2))
if cond([i1 i2;j1 j2;k1 k2])<=b1
Cminim=cond([i1 i2;j1 j2;k1 k2]);
MatrixM=[i1 i2;j1 j2;k1 k2];
end
end
end
end
end
end
end
disp('Cmin=')
disp(Cminim)
disp('Matrix=')
disp(MatrixM)
disp('__________________')

maxTol = 0.41754343;
argmaxTol = [0.28947346, 0.21228088]';

ive=sqrt(2)*Cminim*maxTol*(norm(argmaxTol)/norm(mid(b)))
Xsolv=[argmaxTol-ive,argmaxTol+ive]

rectangle('Position',[Xsolv(1,1) Xsolv(2,1) Xsolv(1,2)-Xsolv(1,1) Xsolv(2,2)-Xsolv(2,1) ] ,'EdgeColor','b')
rectangle('Position',[argmaxTol(1) argmaxTol(2) 0.003 0.003 ],'EdgeColor','r')
text(argmaxTol(1)+0.02,argmaxTol(2),'argmaxTol','FontSize',8)

Title_str='ISLAU Solution'
title(Title_str)
[m, n] =size(Aconst)
xlabel('\it x_1')
ylabel('\it x_2')
title_str_name=strcat(Title_str,' ',num2str(m),' x ',num2str(n))
figure_name_out=strcat(title_str_name,'.png')
print('-dpng', '-r300', figure_name_out), pwd



disp('_____________________________Second_Part_____________________________')
Aconst2=[19 11;10 15;14 12]'
A2=[infsup(17,21) infsup(9,13);infsup(8,12) infsup(13,17);infsup(10,18) infsup(10,14)]'
x2=[0.3;0.2;0.1]
bconst2=Aconst2*x2
b2=[infsup(7.6,10.6);infsup(5.5,9.5)]


Cminim2=cond(inf(A2));
MatrixM2=[0 0;0 0;0 0]';
b12=inf(A2);

for i1=inf(A2(1,1)):(sup(A2(1,1))-inf(A2(1,1))):sup(A2(1,1))
for i2=inf(A2(1,2)):(sup(A2(1,2))-inf(A2(1,2))):sup(A2(1,2))
for i3=inf(A2(1,3)):(sup(A2(1,3))-inf(A2(1,3))):sup(A2(1,3))
for j1=inf(A2(2,1)):(sup(A2(2,1))-inf(A2(2,1))):sup(A2(2,1))
for j2=inf(A2(2,2)):(sup(A2(2,2))-inf(A2(2,2))):sup(A2(2,2))
for j3=inf(A2(2,3)):(sup(A2(2,3))- inf(A2(2,3))):sup(A2(2,3))
if cond([i1 i2 i3;j1 j2 j3])<=b12
Cminim2=cond([i1 i2 i3;j1 j2 j3]);
MatrixM2=[i1 i2 i3;j1 j2 j3];
end
end
end
end
end
end
end
disp('Cmin2=')
disp(Cminim2)
disp('Matrix2=')
disp(MatrixM2)
disp('__________________')

maxTol2 = 0.3639954;
argmaxTol2 = [3.79999850e-01 1.88001053e-01 -5.81838840e-07]';

ive2=sqrt(3)*Cminim2*maxTol2*(norm(argmaxTol2)/norm(mid(b2)))
Xsolv2=[argmaxTol2-ive2,argmaxTol2+ive2]

[Z]=EqnTolR3(inf(A2),sup(A2),inf(b2),sup(b2),1,1)

Title_str='ISLAU Solution'
title(Title_str)
[m, n] =size(Aconst2)
xlabel('\it x_1')
ylabel('\it x_2')
zlabel('\it x_3')

title_str_name=strcat(Title_str,' ',num2str(m),' x ',num2str(n))
figure_name_out=strcat(title_str_name,'.png')
print('-dpng', '-r300', figure_name_out), pwd


Xsolv=[argmaxTol2-ive2/2,argmaxTol2+ive2/2];
A_p = [infsup(17,21) infsup(8,12); infsup(9,13) infsup(13, 17)];
[V,P1,P2,P3,P4]=EqnTolR2(inf(A_p),sup(A_p),inf(b2),sup(b2));

rectangle('Position', [
Xsolv(1,1) Xsolv(2,1) Xsolv(1,2)-Xsolv(1,1) Xsolv(2,2)-Xsolv(2,1)]);
rectangle('Position',[argmaxTol2(1) argmaxTol2(2) 0.001 0.001 ],'EdgeColor','r');

title_str='Projection';
title(title_str);
xlabel('x_1');
ylabel('x_2');

title_str_name=strcat(title_str);
figure_name_out=strcat(title_str_name,'.png');
print('-dpng', '-r300', figure_name_out), pwd