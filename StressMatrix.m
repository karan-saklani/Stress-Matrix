solver3D();
function [data,sigmaPN,sigmaPS] = solver3D()
sigmaPN = NaN;
sigmaPS = NaN;
data = readtable('StressMatrix.xlsx');
delete StressMatrixSolution.xlsx
T = table2array(data(1:end,1:3));
[I1,I2,I3] = invariants(T);
ax1 = table2array(data(1:end,4));
ax2 = table2array(data(1:end,5));
ax3 = table2array(data(1:end,6));
N = table2array(data(1:end,7));
sum1 = sum(isnan(ax1 + ax2 + ax3));
sum2 = sum(isnan(N));
if sum2 == 0
    l = N(1);
    m = N(2);
    n = N(3);
    [Tf,sigmaPN,sigmaPS] = on_plane(T,l,m,n);
    data.Pxyz = Tf(1:end,1);
end
if sum1 == 0
    Tf = on_axis(T,ax1,ax2,ax3);
    data.T1_ax = Tf(1:end,1);
    data.T2_ax = Tf(1:end,2);
    data.T3_ax = Tf(1:end,3);
end
data.I = [I1;I2;I3];
C = [1,-I1,I2,-I3];
s = roots(C);
data.principal_stress = s;
writetable(data,'StressMatrixSolution.xlsx');
end

function [Tf,sigmaPN,sigmaPS] = on_plane(T,l,m,n)
N = [l;m;n];
sigma = T*[l;m;n];
sigmaPx = sigma(1);
sigmaPy = sigma(2);
sigmaPz = sigma(3);
Tf = [sigmaPx;sigmaPy;sigmaPz];
sigmaPN = dot(Tf,N);
sigmaPS = sqrt(sigmaPx^2 + sigmaPy^2 + sigmaPz^2 - sigmaPN^2);
end

function Tf = on_axis(T,ax1,ax2,ax3)
sigmaX = T*ax1;
sigmaY = T*ax2;
sigmaZ = T*ax3;
XX = dot(sigmaX,ax1);
XY = dot(sigmaX,ax2);
XZ = dot(sigmaX,ax3);
YY = dot(sigmaY,ax2);
YZ = dot(sigmaY,ax3);
ZZ = dot(sigmaZ,ax3);
Tf = [XX,XY,XZ;...
    XY,YY,YZ;...
    XZ,YZ,ZZ];
end

function [I1,I2,I3] = invariants(T)
XX = T(1,1);
YY = T(2,2);
ZZ = T(3,3);
XY = T(1,2);
XZ = T(1,3);
YZ = T(2,3);
I1 = XX + YY + ZZ;
I2 = det([XX,XY;XY,YY]) + det([XX,XZ;XZ,ZZ]) + det([YY,YZ;YZ,ZZ]);
I3 = det(T);
end

function s = principal_stress(T)
[I1,I2,I3] = invariants(T);
C = [1,-I1,I2,-I3];
s = roots(C);
end

function Td = deviator_stress(T)
sigm = mean([T(1,1),T(2,2),T(3,3)]);
Td = T - sigm*eye(3,3)
end

function Tf = on_2Daxis(T,theta)
xx = T(1,1);
xy = T(1,2);
yy = T(2,2);
XX = (xx + yy)/2 + (cos(2*theta)*(xx - yy)/2) + (sin(2*theta)*xy);
YY = (xx + yy)/2 - (cos(2*theta)*(xx - yy)/2) - (sin(2*theta)*xy);
XY = -(sin(2*theta)*(xx - yy)/2) + (cos(2*theta)*xy);
Tf = [XX,XY;XY,YY];
end

function [I1,I2,I3] = invariants2D(T)
XX = T(1,1);
YY = T(2,2);
XY = T(1,2);
I1 = XX + YY;
I2 = det([XX,XY;XY,YY]);
I3 = 0;
end

function [c,r] = Mohr(xx,xy,yy)
c = (xx + yy)/2;
r = (sqrt((xx - yy)^2 + 4*xy^2))/4;
end

function [XX,YY,XY] = Transform_2D(xx,yy,xy,theta)
XX = xx*cos(theta)^2 + yy*sin(theta)^2 + 2*xy*sin(theta)*cos(theta);
YY = xx*cos(theta)^2 + yy*sin(theta)^2 - 2*xy*sin(theta)*cos(theta);
XY = -(xx - yy)*sin(theta)*cos(theta) + xy*(cos(theta)^2 - sin(theta)^2);
end
