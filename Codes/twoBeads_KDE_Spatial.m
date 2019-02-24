clear all
close all

dimension = 2;
temperature = 500;
datalength = 1200;
count = 1;

load(['../data/LangevinData/BeadsPos_',num2str(dimension),'_',num2str(temperature),'_',num2str(datalength),'_',num2str(count),'.mat'])

data = x.';
minx = -50;
miny = -20;
maxx = 50;
maxy = 20;

rd = kde2d(data);
db1 = rd(1);
db2 = rd(2);

tempvx = zeros(1,N);
tempvy = zeros(1,N);
for k = 2:N-1
    tempvx(k-1) = (x(1,k+1)-x(1,k-1))/(2*dt);
    tempvy(k-1) = (x(2,k+1)-x(2,k-1))/(2*dt);
end

tempv = sqrt( tempvx.^2 + tempvy.^2 );

vpos = x(:,2:R-1);

rn = numel(tempv');
d = size(vpos',2);
hy = median(abs(tempv'-median(tempv')))/0.6745*(4/(d+2)/rn)^(1/(d+4));
hx = median(abs(vpos'-repmat(median(vpos'),rn,1)))/0.6745*(4/(d+2)/rn)^(1/(d+4));
hx = sqrt(hy*hx);
vb1 = hx(1); vb2 = hx(2);

res1 = 19;
res = res1 + 1;


dx1 = (maxx-minx)/res1;
dy1 = (maxy-miny)/res1;

xpos = minx:dx1:maxx;
ypos = miny:dy1:maxy;

%%  Density
density = zeros(res,res);
normd = 1/sqrt(2*db1^2*pi) * 1/sqrt(2*db2^2*pi);
for i = 1:res
    for j = 1:res
        deltax = (xpos(i) - vpos(1,:)).^2;
        deltay = (ypos(j) - vpos(2,:)).^2;
        density(i,j) = sum( normd .* exp( -deltax / (2*db1^2) ) .* exp( -deltay / (2*db2^2) ) ) / N;
    end
end

density = density / (sum(sum(density))*dx1*dy1);

% figure()
% imagesc([minx maxx], [miny maxy], density');
% set(gca,'ydir','normal');
% colorbar;
% xlabel('X1','FontSize',16);
% ylabel('X2','FontSize',16);
% title('Density','FontSize',16);
% hold on
% scatter(x(1,:),x(2,:),'.w')
% hold on

%%  Current New Method
currentx = zeros(res2,res2);
currenty = zeros(res2,res2);
gridyx = zeros(res2,res2);
gridyy = zeros(res2,res2);
for i = 1:res2
    for j = 1:res2
        for k = 2:Nnew-1
            currentx(i,j) = currentx(i,j) + 1/(Nnew-2) * ( (x(1,k+1)-x(1,k-1))/(2*dt) * 1/sqrt(2*a1^2*pi) * 1/sqrt(2*a2^2*pi) * exp( -(minx+(i-1)*dx2 - x(1,k))^2 / (2*a1^2) ) * exp( -(miny+(j-1)*dy2 - x(2,k))^2 / (2*a2^2) ) );
            currenty(i,j) = currenty(i,j) + 1/(Nnew-2) * ( (x(2,k+1)-x(2,k-1))/(2*dt) * 1/sqrt(2*a1^2*pi) * 1/sqrt(2*a2^2*pi) * exp( -(minx+(i-1)*dx2 - x(1,k))^2 / (2*a1^2) ) * exp( -(miny+(j-1)*dy2 - x(2,k))^2 / (2*a2^2) ) );
        end
        gridyx(i,j) = minx+(i-1)*dx2;
        gridyy(i,j) = miny+(j-1)*dy2;
    end
end
%
% figure()
% imagesc([minx maxx], [miny maxy], density');
% set(gca,'ydir','normal');
% colorbar;
% hold on
% quiver(gridyx,gridyy,currentx,currenty,'w','AutoScaleFactor',1.2,'LineWidth',1);
% xlabel('X1','FontSize',16);
% ylabel('X2','FontSize',16);
% title('Density & Current','FontSize',16);

tic

vxK = zeros(res,res);
vyK = zeros(res,res);
for i = 1:res
    for j = 1:res
        
        deltax = (xpos(i) - vpos(1,:)).^2;
        deltay = (ypos(j) - vpos(2,:)).^2;
        
        % Using Epanechnikov kernel
        temp1pos = find(deltax < vb1);
        temp2pos = find(deltay(temp1pos) < vb2);
        pos = temp1pos(temp2pos);
        
        if  isempty(temp2pos) == 0
            tempK = (1 - deltax(pos) ./ vb1) .* (1 - deltay(pos) ./ vb2);
            
            vxK(i,j) = sum( tempvx(pos) .* tempK ) / sum(tempK);
            vyK(i,j) = sum( tempvy(pos) .* tempK ) / sum(tempK);
            
        end
        
    end
end
toc

Diffusionmatrix = inv([T1/2,0;0,T2/2]);

V_sim = zeros(1,2,res,res);
for i = 1:res
    for j = 1:res
        V_sim(:,:,i,j) = V_sim(:,:,i,j) + [vxK(i,j),vyK(i,j)];
    end
end


%% Calculate the Dissipation
tempDissipation = zeros(res,res);
for i = 1:res
    for j = 1:res
        tempDissipation(i,j) = tempDissipation(i,j) + ( V_sim(:,:,i,j) * Diffusionmatrix * V_sim(:,:,i,j).' ) * density(i,j);
    end
end
