function [portion] = fliptria(hemi,i,ico,featype,F0,id,N,path)
f0 = F0-1;
all_id = 1:size(F0,1);
validID = setdiff(all_id,id);
valid_tria = F0(validID,:);
system(strcat('/',path,'/SphericalRemesh -i /',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.white.aff.s5.vtk -o /',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.white.res.vtk -s /',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk -r /',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.movedv.unfold.real.vtk'));

[v,f] = read_vtk(strcat('/',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.white.res.vtk'));
area1 = zeros(size(valid_tria,1),1);
for l = 1:size(valid_tria,1)
    area1(l,1) = norm(cross(v(valid_tria(l,1),:)-v(valid_tria(l,2),:),v(valid_tria(l,1),:)-v(valid_tria(l,3),:)))/2;
end
tri_area = sum(area1);

total_tria = F0(:,:);
area2 = zeros(size(total_tria,1),1);
for l = 1:size(total_tria,1)
    area2(l,1) = norm(cross(v(total_tria(l,1),:)-v(total_tria(l,2),:),v(total_tria(l,1),:)-v(total_tria(l,3),:)))/2;
end
total_area = sum(area2);

portion = 1-tri_area/total_area;
% 
% id = F0(id,:);
% id = id(:);
% 
% fcolor = zeros(N,1);
% fcolor(id) = 1;
% 
% path = strcat('/',path,'/tria/',num2str(i),'/');
% if ~(exist(path))
%     mkdir(path);
% end
% 
% write_property(strcat('/',path,'/tria/',num2str(i),'/',featype,'.surface.vtk'),v, f, struct('flip',fcolor));
end