% define variables

%% optimization
for i = [3]  %ID of subject
    i
    hemi = 'rh';
    if ismember(i, [15 19 26 4 6 8])
    hemi = 'lh';
    end
    path = 'mnt/Data/liuy108/codes/tlecodes/tmi/v1';
    diary (strcat('/',path,'/result/',num2str(i),'log.txt'));

    featype = 'xyz';
    [v_A,f_A]=read_vtk(strcat('/',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.white.aff.vtk'));
    F_A=f_A+1;
    F_A=double(F_A);
    TR_A=triangulation(F_A,v_A);
    %compute area of tria for GD
    a=sqrt(sum((v_A(F_A(:,1),:)-v_A(F_A(:,2),:)).^2,2));
    b=sqrt(sum((v_A(F_A(:,2),:)-v_A(F_A(:,3),:)).^2,2));
    c=sqrt(sum((v_A(F_A(:,1),:)-v_A(F_A(:,3),:)).^2,2));
    s=0.5*(a+b+c);
    area=sqrt(s.*(s-a).*(s-b).*(s-c));
    [v_T,f_T]=read_vtk(strcat('/',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk'));
    F_T=f_T+1;
    for ico = 0:7 %level of icosahedron
        ico
        nsample=4^ico*10+2; %number of particles
        nsample1 = 4^(ico+1)*10+2;
        [~,f0]=read_vtk(strcat('/',path,'/oursphere/ico',num2str(ico),'.vtk'));
        F0=f0+1;
        if ico ==0
            cmd = strcat('/',path,'/SphericalRemesh -s /',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk -r /',path,'/oursphere/ico',num2str(ico),'.vtk --bary /',path,'/bary/before/',num2str(i),'/',num2str(i),'baryico',num2str(ico),'real.txt');
        else
            cmd = strcat('/',path,'/SphericalRemesh -s /',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk -r /',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.real.vtk --bary /',path,'/bary/before/',num2str(i),'/',num2str(i),'baryico',num2str(ico),'real.txt');
        end

        system(cmd);
     
        system(strcat('/',path,'/SphericalRemesh -i /',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.white.vtk -o /',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.white.res.vtk -s /',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk -r /',path,'/oursphere/ico',num2str(ico),'.vtk'));
        map_B = read_vtk(strcat('/',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.white.res.vtk'));

        [~,d] = dsearchn(v_A,map_B);
        R=1./(1+exp(-(7-d)));
    % R=ones(size(PA,1),1);

%% gradient descend
        niter=2;
        eta = 1;
        alpha=0.5;
        [PA_new] = GD(niter,eta,alpha,nsample,i,ico,hemi,featype,area,map_B,R,path,TR_A,v_A,F_A);
        save(strcat('/',path,'/bary/PA/realPA_new',num2str(i),'.mat'),'PA_new');%save PA_new of xyz for curvature
    % plot(1:k-1,energy,'-*');
%% levelup for xyz
        [id] = levelup(i,ico,hemi,featype,PA_new,nsample,nsample1,path,v_T,F_T,f0,F0);  %id is used in fliptria

%% compute and see flipped triangles of xyz
        N = 4^ico*10+2;
        [proportion_xyz] = fliptria(hemi,i,ico,featype,F0,id,N,path);
        fprintf('The proportion of flipped triangle area with %s is %.4f\n', featype,proportion_xyz);
    end    
    
%% add mean curvature
    featype = 'curvature'
    ico = 7;
    nsample = 4^ico*10+2;
    nsample1 = 4^(ico+1)*10+2;
    system(strcat('/',path,'/SphericalRemesh -i /',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.white.vtk -o /',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.white.res.vtk -s /',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk -r /',path,'/oursphere/ico',num2str(ico),'.vtk'));
    map_B = read_vtk(strcat('/',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.white.res.vtk'));

    [~,d] = dsearchn(v_A,map_B);
    R=1./(1+exp(-(7-d)));
    load(strcat('/',path,'/bary/PA/realPA_new',num2str(i),'.mat'));
    for j = 1:nsample
        moved_v(j,:) = PA_new(j,2:4)*v_A(F_A(PA_new(j,1),:),:);
    end
    X = map_B(R>0.9,:);
    Y = moved_v(R>0.9,:);
    [~,~,tran]=procrustes(X, Y, 'reflection',false);
    v_a = tran.b*v_A*tran.T + tran.c(1,:);

    [PA_newH] = meancurvature(i,ico,hemi,v_a,area,path,TR_A,F_A);
    save(strcat('/',path,'/bary/PA/realPA_newH',num2str(i),'.mat'),'PA_newH');%save PA_new of xyz for curvature
% levelup for mean curvature
    [id] = levelup(i, ico, hemi,featype,PA_newH,nsample,nsample1,path,v_T,F_T,f0,F0); 
%% compute and see flipped triangles of H
    [proportion_H] = fliptria(hemi,i,ico,featype,F0,id,N,path);
    fprintf('The proportion of flipped triangle area with %s is %.4f\n', featype,proportion_H);
%% visualization of our method
%     colortype = 'parc';
%     featype = 'xyz';
%     pc_B = [];
%     pc_B = vis(colortype,featype,i,hemi,ico,path);% pc_B is compute in vis for 'checkerboard' and load for 'parc'
%     save(strcat('/',path,'/color/pcB/',colortype,'pc_Bxyz',num2str(i),'.mat'),'pc_B');%save for baseline checkerboard
%     featype = 'curvature';
%     pc_B = [];
%     pc_B = vis(colortype,featype,i,hemi,ico,path);% pc_B is compute in vis for 'checkerboard' and load for 'parc'
%     save(strcat('/',path,'/color/pcB/',colortype,'pc_Bcurvature',num2str(i),'.mat'),'pc_B');%save for baseline checkerboard
%     clc;clear;
end
