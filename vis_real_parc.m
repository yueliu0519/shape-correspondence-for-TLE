colortype = 'parc';
featype = 'xyz';
path = 'mnt/Data/liuy108/codes/tlecodes/tmi/v1';

ico =1;
slist = 7;%[2 3 4 5 6 7 8 9 11 13 14 15 18 19 20 21 22 23 24 25 26 27 30 31 32];

for i = slist
    hemi = 'rh';
    if ismember(i, [15 19 26 4 6 8])
    hemi = 'lh';
    end
    %% on particle sphere
    N = 4^ico*10+2;
    [~,f_S]=read_vtk(strcat('/',path,'/oursphere/ico',num2str(ico),'.vtk'));
    F_S=f_S+1;

    cmd1 = strcat('/',path,'/SphericalRemesh -s /',path,'/oursphere/ico',num2str(ico),'.vtk -r /',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk --bary /',path,'/vis/bary/real/',num2str(i),'baryico',num2str(ico),'.txt');
    system(cmd1);
    cmd2 = strcat('/',path,'/SphericalRemesh -s /',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.movedv.unfold.real.vtk -r /',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.sphere.vtk --bary /',path,'/vis/bary/before/',num2str(i),'baryico',num2str(ico),'.real.txt');
    system(cmd2);
    P1 = load(strcat('/',path,'/vis/bary/real/',num2str(i),'baryico',num2str(ico),'.txt'));
    P2 = load(strcat('/',path,'/vis/bary/before/',num2str(i),'baryico',num2str(ico),'.real.txt'));
    P1(:,1)=P1(:,1)+1;
    P2(:,1)=P2(:,1)+1; 
    map_1=zeros(size(P1,1),1);
    map_2=zeros(size(P2,1),1);
    nsample1 = size(P1,1);
    nsample2 = size(P2,1);

    pc_B = load(strcat('/',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.resdktparc.moved.ico.',num2str(ico),'.txt'));
%     pc_A = pc_B;

    %% for original sphere
    for j=1:nsample1
        BC_1=P1(j,2:4);
        ID_1=P1(j,1);
        % update features
        [~,M] = max(BC_1);
        map_1(j,:)=pc_B(F_S(ID_1,M),:);
    end
    [~, f_A] = read_vtk(strcat('/',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.movedv.unfold.real.vtk'));
    F_A = f_A+1;
    for j=1:nsample2
        BC_2=P2(j,2:4);
        ID_2=P2(j,1);
        % update features
        [~,M] = max(BC_2);
        map_2(j,:)=pc_B(F_A(ID_2,M),:);
    end
    %% save for quantitative test in evaluate
    fid = fopen(strcat('/mnt/Data/liuy108/codes/tlecodes/tmi/v1/tmp/ours/',featype,'/',num2str(i),'.parc.dkt.before.real.pretopost.',num2str(ico),'.txt'),'w');
    fprintf(fid,'%g\n',map_2);
    fclose(fid);
    fid = fopen(strcat('/mnt/Data/liuy108/codes/tlecodes/tmi/v1/tmp/ours/',featype,'/',num2str(i),'.parc.dkt.after.real.pretopost.',num2str(ico),'.txt'),'w');
    fprintf(fid,'%g\n',map_1);
    fclose(fid);
    % save(strcat('/',path,'/tmp/',num2str(i),'map_Bwhite.mat'),'map_1');
    % save(strcat('/',path,'/tmp/',num2str(i),'map_Awhite.mat'),'map_2');
    %% save visualization
    [v_2, f_2] = read_vtk(strcat('/',path,'/data/pre/',hemi,'/',num2str(i),'/',hemi,'.white.aff.s5.vtk'));
    F_2 = f_2 + 1;
    [v_1, f_1] = read_vtk(strcat('/',path,'/data/post/',hemi,'/',num2str(i),'/',hemi,'.white.s5.vtk'));
    F_1 = f_1 + 1;
    %% %for quanlitative test
    write_property(strcat('/',path,'/color/',featype,'/',colortype,'/original/white.after.real',num2str(i),'.ico',num2str(ico),'.pretopost.vtk'),v_1,f_1,struct('color',map_1));
    write_property(strcat('/',path,'/color/',featype,'/',colortype,'/original/white.before.real',num2str(i),'.ico',num2str(ico),'.pretopost.vtk'),v_2,f_2,struct('color',map_2));
end