path = 'mnt/Data/liuy108/codes/tlecodes/tmi/v1';
featype = 'curvature';
colortype = 'parc';
ico = 7;
d = 25;
total_mean_freesurfer = zeros(d,32);
total_mean_hsd = zeros(d,32);
total_mean_icp = zeros(d,32);
total_mean_ours = zeros(d,32);
total_mean_cpd = zeros(d,32);
total_hd_cpd = zeros(d,32);
total_hd_freesurfer = zeros(d,32);
total_hd_hsd = zeros(d,32);
total_hd_icp = zeros(d,32);
total_hd_ours = zeros(d,32);
for sub = [2 3 4 5 6 7 8 9 11 13 14 18 19 20 21 22 23 24 25 26 27 30 31]
    hemi = 'rh';
    if ismember(sub, [15 19 26 4 6 8])
    hemi = 'lh';
    end
    [v_b,f_b] = read_vtk(strcat('/',path,'/color/',featype,'/',colortype,'/original/white.before.sim',num2str(sub),'.ico',num2str(ico),'.pretopost.vtk'));
    F_b = f_b+1;
    [v_a,f_a] = read_vtk(strcat('/',path,'/color/',featype,'/',colortype,'/original/white.after.sim',num2str(sub),'.ico',num2str(ico),'.pretopost.vtk'));
    F_a = f_a+1;
    [v_a1,f_a1] = read_vtk(strcat('/',path,'/baseline/resultfs/',hemi,'.',num2str(sub),'.after.white.fs.sim.pretopost.new.vtk'));
    F_a1 = f_a1+1;
    [v_a2,f_a2] = read_vtk(strcat('/',path,'/baseline/resulthsd/',hemi,'.',num2str(sub),'.after.white.hsd.sim.pretopost.new.vtk'));
    F_a2 = f_a2+1;
    [v_a3,f_a3] = read_vtk(strcat('/',path,'/baseline/resultcpd/',num2str(sub),'after.white.cpd.parc.sim.posttopre.vtk'));
    F_a3 = f_a3+1;
    [v_a4,f_a4] = read_vtk(strcat('/',path,'/baseline/resulticp/',num2str(sub),'after.white.icp.parc.sim.posttopre.vtk'));
    F_a4 = f_a4+1;

    map_A = load(strcat('/',path,'/tmp/ours/',featype,'/',num2str(sub),'.parc.dkt.after.sim.pretopost.txt'));
    map_B = load(strcat('/',path,'/tmp/ours/',featype,'/',num2str(sub),'.parc.dkt.before.sim.pretopost.txt'));
    map_fs = load(strcat('/',path,'/tmp/fs/',num2str(sub),'.sim.post.dktparc.txt'));
    map_hsd = load(strcat('/',path,'/tmp/hsd/',num2str(sub),'.sim.post.dktparc.txt'));
    map_icp = load(strcat('/',path,'/tmp/icp/',num2str(sub),'.parc.sim.posttopre.txt'));
    map_cpd = load(strcat('/',path,'/tmp/cpd/',num2str(sub),'.parc.sim.posttopre.txt'));
    %before
    nn_b = cell(length(v_b), 1);
    for j = 1: size(F_b, 1)
        nn_b{F_b(j,1)} = [nn_b{F_b(j,1)}, F_b(j, 2: 3)];
        nn_b{F_b(j,2)} = [nn_b{F_b(j,2)}, F_b(j, 1:2:3)];
        nn_b{F_b(j,3)} = [nn_b{F_b(j,3)}, F_b(j, 1: 2)];
    end
    for j = 1: length(v_b)
        nn_b{j} = unique(nn_b{j});
    end
    %ours
    nn_a = cell(length(v_a), 1);
    for j = 1: size(F_a, 1)
        nn_a{F_a(j,1)} = [nn_a{F_a(j,1)}, F_a(j, 2: 3)];
        nn_a{F_a(j,2)} = [nn_a{F_a(j,2)}, F_a(j, 1:2:3)];
        nn_a{F_a(j,3)} = [nn_a{F_a(j,3)}, F_a(j, 1: 2)];
    end
    for j = 1: length(v_a)
        nn_a{j} = unique(nn_a{j});
    end
    %fs
    nn_a1 = cell(length(v_a1), 1);
    for j = 1: size(F_a1, 1)
        nn_a1{F_a1(j,1)} = [nn_a1{F_a1(j,1)}, F_a1(j, 2: 3)];
        nn_a1{F_a1(j,2)} = [nn_a1{F_a1(j,2)}, F_a1(j, 1:2:3)];
        nn_a1{F_a1(j,3)} = [nn_a1{F_a1(j,3)}, F_a1(j, 1: 2)];
    end
    for j = 1: length(v_a1)
        nn_a1{j} = unique(nn_a1{j});
    end
    %hsd
    nn_a2 = cell(length(v_a2), 1);
    for j = 1: size(F_a2, 1)
        nn_a2{F_a2(j,1)} = [nn_a2{F_a2(j,1)}, F_a2(j, 2: 3)];
        nn_a2{F_a2(j,2)} = [nn_a2{F_a2(j,2)}, F_a2(j, 1:2:3)];
        nn_a2{F_a2(j,3)} = [nn_a2{F_a2(j,3)}, F_a2(j, 1: 2)];
    end
    for j = 1: length(v_a2)
        nn_a2{j} = unique(nn_a2{j});
    end
    %cpd
    nn_a3 = cell(length(v_a3), 1);
    for j = 1: size(F_a3, 1)
        nn_a3{F_a3(j,1)} = [nn_a3{F_a3(j,1)}, F_a3(j, 2: 3)];
        nn_a3{F_a3(j,2)} = [nn_a3{F_a3(j,2)}, F_a3(j, 1:2:3)];
        nn_a3{F_a3(j,3)} = [nn_a3{F_a3(j,3)}, F_a3(j, 1: 2)];
    end
    for j = 1: length(v_a3)
        nn_a3{j} = unique(nn_a3{j});
    end

    %icp
    nn_a4 = cell(length(v_a4), 1);
    for j = 1: size(F_a4, 1)
        nn_a4{F_a4(j,1)} = [nn_a4{F_a4(j,1)}, F_a4(j, 2: 3)];
        nn_a4{F_a4(j,2)} = [nn_a4{F_a4(j,2)}, F_a4(j, 1:2:3)];
        nn_a4{F_a4(j,3)} = [nn_a4{F_a4(j,3)}, F_a4(j, 1: 2)];
    end
    for j = 1: length(v_a4)
        nn_a4{j} = unique(nn_a4{j});
    end
    %label ID
    
    parc_rh = [3 4 6 9 11 12 13 14 15 18 19 20 21 22 23 24 25 26 27 28 29 30 32 33 36];
    q = 1;
    for p = parc_rh
        % prac 1
        ind1 = find(map_B == p);
        ind2 = find(map_A == p);
        ind3 = find(map_fs == p);
        ind4 = find(map_hsd == p);
        ind5 = find(map_cpd == p);
        ind6 = find(map_icp == p);
        %before
        k = 1; 
        bdry_b = [];
        for i = ind1'
            if numel(unique(map_B(nn_b{i})))>1
                bdry_b(k,:) = v_b(i,:);
                k = k + 1;
            end
        end
        %ours
        k = 1;  
        bdry_a = [];
        for i = ind2'
            if numel(unique(map_A(nn_a{i})))>1
                bdry_a(k,:) = v_a(i,:);
                k = k + 1;
            end
        end
        %freesurfer
        k = 1;  
        bdry_a1 = [];
        for i = ind3'
            if numel(unique(map_fs(nn_a1{i})))>1
                bdry_a1(k,:) = v_a1(i,:);
                k = k + 1;
            end
        end
        %sphere
        k = 1;  
        bdry_a2 = [];
        for i = ind4'
            if numel(unique(map_hsd(nn_a2{i})))>1
                bdry_a2(k,:) = v_a2(i,:);
                k = k + 1;
            end
        end
        %cpd
        k = 1;  
        bdry_a3 = [];
        for i = ind5'
            if numel(unique(map_cpd(nn_a3{i})))>1
                bdry_a3(k,:) = v_a3(i,:);
                k = k + 1;
            end
        end
        %icp
        k = 1;  
        bdry_a4 = [];
        for i = ind6'
            if numel(unique(map_icp(nn_a4{i})))>1
                bdry_a4(k,:) = v_a4(i,:);
                k = k + 1;
            end
        end
        %ours number
        [~,dist_ab]=dsearchn(bdry_a,bdry_b);
        [~,dist_ba]=dsearchn(bdry_b,bdry_a);
        mean_ours(q,1) = (mean(dist_ab)+mean(dist_ba))/2;
        hd_ours(q,1) = max(max(dist_ab),max(dist_ba));
        %fs
        [~,dist_ab1]=dsearchn(bdry_a1,bdry_b);
        [~,dist_ba1]=dsearchn(bdry_b,bdry_a1);
        mean_freesurfer(q,1) = (mean(dist_ab1)+mean(dist_ba1))/2;
        hd_freesurfer(q,1) = max(max(dist_ab1),max(dist_ba1));
        %hsd
        [~,dist_ab2]=dsearchn(bdry_a2,bdry_b);
        [~,dist_ba2]=dsearchn(bdry_b,bdry_a2);
        mean_hsd(q,1) = (mean(dist_ab2)+mean(dist_ba2))/2;
        hd_hsd(q,1) = max(max(dist_ab2),max(dist_ba2));
        %cpd
        [~,dist_ab3]=dsearchn(bdry_a3,bdry_b);
        [~,dist_ba3]=dsearchn(bdry_b,bdry_a3);
        mean_cpd(q,1) = (mean(dist_ab3)+mean(dist_ba3))/2;
        hd_cpd(q,1) = max(max(dist_ab3),max(dist_ba3));
        %icp
        [~,dist_ab4]=dsearchn(bdry_a4,bdry_b);
        [~,dist_ba4]=dsearchn(bdry_b,bdry_a4);
        mean_icp(q,1) = (mean(dist_ab4)+mean(dist_ba4))/2;
        hd_icp(q,1) = max(max(dist_ab4),max(dist_ba4));

        q = q + 1;
    end
    total_mean_freesurfer(:,sub) = mean_freesurfer;
    total_mean_hsd(:,sub) = mean_hsd;
    total_mean_icp(:,sub) = mean_icp;
    total_mean_ours(:,sub) = mean_ours;
    total_mean_cpd(:,sub) = mean_cpd;
    total_hd_freesurfer(:,sub) = hd_freesurfer;
    total_hd_hsd(:,sub) = hd_hsd;
    total_hd_icp(:,sub) = hd_icp;
    total_hd_ours(:,sub) = hd_ours;
    total_hd_cpd(:,sub) = hd_cpd;
end

%% plot
% load('/mnt/Data/liuy108/codes/tlecodes/demo2/hd_ours.mat');
% load('/mnt/Data/liuy108/codes/tlecodes/demo2/hd_freesurfer.mat');
% load('/mnt/Data/liuy108/codes/tlecodes/demo2/hd_sphere.mat');
% for i = 1:length(parc_rh)
%     hd(i,:) = [hd_ours(i) hd_freesurfer(i) hd_sphere(i)];
% end 
% figure;
% bar(parc_rh, hd);