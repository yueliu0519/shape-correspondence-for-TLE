function [id] = levelup(i,ico,hemi,featype,PA_new,nsample,nsample1,path,v_T,F_T,f0,F0)

if strcmp(featype, 'curvature')
    for j = 1:nsample
        mv(j,:) = PA_new(j,2:4)*v_T(F_T(PA_new(j,1),:),:);
    end
    mv = mv ./ vecnorm(mv,2,2);
    [mv,id]=unfold(mv,F0);
    write_vtk(strcat('/',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.movedv.unfold.real.vtk'),mv,f0);
    return
end
if ico <= 6

    [v1,f1]=read_vtk(strcat('/',path,'/oursphere/ico',num2str(ico+1),'.vtk'));
    F1=f1+1;

    for j = 1:nsample
        moved_v(j,:) = PA_new(j,2:4)*v_T(F_T(PA_new(j,1),:),:);
    end   

    moved_v = moved_v ./ vecnorm(moved_v,2,2);

    moved_v = unfold(moved_v,F0);
    %save unfolded sphere
    write_vtk(strcat('/',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.movedv.unfold.real.vtk'),moved_v,f0);

    neighbor = cell(length(v1), 1);

    for j = 1: size(F1, 1)
        neighbor{F1(j,1)} = [neighbor{F1(j,1)}, F1(j, 2: 3)];
        neighbor{F1(j,2)} = [neighbor{F1(j,2)}, F1(j, 1:2:3)];
        neighbor{F1(j,3)} = [neighbor{F1(j,3)}, F1(j, 1: 2)];
    end
    for j = 1: length(v1)
        neighbor{j} = unique(neighbor{j});
        neighbor{j}(neighbor{j}>nsample) = [];
    end

    % subdivision
    v = zeros(nsample1,3);
    v(1:nsample,:) = moved_v;
    for j = (nsample+1):nsample1
        v(j,:) = 0.5*(moved_v(neighbor{j}(1),:)+moved_v(neighbor{j}(2),:));
    end
    v = v ./ vecnorm(v,2,2);
    %save upleveled sphere
    write_vtk(strcat('/',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico+1),'.real.vtk'),v,f1);
    id = [];
end
if ico == 7

    for j = 1:nsample
        moved_v(j,:) = PA_new(j,2:4)*v_T(F_T(PA_new(j,1),:),:);
    end
    moved_v = moved_v ./ vecnorm(moved_v,2,2);

    [moved_v,id]=unfold(moved_v,F0);

    %save unfolded sphere
    write_vtk(strcat('/',path,'/oursphere/',featype,'/',num2str(i),'/',num2str(i),'ico',num2str(ico),'.movedv.unfold.real.vtk'),moved_v,f0);
    
    %see flipped triangles
    
end
end