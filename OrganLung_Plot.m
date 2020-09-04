%% Simulation of the dynamics of tuberculosis lesions in mice lungs
% Computational Biophysics (Enginyeria Física, Universitat Politècnica de Catalunya)
% Alexandre Justo, José Javier Ruiz

%% 3D and 2D representations
x_max = 15;             % 1. Number of spatial cells (x)   1 mm/cell
y_max = 15;             % 2. Number of spatial cells (y)   1 mm/cell
z_max = 15;             % 3. Number of spatial cells (z)   1 mm/cell
init_radius = 0.15/2;	% 4. Mean initial radius of initial lesions [mm]
init_age = 14;          % 5. Mean initial age of initial lesions [days]
v_mean = 0.15;          % 6. Mean value for growth rate [day]^(-1)
R1_mean = 1;            % 7. a) Mean maximum radius before merging [mm]
R2_mean = 10;           % 7. b) Mean maximum radius after merging [mm]
K_mean = 5;             % 8. Reproduction constant [mm]^(-1)
ini_lesions = 50;       % Number of initial lesions
dt = 0.5;               % Time step of half a day
n_steps = 23;           % Number of time steps, final integration time of (init_age + (n_steps-1)*dt) days
do = 1;                 % Constant for determining the distance at which a daughter lesion appears [mm]
init_size = 15;         % Side of the centered cube where all the initial lesions are located

[time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size);

[x,y,z] = sphere;
for j = 1:length(lesion03(:,1))
    r = lesion03(j,4);
    figure(1)
    subplot(2,2,1)
    surf(r*x+lesion03(j,1),r*y+lesion03(j,2),r*z+lesion03(j,3))
    title(['Lesions 3D, day ' num2str(time(round(1+0*n_steps/3)))])
    xlim ([-1 x_max+1]); ylim ([-1 y_max+1]); zlim ([-1 z_max+1]);
    xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]')
    hold on
    figure(2)
    subplot(2,2,1)
    if abs(lesion03(j,3)-z_coord) < lesion03(j,4)
        R = sqrt(r^2-(z_coord-lesion03(j,3))^2);
        pos = [(lesion03(j,1)-R) (lesion03(j,2)-R) 2*R 2*R]; 
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0 .5 .5])
        hold on
    end
    title(['Lesions 2D, slide z=' num2str(z_coord) ', day ' num2str(time(round(1+0*n_steps/3)))]);
    xlabel('x [mm]'); ylabel('y [mm]')
    xlim ([-1 x_max+1]); ylim ([-1 y_max+1]); zlim ([-1 z_max+1]);
end

for j = 1:length(lesion13(:,1))
    r = lesion13(j,4);
    figure(1)
    subplot(2,2,2)
    surf(r*x+lesion13(j,1),r*y+lesion13(j,2),r*z+lesion13(j,3))
    title(['Lesions 3D, day ' num2str(time(round(1*n_steps/3)))])
    xlim ([-1 x_max+1]); ylim ([-1 y_max+1]); zlim ([-1 z_max+1]);
    xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]')
    hold on
    figure(2)
    subplot(2,2,2)
    if abs(lesion13(j,3)-z_coord) < lesion13(j,4)
        R = sqrt(r^2-(z_coord-lesion13(j,3))^2);
        pos = [(lesion13(j,1)-R) (lesion13(j,2)-R) 2*R 2*R]; 
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0 .5 .5])
        hold on
    end
    title(['Lesions 2D, slide z=' num2str(z_coord) ', day ' num2str(time(round(1*n_steps/3)))]);
    xlabel('x [mm]'); ylabel('y [mm]')
    xlim ([-1 x_max+1]); ylim ([-1 y_max+1]); zlim ([-1 z_max+1]);
end

for j = 1:length(lesion23(:,1))
    r = lesion23(j,4);
    figure(1)
    subplot(2,2,3)
    surf(r*x+lesion23(j,1),r*y+lesion23(j,2),r*z+lesion23(j,3))
    title(['Lesions 3D, day ' num2str(time(round(2*n_steps/3)))])
    xlim ([-1 x_max+1]); ylim ([-1 y_max+1]); zlim ([-1 z_max+1]);
    xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]')
    hold on
    figure(2)
    subplot(2,2,3)
    if abs(lesion23(j,3)-z_coord) < lesion23(j,4)
        R = sqrt(r^2-(z_coord-lesion23(j,3))^2);
        pos = [(lesion23(j,1)-R) (lesion23(j,2)-R) 2*R 2*R]; 
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0 .5 .5])
        hold on
    end
    title(['Lesions 2D, slide z=' num2str(z_coord) ', day ' num2str(time(round(2*n_steps/3)))]);
    xlabel('x [mm]'); ylabel('y [mm]')
    xlim ([-1 x_max+1]); ylim ([-1 y_max+1]); zlim ([-1 z_max+1]);
end

for j = 1:length(lesion33(:,1))
    r = lesion33(j,4);
    figure(1)
    subplot(2,2,4)
    surf(r*x+lesion33(j,1),r*y+lesion33(j,2),r*z+lesion33(j,3))
    title(['Lesions 3D, day ' num2str(time(round(3*n_steps/3)))])
    xlim ([-1 x_max+1]); ylim ([-1 y_max+1]); zlim ([-1 z_max+1]);
    xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]')
    hold on
    figure(2)
    subplot(2,2,4)
    if abs(lesion33(j,3)-z_coord) < lesion33(j,4)
        R = sqrt(r^2-(z_coord-lesion33(j,3))^2);
        pos = [(lesion33(j,1)-R) (lesion33(j,2)-R) 2*R 2*R]; 
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0 .5 .5])
        hold on
    end
    title(['Lesions 2D, slide z=' num2str(z_coord) ', day ' num2str(time(round(3*n_steps/3)))]);
    xlabel('x [mm]'); ylabel('y [mm]')
    xlim ([-1 x_max+1]); ylim ([-1 y_max+1]); zlim ([-1 z_max+1]);
end

figure(3);
subplot(3,2,1);
plot(time,vec_number_of_total_lesions);
hold on
subplot(3,2,2);
plot(time,vec_number_of_lesions_plane);
hold on
subplot(3,2,3);
plot(time,vec_total_volume);
hold on
subplot(3,2,4);
plot(time,vec_total_area);
hold on
subplot(3,2,5);
plot(time,vec_mean_volume);
hold on
subplot(3,2,6);
plot(time,vec_mean_area);
hold on

figure(4);
subplot(3,2,1);
plot(time,vec_number_of_total_lesions);
hold on
subplot(3,2,2);
plot(time,vec_number_of_lesions_plane);
hold on
subplot(3,2,3);
plot(time,vec_total_volume);
hold on
subplot(3,2,4);
plot(time,vec_total_area);
hold on
subplot(3,2,5);
plot(time,vec_mean_volume);
hold on
subplot(3,2,6);
plot(time,vec_mean_area);
hold on

figure(5);
subplot(3,2,1);
plot(time,vec_number_of_lesions_plane,'ok');
hold on
subplot(3,2,2);
plot(time,vec_mean_area,'ok');
hold on
subplot(3,2,3);
plot(time,vec_number_of_lesions_plane,'ok');
hold on
subplot(3,2,4);
plot(time,vec_mean_area,'ok');
hold on
subplot(3,2,5);
plot(time,vec_number_of_lesions_plane,'ok');
hold on
subplot(3,2,6);
plot(time,vec_mean_area,'ok');
hold on

%% See how the number of initial lesions affects
ini_lesions = 100;

[time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size);

figure(3);
subplot(3,2,1);
plot(time,vec_number_of_total_lesions);
hold on
subplot(3,2,2);
plot(time,vec_number_of_lesions_plane);
hold on
subplot(3,2,3);
plot(time,vec_total_volume);
hold on
subplot(3,2,4);
plot(time,vec_total_area);
hold on
subplot(3,2,5);
plot(time,vec_mean_volume);
hold on
subplot(3,2,6);
plot(time,vec_mean_area);
hold on

ini_lesions = 150;

[time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size);

figure(3);
subplot(3,2,1);
plot(time,vec_number_of_total_lesions);
hold off
legend('50 initial lesions','100 initial lesions','150 initial lesions','Location','northwest')
xlabel('Time [days]'); title('Number of total lesions'); xlim([time(1)+dt time(end)])
subplot(3,2,2);
plot(time,vec_number_of_lesions_plane);
hold off
legend('50 initial lesions','100 initial lesions','150 initial lesions','Location','northwest')
xlabel('Time [days]'); title(['Number of lesions in slide z=' num2str(z_coord)]); xlim([time(1)+dt time(end)])
subplot(3,2,3);
plot(time,vec_total_volume);
hold off
legend('50 initial lesions','100 initial lesions','150 initial lesions','Location','northwest')
xlabel('Time [days]'); title('Total affected volume [mm^3]'); xlim([time(1)+dt time(end)])
subplot(3,2,4);
plot(time,vec_total_area);
hold off
legend('50 initial lesions','100 initial lesions','150 initial lesions','Location','northwest')
xlabel('Time [days]'); title('Total affected area [mm^2]'); xlim([time(1)+dt time(end)])
subplot(3,2,5);
plot(time,vec_mean_volume);
hold off
legend('50 initial lesions','100 initial lesions','150 initial lesions','Location','northwest')
xlabel('Time [days]'); title('Mean volume per lesion [mm^3]'); xlim([time(1)+dt time(end)])
subplot(3,2,6);
plot(time,vec_mean_area);
hold off
legend('50 initial lesions','100 initial lesions','150 initial lesions','Location','northwest')
xlabel('Time [days]'); title('Mean area per lesion [mm^2]'); xlim([time(1)+dt time(end)])

ini_lesions = 50;

%% See how the initial mean distance affects

init_size = 10;

[time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size);

figure(4);
subplot(3,2,1);
plot(time,vec_number_of_total_lesions);
hold on
subplot(3,2,2);
plot(time,vec_number_of_lesions_plane);
hold on
subplot(3,2,3);
plot(time,vec_total_volume);
hold on
subplot(3,2,4);
plot(time,vec_total_area);
hold on
subplot(3,2,5);
plot(time,vec_mean_volume);
hold on
subplot(3,2,6);
plot(time,vec_mean_area);
hold on

init_size = 6;

[time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size);

figure(4);
subplot(3,2,1);
plot(time,vec_number_of_total_lesions);
hold off
legend('15x15x15','10x10x10','6x6x6','Location','northwest')
xlabel('Time [days]'); title('Number of total lesions'); xlim([time(1)+dt time(end)])
subplot(3,2,2);
plot(time,vec_number_of_lesions_plane);
hold off
legend('15x15x15','10x10x10','6x6x6','Location','northwest')
xlabel('Time [days]'); title(['Number of lesions in slide z=' num2str(z_coord)]); xlim([time(1)+dt time(end)])
subplot(3,2,3);
plot(time,vec_total_volume);
hold off
legend('15x15x15','10x10x10','6x6x6','Location','northwest')
xlabel('Time [days]'); title('Total affected volume [mm^3]'); xlim([time(1)+dt time(end)])
subplot(3,2,4);
plot(time,vec_total_area);
hold off
legend('15x15x15','10x10x10','6x6x6','Location','northwest')
xlabel('Time [days]'); title('Total affected area [mm^2]'); xlim([time(1)+dt time(end)])
subplot(3,2,5);
plot(time,vec_mean_volume);
hold off
legend('15x15x15','10x10x10','6x6x6','Location','northwest')
xlabel('Time [days]'); title('Mean volume per lesion [mm^3]'); xlim([time(1)+dt time(end)])
subplot(3,2,6);
plot(time,vec_mean_area);
hold off
legend('15x15x15','10x10x10','6x6x6','Location','northwest')
xlabel('Time [days]'); title('Mean area per lesion [mm^2]'); xlim([time(1)+dt time(end)])

init_size = 15;

%% Test FIRST hypothesis (slower growth rate)
v_mean = 0.1*v_mean;

[time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size);

figure(5);
subplot(3,2,1);
plot(time,vec_number_of_lesions_plane,'ob');
hold off
legend('Original model','Small inflammatory response','Location','northwest')
xlabel('Time [days]'); title(['Number of lesions in slide z=' num2str(z_coord)]); xlim([time(1)+dt time(end)])
subplot(3,2,2);
plot(time,vec_mean_area,'ob');
hold off
legend('Original model','Small inflammatory response','Location','northwest')
xlabel('Time [days]'); title('Mean area per lesion [mm^2]'); xlim([time(1)+dt time(end)])

v_mean = 10*v_mean;

%% Test SECOND hypothesis (endogeneous reinfection)

[time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung2(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size);

figure(5);
subplot(3,2,3);
plot(time,vec_number_of_lesions_plane,'ob');
hold off
legend('Original model','No endogeneous reinfection','Location','northwest')
xlabel('Time [days]'); title(['Number of lesions in slide z=' num2str(z_coord)]); xlim([time(1)+dt time(end)])
subplot(3,2,4);
plot(time,vec_mean_area,'ob');
hold off
legend('Original model','No endogeneous reinfection','Location','northwest')
xlabel('Time [days]'); title('Mean area per lesion [mm^2]'); xlim([time(1)+dt time(end)])

%% Test THIRD hypothesis (coalescence)

[time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung3(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size);

figure(5);
subplot(3,2,5);
plot(time,vec_number_of_lesions_plane,'ob');
hold off
legend('Original model','No coalescence','Location','northwest')
xlabel('Time [days]'); title(['Number of lesions in slide z=' num2str(z_coord)]); xlim([time(1)+dt time(end)])
subplot(3,2,6);
plot(time,vec_mean_area,'ob');
hold off
legend('Original model','No coalescence','Location','northwest')
xlabel('Time [days]'); title('Mean area per lesion [mm^2]'); xlim([time(1)+dt time(end)])