%% Simulation of the dynamics of tuberculosis lesions in mice lungs
% Computational Biophysics (Enginyeria F�sica, Universitat Polit�cnica de Catalunya)
% Alexandre Justo, Jos� Javier Ruiz

function [time,z_coord,lesion03,lesion13,lesion23,lesion33,...
    vec_number_of_total_lesions,vec_number_of_lesions_plane,...
    vec_total_volume,vec_total_area,vec_mean_volume,vec_mean_area] = ...
    OrganLung(x_max,y_max,z_max,init_radius,init_age,v_mean,R1_mean,...
    R2_mean,K_mean,ini_lesions,dt,n_steps,do,init_size)


%% INITIALIZATION

% Variables and initial values

% Properties of a lesion
lesion = zeros(1,8);
% 1: x coordinate; 2: y coordinate; 3: z coordinate; 4: radius; 5: age;
% 6: growth rate; 7: max radius; 8: reproduction constant.

% Vector where number of lesions in the chosen slide (2D) is saved for
% every time step.
vec_number_of_lesions_plane = zeros(1,n_steps);

% Vector where number of total lesions 3D is saved for every time step.
vec_number_of_total_lesions = zeros(1,n_steps);

% Vector of time
time = init_age + (0:dt:((n_steps-1)*dt));

% Chose randomly a 2D slide:
z_coord = (z_max-init_size)/2+((z_max+init_size)/2-(z_max-init_size)/2)*rand(1);

% Auxiliary variable to count the number of initial lesions in the 2D slide
lesions_init_z = 0;

% Vector where the total affected volume 3D is saved for every time step
vec_total_volume = zeros(1,n_steps);

% Vector where the total affected area 2D in the chosen slide is saved for
% every time step
vec_total_area = zeros(1,n_steps);

% Vector where the mean volume 3D per lesion is saved for every time step
vec_mean_volume = zeros(1,n_steps);

% Vector where the mean area 2D per lesion in the chosen slide is saved for each value of time
vec_mean_area = zeros(1,n_steps);

% Set initial properties of lesions
for i = 1:ini_lesions                     
    lesion(i,1) = (x_max-init_size)/2+((x_max+init_size)/2-(x_max-init_size)/2)*rand(1);                % x-coordinate
    lesion(i,2) = (y_max-init_size)/2+((y_max+init_size)/2-(y_max-init_size)/2)*rand(1);                % y-coordinate
    lesion(i,3) = (z_max-init_size)/2+((z_max+init_size)/2-(z_max-init_size)/2)*rand(1);                % z-coordinate
    lesion(i,4) = init_radius*normrnd(1,0.2);	% Initial radius
    lesion(i,5) = init_age*normrnd(1,0.2);      % Initial age
    lesion(i,6) = v_mean*normrnd(1,0.2);        % Growth rate
    lesion(i,7) = R1_mean*normrnd(1,0.2);       % Maximum radius
    lesion(i,8) = K_mean*normrnd(1,0.2);        % Reproduction constant
    if abs(lesion(i,3)-z_coord) < lesion(i,4)
       lesions_init_z = lesions_init_z + 1;     % Count the initial number of lesions in the 2D slide
    end
end

% Update vectors containing the number of lesions:
vec_number_of_lesions_plane(1) = lesions_init_z;
vec_number_of_total_lesions(1) = ini_lesions;

% Output for representing
lesion03(:,:) = lesion(:,:);

%% MAIN LOOP

% Temporal evolution
for i=2:n_steps
    
number_of_lesions_plane = 0;	% Auxiliary variable to count the number of lesions in the chosen slide
less_lesions_number = 0;        % Auxiliary variable to count the number of sticking lesions that must be discounted from the total number of lesions

for j = 1:length(lesion(:,1)) 
    % Lesion age increases:
    lesion(j,5) = lesion(j,5)+dt;

    % Lesion logistic growth:
    lesion(j,4) = lesion(j,4)*(1+dt*lesion(j,6)*(1-(lesion(j,4)/lesion(j,7))^2));
end

index_vec = []; % Auxiliary vector to delete a lesion which has been merged with another one

% Cases of merging and sticking, all lesion sizes and locations must be
% compared. With this double loop any pair of lesions is repeated and any
% lesion is compared with itself.
for m = 1:length(lesion(:,1))
    for n = (m+1):length(lesion(:,1))
        % Comparing lesion m with lesion n:
        % Compute the distance between their centers = dmn:
        xmn = lesion(m,1)-lesion(n,1);
        ymn = lesion(m,2)-lesion(n,2);
        zmn = lesion(m,3)-lesion(n,3);
        dmn = sqrt(xmn^2+ymn^2+zmn^2);
        
        % Their respective radii are:
        rm = lesion(m,4); rn = lesion(n,4);
        
        % If dmn < rm+rn there appears merging or sticking:
        
        % rm~=0 && rn~=0 is imposed to not consider lesions which have been
        % merged and that must be deleted later. If they were deleted during
        % this loop, then there would be a problem while calling a
        % component of the "lesion" matrix, since its size would have changed
        if dmn < (rm+rn) && rm ~= 0 && rn ~= 0
            
            % If there is merging, update the properties of lesion m
            % and delete the properties of lesion n
            if dmn < max(rm,rn)
                
                % The location of the merged lesion is located at the
                % mass center of both lesions:
                lesion(m,1) = 1/(rm^3+rn^3)*(rm^3*lesion(m,1)+rn^3*lesion(n,1));
                lesion(m,2) = 1/(rm^3+rn^3)*(rm^3*lesion(m,2)+rn^3*lesion(n,2));
                lesion(m,3) = 1/(rm^3+rn^3)*(rm^3*lesion(m,3)+rn^3*lesion(n,3));
                
                % Set its properties
                lesion(m,4) = (rm^3+rn^3)^(1/3);            % Radius
                lesion(m,5) = max(lesion(m,5),lesion(n,5));	% The age is the maximum age of the two merging lesions
                lesion(m,6) = v_mean*normrnd(1,0.2);        % Growth velocity
                lesion(m,7) = R2_mean*normrnd(1,0.2);       % The maximum radius can be larger than the one of a lesion which has appear without merging
                lesion(m,8) = K_mean*normrnd(1,0.2);        % Reproduction constant
                
                % Set the properties of lesion n to 0, then, in the
                % loop it will not be considered. This lesion will be
                % deleted later
                lesion(n,:) = zeros(1,8);
                
                % Save the lesion that must be deleted later
                index_vec = [index_vec n];
                
            % If there is no merging, but sticking:
            else
                % Both lesions evolve individually, but one lesion must
                % be discounted from the total number:
                less_lesions_number = less_lesions_number+1;
                
                % Check if there is sticking in the chosen slide:
                if abs(lesion(n,3)-z_coord) < lesion(n,4) && abs(lesion(m,3)-z_coord) < lesion(m,4)
                    % Then, one lesion must be discounted from the number
                    % of lesions present in the slide:
                    number_of_lesions_plane = number_of_lesions_plane-1;
                end
            end
        end
    end
end

% Remove lesions which have been merged with another lesion
lesion(index_vec,:) = []; % Remove rows which are 0

% Total number of lesions (and sublesions):
len = length(lesion(:,1));

% Depending on the age of each lesion, daughter lesions can appear
% In this loop is also computed the number of lesions (and sublesions)
% present in the chosen slide
for j = 1:len
    
    % If the age of the lesion is between 14 and 28 days, then there may
    % appear daughter lesions:
    if lesion(j,5) > 14 && lesion(j,5) < 28
        
        % Number of daughter lesions produced:
        Ndaughter = fix(lesion(j,8)*lesion(j,4)*(2-lesion(j,5)/14));
        
        % Set the initial properties of each daughter lesion
        if Ndaughter ~= 0
            % Add Ndaughter lesions at the end of the matrix:
            final=length(lesion(:,1));
            for w = 1:Ndaughter
                le=final+w;
                % Set the properties of the daughter lesion
                lesion(le,1) = lesion(j,1)+do*gamrnd(10,0.1);	% x-coordinate
                lesion(le,2) = lesion(j,2)+do*gamrnd(10,0.1);	% y-coordinate
                lesion(le,3) = lesion(j,3)+do*gamrnd(10,0.1);	% z-coordinate
                lesion(le,4) = init_radius*normrnd(1,0.2);      % Initial radius
                lesion(le,5) = init_age*normrnd(1,0.2);         % Initial age
                lesion(le,6) = v_mean*normrnd(1,0.2);           % Growth velocity
                lesion(le,7) = R1_mean*normrnd(1,0.2);          % Maximum radius
                lesion(le,8) = K_mean*normrnd(1,0.2);           % Reproduction constant
                
                % Check if the new daughter lesion is present in the chosen
                % slide
                if  abs(lesion(le,3)-z_coord) < lesion(le,4)
                    % Then, add one lesion present in the chosen 2D slide
                    number_of_lesions_plane = number_of_lesions_plane+1;
                end
            end     
        end  
    end
    
    % Compute the total number of lesions and sublesions present in the 2D
    % slide:
    if  abs(lesion(j,3)-z_coord)<lesion(j,4)
        number_of_lesions_plane = number_of_lesions_plane+1;
    end
end
   
% Update the number of lesions present in the chosen 2D slide
vec_number_of_lesions_plane(i) = number_of_lesions_plane;

% Update the total number of lesions, which is the total number of lesions
% and sublesions minus the number of lesions that have been merged with
% another one (actually, the number of lesions that are merged are two, but
% only one lesion must be discounted, as considered before)
vec_number_of_total_lesions(i) = length(lesion(:,1))-less_lesions_number;

vec_volume_lesion = zeros(1,length(lesion(:,1)));
% Auxiliary matrix to compute the total affected area
U = [];

for j = 1:length(lesion(:,1))
    r = lesion(j,4);
    vec_volume_lesion(j) = 4/3*pi*r^3;

    if abs(lesion(j,3)-z_coord) < lesion(j,4)
        R = sqrt(r^2-(z_coord-lesion(j,3))^2);
        U = [U;[lesion(j,1),lesion(j,2),R]];
    end
end

if ~isempty(U)
    % Compute the total affected area in the chosen slide by means of
    % Monte Carlo method
    Total_affected_area = MC(U,10^8);
else
    Total_affected_area = 0;
end

vec_total_volume(i) = sum(vec_volume_lesion);
vec_total_area(i) = Total_affected_area;

vec_mean_volume(i) = mean(vec_volume_lesion);
vec_mean_area(i) = Total_affected_area/length(lesion(:,1));

% Outputs for representing
if i == round(1*n_steps/3)
    lesion13(:,:) = lesion(:,:);
elseif i == round(2*n_steps/3)
    lesion23(:,:) = lesion(:,:);
elseif i == round(3*n_steps/3)
    lesion33(:,:) = lesion(:,:);
end
i
end
end