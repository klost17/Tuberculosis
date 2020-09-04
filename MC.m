function area=MC(U,N)
% U= matrix with x and y positions, and radius in the chosen slide, each
% row corresponds to one lesion
% N= number of random points released

num=length(U(:,1)); %number of lesions
favourable_cases=0; % number of points that fall into an affected area

% for each point released, check if it has fallen into an affected area
for i=1:N
    % coordinates of the random point released
    xrand=15*rand;
    yrand=15*rand;
    
    % check if the released point falls into an affected area.
    % loop for all the lesions
    j=1;
    while j<=num 
        x=U(j,1); % x position of the lesion's center
        y=U(j,2); % y position of the lesion's center
        d=sqrt((x-xrand)^2+(y-yrand)^2); % distance between the random point and the center of the lesion
        
        % if the distance between the random point and the center is lesser
        % than the radius of the lesion, then it falls into an affected
        % area:
        if d<U(j,3) 
            favourable_cases=favourable_cases+1;
            j=num+1; % break the loop, since we know that the point has fallen into an affected region
            
        % otherwise, consider another lesion
        else
            j=j+1;
        end
    end
end

% total affected area:
area=favourable_cases/N*15*15;   
end
