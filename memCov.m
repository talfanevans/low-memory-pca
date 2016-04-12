%% Compute incrmental covariance matrix to avoid memory overloading

function C_tot =  memCov(X)

% First of all try using the regular method. If this doesn't work try the
% iterative method
try C_tot = cov(X);
    
    % Try the iterative method
catch %#ok<CTCH>
    
    % [Observations,Variables]
    [O,V] = size(X);
    
    % The algorithm starts with as few segments as possible (2). If it
    % encounters a memory error, it increases the number of segments and
    % tries again until it succeeds.
    num_segments = 4;
    out_flag = 0; % This is set to 1 when the algorithm has succesfully completed
    
    while out_flag == 0
        num_segments = num_segments*2;
        
        try
            % If intended number of segments is fewer than the number of observations,
            % igore the iterative algorithm and just do it all in one go
            if O < num_segments
                num_segments = 1;
            end
            
            clc
            fprintf('Trying %i segments\n...',num_segments)
            
            % Indices of each of the segments
            seg_index = [1:floor(O/num_segments):floor(O/num_segments)*num_segments,O+1];
            
            C_tot = zeros(V);
            % Iteratively combine each new segment to the running total
            for i = 1:num_segments
                %fprintf('\b\b\b%i%%',round(100*i/num_segments))
                clc
                fprintf('%i%%',round(100*i/num_segments))
                
                % Number of elements in the new segment and the running total respectively
                n_seg = seg_index(i+1) - seg_index(i);
                n_tot = seg_index(i)-1;
                
                % Calculate the covariance matrix and means vector of the new segment
                C_seg =   cov(X(seg_index(i):(seg_index(i+1)-1),:));
                Mu_seg = mean(X(seg_index(i):(seg_index(i+1)-1),:));
                
                if i == 1
                    % On the first iteration the covariance matrix and the mean vector
                    % are the same as the covariance and mean for the first segment of
                    % the total data matrix X
                    C_tot = C_seg;
                    Mu_tot = Mu_seg;
                else
                    % Update the new covariance matrix using the running total and new segment
                    M = Mu_tot - Mu_seg;
                    [MA,MB] = meshgrid(M,M);
                    C_tot =      C_tot*     (n_tot-1)/(n_seg+n_tot-1) + C_seg*     (n_seg-1)/(n_seg+n_tot-1) + (MA.*MB)*                                   (n_tot*n_seg)/((n_seg+n_tot)*(n_seg+n_tot-1));
                    % Option for calculating element-wise is much slower
                    %{
                    for j = 1:V
                for k = 1:V
                    C_tot(j,k) = C_tot(j,k)*(n_tot-1)/(n_seg+n_tot-1) + C_seg(j,k)*(n_seg-1)/(n_seg+n_tot-1) + (Mu_tot(j)-Mu_seg(j))*(Mu_tot(k)-Mu_seg(k))*(n_tot*n_seg)/((n_seg+n_tot)*(n_seg+n_tot-1));
                end
            end
                    %}
                    
                    % Update the new means vector using the running total and new segment
                    Mu_tot = (Mu_tot*n_tot + Mu_seg*n_seg)/(n_tot+n_seg);
                    
                    clc
                    fprintf('n_seg = %.2d\nn_tot = %.2d\nN = %i\n',n_seg,n_tot,O);
                    
                end
                
            end
            
            out_flag = 1;            
            
        end
        
    end
    
end
