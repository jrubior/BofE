function P_matrices = findPermutationMatrices(L0init, info)
    % Initialize array to store permutation matrices
    P_matrices = cell(1, info.nshocks);
    
    % For each shock
    for j = 1:info.nshocks
        % Create identity matrix as starting point
        P_j = eye(size(L0init, 1));
        
        % Find which elements were changed in sign
        index = find(info.Ss{j,1} * L0init(:,j) < 0);
        if ~isempty(index)
            varindex = find((sum(info.Ss{j,1}(index,:),1) ~= 0) == 1);
            
            % For each variable that had sign change
            for i = 1:length(varindex)
                % Create a diagonal matrix with -1 at the position where sign changed
                P_j(varindex(i), varindex(i)) = -1;
            end
        end
        
        % Store this permutation matrix
        P_matrices{j} = P_j;
    end