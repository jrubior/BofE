function X = ensure3D(X)
    if iscell(X) && ~isempty(X)
        X = cat(3, X{:});
    end
end