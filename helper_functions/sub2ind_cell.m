function cellLinearIndex = sub2ind_cell(cellArr, cellInds)
k     = [1, cumprod(size(cellArr))];
cellLinearIndex = sum(k(1:length(cellInds)) .* (cellInds - 1)) + 1;
end