function Data = tensor2Data(dataTensor, otherDataVars)
[T, N, C] = size(dataTensor);
XN = reshape(permute(dataTensor,[1 3 2]), [], N);
if length(otherDataVars)<C
   otherDataVars = repmat(otherDataVars(1), C,1);
end
for c = 1:C
    Data(c) = otherDataVars(c);
end
for c = 1:C
    mskConds = false(1, C);
    mskConds(c) = true;
    mskConds = repmat(mskConds, T, 1);
    mskConds = mskConds(:);
    Data(c).A = XN(mskConds,:);
end
end