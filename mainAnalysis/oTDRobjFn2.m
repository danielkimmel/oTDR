function [f, gradBetai] = oTDRobjFn2(Beta, R, codedParams, condRep)

    numSignals = length(codedParams);
    numNeus = size(R,2);
    numConds = size(condRep,2);
    signalsIx = 1:numSignals;
    numTimes = size(R,1)./numConds;
    Rscale = R.*kron(sqrt(condRep'),ones(numTimes,1));
    for i = signalsIx
    Betai = Beta(:, i);
    Ci = kron(codedParams{i}, ones(numTimes,1));
    Cscale = repmat(Ci,1, numNeus).*kron(sqrt(condRep'),ones(numTimes,1));
    A = bsxfun(@times, Cscale, Betai');
    niSet = signalsIx(signalsIx~=i);
    z = 1;
    CbarScale = [];
    Abar = [];
    for j = niSet
        Betaj = Beta(:, j);
        Cj = kron(codedParams{j}, ones(numTimes,1));
        CbarScale(:,:,z) = repmat(Cj,1, numNeus).*kron(sqrt(condRep'),ones(numTimes,1));
        Abar(:,:,z) = bsxfun(@times, CbarScale(:,:,z), Betaj');        
        z = z+1;
    end    
    RscaleNew = Rscale-sum(Abar, 3);
    gradBetai(:,i) = 2*(-diag(RscaleNew'*Cscale)+Betai.*diag(Cscale'*Cscale))./norm(Rscale,'fro')^2;
    fi(i) = norm(RscaleNew-A,'fro')^2./norm(Rscale,'fro')^2;
    end
    f = mean(fi);
end