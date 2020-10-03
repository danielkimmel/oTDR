
function [f, gradQi, Rscale, Cscale, Di] = oTDRobjFn(Q, R, codedParams, trCountmtx)
    numNeus = size(R{1},2);
    numConds = size(trCountmtx,2);
    Di = zeros(length(R),1);
    for i = 1:length(R)
    Qi = Q(:, i);
    numTimes = size(R{i},1)./numConds;
    Ci = kron(codedParams{i}, ones(numTimes,1));
    Rscale = R{i}.*kron(sqrt(trCountmtx'),ones(numTimes,1));
    Cscale = repmat(Ci,1, numNeus).*kron(sqrt(trCountmtx'),ones(numTimes,1));
    A = bsxfun(@times, Cscale, Qi');
    D = (Rscale(:)'*A(:))./norm(A(:))^2;
    Di(i) = D;
    %%%%%% solve the offset component
    Rscale0 = Rscale - D*A;
    CbarScale0 = repmat(kron(ones(numConds, 1), ones(numTimes,1)),1, numNeus).*kron(sqrt(trCountmtx'),ones(numTimes,1));
    
    Q0 = diag(Rscale0'*CbarScale0)./diag(CbarScale0'*CbarScale0);
    Abar0 = bsxfun(@times, CbarScale0, Q0');
    RscaleNew = Rscale-Abar0;
    gradQ = 2*(-diag(RscaleNew'*Cscale)*D+Qi.*D.*diag(Cscale'*Cscale).*D);
    gradQi(:,i) = gradQ./norm(Rscale,'fro')^2;
    fi(i) = norm(RscaleNew-A*D,'fro')^2./norm(Rscale,'fro')^2;
%     fi(i) = TDRobjectiveFn(permute(R{i}, [3, 2, 1]), [codedParams{i} ones(numConds,1)], permute([D*Qi, Q0], [3 1 2]), trCountmtx);
    end
    f = mean(fi);
end


% %% objective ignore gradients
% function [f, gradQi, Rscale, Cscale, Di] = oTDRobjFn(Q, R, codedParams, trCountmtx)
%     numSignals = length(codedParams);
%     numNeus = size(R{1},2);
%     numConds = size(trCountmtx,2);
%     signalsIx = 1:numSignals;
%     for i = signalsIx
%     Qi = Q(:, i);
%     numTimes = size(R{i},1)./numConds;
%     Ci = kron(codedParams{i}, ones(numTimes,1));
%     Rscale = R{i}.*kron(sqrt(trCountmtx'),ones(numTimes,1));
%     Cscale = repmat(Ci,1, numNeus).*kron(sqrt(trCountmtx'),ones(numTimes,1));
%     A = bsxfun(@times, Cscale, Qi');
%     D = (Rscale(:)'*A(:))./norm(A(:))^2;
%     Di(i) = D;
%     niSet = signalsIx(signalsIx~=i);
%     z = 1;
%     CbarScale = [];
%     Abar = [];
%     Dbar = [];
%     AbarD = [];
%     for j = niSet
%         Qbar = Q(:, j);
%         Cj = kron(codedParams{j}, ones(numTimes,1));
%         CbarScale(:,:,z) = repmat(Cj,1, numNeus).*kron(sqrt(trCountmtx'),ones(numTimes,1));
%         Abar(:,:,z) = bsxfun(@times, CbarScale(:,:,z), Qbar');
%         Dbar(z) = (Rscale(:)'*reshape(Abar(:,:,z), [],1))./norm(reshape(Abar(:,:,z), [],1))^2;
%         AbarD(:,:,z) = Abar(:,:,z)*(Dbar(z));
%         
%         z = z+1;
%     end
%     %%%%%% solve the offset component
%     Rscale0 = Rscale-sum(AbarD, 3)-D*A;
%     CbarScale0 = repmat(kron(ones(numConds, 1), ones(numTimes,1)),1, numNeus).*kron(sqrt(trCountmtx'),ones(numTimes,1));
%     
%     Q0 = diag(Rscale0'*CbarScale0)./diag(CbarScale0'*CbarScale0);
%     Abar0 = bsxfun(@times, CbarScale0, Q0');
% 
%     
% %%%%%%%%%% just a sanity check that Qdc in the compact form is equal to Q0
% % % % % %     V = Q*D;
% % % % % %     Vbar1 = Qbar(:,1)*Dbar1;
% % % % % %     Vbar2 = Qbar(:,2)*Dbar2;
% % % % % %     for i = 1:numNeus
% % % % % %         r = R(:,i);
% % % % % %         CC = [ones(size(C)) C  Cbar];
% % % % % %         condMtx = diag(trCountmtx(i,:).^0.5);
% % % % % %         r = condMtx*r;
% % % % % %         CC = condMtx*CC;
% % % % % %         rnew = r- CC(:,2)*V(i)- CC(:,3)*Vbar1(i)- CC(:,4)*Vbar2(i);
% % % % % %         Q0(i,1) = rnew'*CC(:,1)./(CC(:,1)'*CC(:,1));
% % % % % %     end
% % % % % %     norm(Qdc-Q0)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     RscaleNew = Rscale-sum(AbarD, 3)-Abar0;
%     gradQ = 2*(-diag(RscaleNew'*Cscale)*D+Qi.*D.*diag(Cscale'*Cscale).*D);
%     gradQi(:,i) = gradQ./norm(Rscale,'fro')^2;
%     fi(i) = norm(RscaleNew-A*D,'fro')^2./norm(Rscale,'fro')^2;
%     end
%     f = mean(fi);
% %     f = sum(f);
%     %% loop over this just as a sanity check for the compact form above
% % % % % %     for i = 1:numNeus
% % % % % % %         rs = Rscale(:,i);
% % % % % % %         Ps = [Cscale(:,i) CbarScale(:,i,1) CbarScale(:,i,2)]*diag([D Dbar1 Dbar2]);
% % % % % % %         b = [Q(i) Qbar(i,:)]';
% % % % % % %         gradQ_i(i,:) = 2*(Ps'*Ps*b-Ps'*rs);
% % % % % %         %%%%%%%%%%%%%% below is equivalent to above (look up notes page 10)
% % % % % %         rsNew(:,i) = Rscale(:,i)- CbarScale(:,i,1)*Dbar1*Qbar(i,1)- CbarScale(:,i,2)*Dbar2*Qbar(i,2);
% % % % % %         PPs(:,i) = Cscale(:,i)*D;
% % % % % %         gradQ_i(i,1) = 2*(PPs(:,i)'*PPs(:,i)*Q(i)-PPs(:,i)'*rsNew(:,i));
% % % % % % 
% % % % % %     end
% % % % % %     norm(2*(diag(PPs'*PPs).*Q-diag(PPs'*rsNew))-gradQ_i) % should be zero
% % % % % %     gradQ_i = gradQ_i(:,1)./norm(RscaleNew,'fro')^2;
% % % % % %     norm(gradQ - gradQ_i) % should be zero
% 
% end