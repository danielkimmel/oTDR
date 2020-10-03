function dataTensor = data2Tensor(Data)
 [T, N] = size(Data(1).A);
  C = length(Data);
  dataTensor = permute(reshape(vertcat(Data.A)', N, T, C), [2 1 3]);



end