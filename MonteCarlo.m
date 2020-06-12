% MonteCarlo Generate Required Data
function result = MonteCarlo(M, T, S)
  % create a three dimessional matrix with 2 rows T+S+1 columns and M pages and filled with zero
  result = zeros(2, T+S+1, M);
  % traverse pages
  for i = 1: M
    % save data in a result page
    result(:, :, i) = GenerateSeriesData(T+S+1);
  end
end

% GenerateSeriesData Using DGP to Generate Series Data by a Given Series Length
function data = GenerateSeriesData(length)
  for j = 1: length
    % DGP, generate random series data with length T+S+1
    if j == 1
      data = [[0, 0]', ];
    else
      SigmaU = [1, 0.3; 0.3, 1];
      % generate u_t while using, instead of generate a series of u_t in advance
      data(:, j) = [0.9, 0; 0.5, 0.5] * data(:, j-1) + mvnrnd([0, 0], SigmaU)';
    end
  end
end
