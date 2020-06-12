% Generate Parameter Series
function Parameters = GenerateParameterSeries(Y, Z, residual, n)
  Parameters = [];
  % Repeate N times
  for i = 1:n
    Y_New = Y - Substitute(residual);
    Parameters(:, :, end+1) = Y_New * Z' * inv(Z*Z');
    end
  end


% Substitute Algo
% Random column Substitute
function residual = Substitute(residual)
  residualStar = [];
  for i = residual
    residualStar(:, end+1)  = i -  mean(residual, 2);
  end
  % substitute random times ([1, length(residual)])
  for i = 1:randi(length(residual))
  % substituted column index is a uniform distribution,
  % intended to substitute each column equally
    SubstitutedIndex = round(unifrnd(1, length(residual)));
    residual(:, SubstitutedIndex) = residualStar(:, SubstitutedIndex);
  end
end

