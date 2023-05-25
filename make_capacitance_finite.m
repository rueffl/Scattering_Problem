function C = make_capacitance_finite(N,lij)
% MAKE_CAPACITANCE construct the capacitance matrix
%   N:      number of resonators
%   lij:    spacing between the resonators

    if N == 2
        D1 = cat(2,[1/lij(1)],[1/lij(1)]);
    elseif N > 2
        D1 = cat(2,[1/lij(1)],1./lij(1:end-1)+1./lij(2:end),[1/lij(end-1)]);
    end
    D2 = -1./lij;
    C = diag(D1) + diag(D2,1) + diag(D2,-1);

%     C = zeros(N);
%     for i = 1:N
%         for j = 1: N
%             if i == j - 1
%                 C(i, j) = C(i, j) - 1 / lij(i);
%             end
%             if i == j
%                 if j-1 <1
%                     C(i, j) = C(i, j) + (1 / lij(N) + 1 / lij(j));
%                 else
%                     C(i, j) = C(i, j) + (1 / lij(j - 1) + 1 / lij(j));
%                 end
%             end
%             if i == j + 1
%                 C(i, j) =  C(i, j) - 1 / lij(j);
%             end
%         end
%     end

end