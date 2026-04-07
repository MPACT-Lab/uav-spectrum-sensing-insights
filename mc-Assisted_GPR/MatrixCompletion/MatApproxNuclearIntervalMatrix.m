function [NewMAT, Error] = MatApproxNuclearIntervalMatrix(A, Init, Zone, lambda, N, std_dev_tol, display)
% This function performs matrix approximation when including only some
% entries. In other words, it solves the problem NewMAT = argmin
% ||(NewMAT-A).*Zone||_F, s.t. ||NewMAT||_*<=lambda. Because of convexity
% the solution is global/optimal.
%
% Coded by Gil Shabat, 2012
if nargin<7
    display=1;
end
NewMAT = Init;
S=sum2(Zone);
ii=0;
ErrDeviation = (NewMAT-A).*Zone;
ErrInterval = abs(ErrDeviation);
%Error = 0;
%disp('tolerance')
%disp(min(min(tol)))
if display
    fprintf('target nuclear norm: %d\n', lambda);
end
while (ii<=N) && ((any(any(ErrInterval > std_dev_tol))) || ii==0) % ErrInterval, tol both zero in non zone
    prevDeviation = ErrDeviation;
    % calculate new matrix
    NewMAT = NewMAT.*(1-Zone)+Zone.*A;
    %NewMAT = NewMAT.*(1-Zone) + Zone.* ((ErrDeviation>std_dev_tol).*(A+std_dev_tol/2)+(ErrDeviation<-std_dev_tol).*(A-std_dev_tol/2)+(ErrInterval<std_dev_tol).*NewMAT);
    NewMAT = FindNuclearNormApprox(NewMAT, lambda); % not garanteed to 
    if display
        fprintf('got: %d\n', sum(svd(NewMAT)));
    end
    ErrDeviation = (NewMAT-A).*Zone;
    ErrInterval = abs(ErrDeviation);
    if display
        fprintf('ii=%d, %g \n',ii,max2(ErrInterval-std_dev_tol)); 
    end
    Error = ErrInterval;
    ii=ii+1;
    if ii>2
        Diff = max2(abs(prevDeviation-ErrDeviation));
        if Diff<0.0001
            if display
                fprintf('Diff too small \n');
            end
            return;
        end
    end
    %fprintf('\n');
end

end

