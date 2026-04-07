function [CompletedMat,status] = complete_matrix_paper(mat,B,n_iter,matrix_change_tolerance,norm_stability_threshold,display,method)
    % status = 0 (no error, matrix_change_tolerance acheived), 1 matrix
    % changed more than the tolerance
    if method=="frobenius"
        [CompletedMat, status] = MatrixCompletion(mat, B,n_iter, 'spectral', norm_stability_threshold, matrix_change_tolerance, display);
    elseif method=="interval"
        [CompletedMat, status] = MatrixCompletionIntervalConstraint(mat, B,n_iter, 'nuclear', norm_stability_threshold, matrix_change_tolerance, display, B.*mat_variance*1.0);
    elseif method=="intervalMatrix"
        [CompletedMat, status] = MatrixCompletionIntervalMatrix(mat, B,n_iter, 'nuclear', norm_stability_threshold, B.*sqrt(mat_variance)*matrix_change_tolerance, display);
    elseif method=="weighted"
        [CompletedMat, status] = MatrixCompletionErrorWeighted(mat,B,n_iter,'nuclear', norm_stability_threshold, matrix_change_tolerance, display, B.*mat_weight);
    else
        error("unknown matrix completion. valid: frobenius, interval, weighted")
    end
    if(display)
        fprintf('\n Corrupted matrix nuclear norm (initial): %g \n',sum(svd(mat.*B)));
        fprintf('Restored matrix nuclear norm (final): %g \n',sum(svd(CompletedMat)));
        Diff_sq = (CompletedMat - mat).^2;
        fprintf('MSE on known entries: %g \n',sqrt(sum2(Diff_sq.*B)/sum(B(:)) ));
    end
end

