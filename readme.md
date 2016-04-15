memCov.m can be used in exactly the same way as the inbuilt MatLab function. Firstly, it will try to compute the covariance matrix of X using the inbuilt MatLab function cov.m.

If the process runs out of memory, memCov.m will compute the covariance matrix iteratively, computing small chunks of X sequentially and summing the result.

The method is summarised in the attached notes, 'Combining_covariance_matrices.pdf'.