function check_compile_c_components

mock_data = randn(2,1000);
median_min_dist = 0.0221;
kernel_width =  median_min_dist*5;
optimizaiton_para = median_min_dist*1.5;
is_any_code_newly_compiled = 0;

fprintf('Checking whether c compoents in SPADE are properly compiled: \n')
try
    fprintf('   compute_min_dist_downsample.c ')
    [min_dist,NN_ind] = compute_min_dist_downsample(mock_data(:,1:10),mock_data);
    fprintf(' Done\n')
    fprintf('   compute_local_density.c ')
    [local_density] = compute_local_density(mock_data, kernel_width, optimizaiton_para); 
    fprintf(' Done\n')
    fprintf('   compute_min_dist_upsample.c ')
    [min_dist,NN_index] = compute_min_dist_upsample(mock_data, [1 -1 -1 1; 1 1 -1 -1]); 
    fprintf(' Done\n')
catch
    fprintf(' ERROR, not properly complied, compiling now ... ');
    mex compute_min_dist_downsample.c
    mex compute_local_density.c
    mex compute_min_dist_upsample.c
    mex faithful_downsample.c          
    mex L1_dist_one_to_many.c
    fprintf('Done\n')
    check_compile_c_components
end





