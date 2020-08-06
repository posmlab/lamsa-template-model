function [axes,pac] = sensitive_axes(fun,x0,grad_pct,mesh_pct,pca_pct,pca_N)
%Determines the most sensitive axis by the gradient, then performs a PCA on
%the projection of the function onto the tangent plane
    if nargin == 4 %Set default value for pca_N
        pca_N = 5;
    end
    
    %% Find gradient and find tangent plane
    %Evaluate the gradient
    f0 = feval(fun,x0);
    grad = relative_gradient(fun,x0,grad_pct)
    %If the gradient is small, call it zero
    if 2*grad_pct*grad/f0<1e-4 %If the total variation is less than some small percentage
        grad = zeros(length(x0),1);
    end
    %Determine the matrix U that maps to the tangent plane
    nullBasis = null(grad')
    
    %% Create mesh in input space
    dim = size(nullBasis,2);
    %build axis vectors
    sz = pca_N*ones(1,dim);
    perf = zeros(sz);
    
    %Create axes in input space
    input_ax_vals = zeros(length(x0),pca_N);
    for i = 1:length(x0)
        low_bound = -mesh_pct;
        up_bound = +mesh_pct;
        input_ax_vals(i,:) = linspace(low_bound,up_bound,pca_N);
    end
    
    %Assemble a string to plug into eval(). We put the columns of input_mesh as the inputs to ndgrid.
    eval_str = ['[p{' num2str(1) '}'];
    for i = 2:length(x0)
        eval_str = [eval_str ',p{' num2str(i) '}'];
    end
    eval_str = [eval_str '] = ndgrid(input_ax_vals(1,:)'];
    
    for i = 2:length(x0)
        eval_str = [eval_str ',input_ax_vals(' num2str(i) ',:)'];
    end
    eval_str = [eval_str ');'];
    
    eval(eval_str); %Run the constructed string
    
    %Create an array with its columns being gridpoints in column vector form
    input_mesh_pts = zeros(length(x0),0);
    for i = 1:numel(p{1})
        for j = 1:length(x0)
            input_mesh_pts(j,i) = p{j}(i);
        end
    end

    %% Create mesh in nullspace
    null_proj_pts = zeros(dim,numel(p{1}));
    %Find min and max values on each nullspace axis
    for i = 1:dim
        for j = 1:numel(p{1})
            null_proj_pts(i,j) = dot(input_mesh_pts(:,j),nullBasis(:,i));
        end
    end
    min_vals = min(null_proj_pts,[],2); %Column vector of the minimum value of each nullspace dimension
    max_vals = max(null_proj_pts,[],2); %Maximum values
    
    %Create axes in nullspace
    for i = 1:dim
        null_ax_vals(i,:) = linspace(min_vals(i),max_vals(i),pca_N);
    end
    
    %Assemble a string to plug into eval(). We put the columns of null_ax_vals as the inputs to ndgrid.
    null_str = ['[q{' num2str(1) '}'];
    for i = 2:dim
        null_str = [null_str ',q{' num2str(i) '}'];
    end
    null_str = [null_str '] = ndgrid(null_ax_vals(1,:)'];
    
    for i = 2:dim
        null_str = [null_str ',null_ax_vals(' num2str(i) ',:)'];
    end
    null_str = [null_str ');'];
    
    eval(null_str); %Run the constructed string
    
    %Create an array with its columns being gridpoints in column vector form
    null_mesh_pts = zeros(dim,0);
    for i = 1:numel(q{1})
        for j = 1:dim
            null_mesh_pts(j,i) = q{j}(i);
        end
    end
    
    %% Evaluate function on the mesh and select coordinates for pca
    num_pca_points = 0;
    crds_to_pca = zeros(dim,0);
    for i = 1:numel(q{1})
        progress = i/(numel(q{1}));
          F = feval(fun,x0.*(1+nullBasis*null_mesh_pts(:,i))')
          f0
          pct_diff = abs(F/f0-1)
        if pct_diff<pca_pct
            num_pca_points = num_pca_points + 1;
            crds_to_pca(:,num_pca_points) = null_mesh_pts(:,i);
        end
    end
    %Now perform PCA on the resulting points
    s = size(crds_to_pca);
    if s(2) == 0
        error('pca_percent too low! No points to perform pca on.');
    end
    [pca_u_basis,~,~,~,explained] = pca(crds_to_pca');
    pca_x_basis = nullBasis*pca_u_basis;
    axes = [grad/norm(grad) flip(pca_x_basis,2)];
    coeff=pca_u_basis.^2;
    explained=explained*.01;
    puck=coeff.*explained';
    percentcoeff=nullBasis*puck;
    buck=[0];
    for i=1:length(x0)
        buck_new=sum(percentcoeff(i,:));
        buck=[buck,buck_new];
    end
    buck=buck(2:end);
    pac=(1-buck.^2)/sum(1-buck.^2)*100;
end