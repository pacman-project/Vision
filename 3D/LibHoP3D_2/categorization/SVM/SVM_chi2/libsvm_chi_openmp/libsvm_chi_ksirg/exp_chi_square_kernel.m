function kernel_matrix =chi_square_kernel(vectors,gamma)
%computes exponential chi-square kernel for svm , include sample serial number as first column

[nr_vec dim]=size(vectors);

%kernel_matrix=diag(ones(nr_vec,1));

kernel_matrix=zeros(nr_vec,nr_vec);
for i=1:nr_vec
    
    vec_i=vectors(i,:);
    for j=1:i
        val = exp_chi_square_val(vec_i,vectors(j,:),gamma);
        kernel_matrix(i,j)=val;
        kernel_matrix(j,i)=val;
    end
end

kernel_matrix=[(1:nr_vec)' , kernel_matrix];
end

function chi_val=exp_chi_square_val(vec1,vec2, gamma)
chi_val=(vec1-vec2)./(vec1+vec2+eps)*(vec1-vec2)';
chi_val = exp(-gamma*chi_val);
end
