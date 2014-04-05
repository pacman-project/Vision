function kernel_matrix =chi_square_kernel(vectors)
%computes chi-square kernel for svm

[nr_vec dim]=size(vectors);

%kernel_matrix=diag(ones(nr_vec,1));

kernel_matrix=zeros(nr_vec,nr_vec);

if( issparse(vectors))
    vectors= vectors';
    for i=1:nr_vec
        
        vec_i=vectors(:,i);
        for j=1:i
            val = chi_square_val(vec_i',vectors(:,j)');
            kernel_matrix(i,j)=val;
            kernel_matrix(j,i)=val;
        end
    end
else
    
    for i=1:nr_vec
        
        vec_i=vectors(i,:);
        for j=1:i
            val = chi_square_val(vec_i,vectors(j,:));
            kernel_matrix(i,j)=val;
            kernel_matrix(j,i)=val;
        end
    end
end
kernel_matrix=[(1:nr_vec)' , kernel_matrix];
end

function chi_val=chi_square_val(vec1,vec2)
% 1-sum( (xi-xj)^2/0.5*(xi+xj) )

chi_val=2*(vec1-vec2)./(vec1+vec2+eps)*(vec1-vec2)';

chi_val=1-chi_val;
%chi_val=exp(chi_val);
% a=(vec1-vec2).^2;
% b=(vec1+vec2);
% resHelper = 1-sum(2*a./(b + eps));
end
