% this program takes a matrix A and a perticular position k and returns the
% matrix after deleting the k-th row and k-th column.
function [A_jj] = DropColumn(A, k)
  d = size(A,2);
  if (k<1 || d < k) ;
      error('Error in function mypartition: partition location out of range');
  else 
     %A_jj = A;
     %A_jj(:,k) = []; % deleting k-th column
     
     %change for coder
     idx = zeros(size(A,2)-1,1);
     if k==1
         idx(:,1) = 2:size(A,2);
     elseif k==size(A,2)
         idx(:,1) = 1:(size(A,2)-1);
     else
         idx(:,1) = [1:(k-1) (k+1):(size(A,2))];
     end
     A_jj = A(:,idx);    
  end
end
