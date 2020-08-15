%% compute the distances btw each pair HMMs
%% here the distance using symtric KL

ND = length(HMMs{1});
N = ND*(ND-1)/2;
dist = zeros(ND,ND);
k = 1;
it =1;

hmms = HMMs{it};
data = Adata{it};

for i =1:ND
    hmm1 = hmms{i};
    data1 = data{i};
    for j =(i+1):ND
        
        hmm2 = hmms{j};
        data2 = data{j};
        
        dist(i,j) = 0.5*(vbhmm_kld(hmm1, hmm2, data1) + vbhmm_kld(hmm2, hmm1, data2));
        dist(j,i) = 0.5*(vbhmm_kld(hmm1, hmm2, data1) + vbhmm_kld(hmm2, hmm1, data2));
        pur(k) = dist(i,j);
        k = k+1;
    end
end


percent=2.0;
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);

position=round(N*percent/100);
sda=sort(pur(:));
dc=sda(position);

fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);


run cluster_dp

% 
[gamma_sorted,ordgamma]=sort(gamma,'descend');

tt=plot(1:ND,gamma_sorted(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
xlabel ('n')
ylabel ('\gamma')
%         
