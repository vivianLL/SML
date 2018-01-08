clear;
clc;
warning('off'); %不在控制台显示警告
dimension=30;
filename=dir('pictures/*.jpeg');%读取文件
pic_amount=length(filename);
pi=zeros(8,pic_amount);
mu=zeros(dimension,8,pic_amount);
sigma=zeros(dimension,dimension,8,pic_amount);

for i=1:pic_amount
    fprintf('正在处理第%d张图片...\n',i);
    data=getDCT(dimension,['pictures/',filename(i).name]);
    [estimate, varargout] = gmm(data) ;
    pi(:,i)=estimate.weight;%8*1
    mu(:,:,i)=estimate.mu;%30*8
    sigma(:,:,:,i)=estimate.sigma;%30*30*8

end
%process all pictures,8*N个高斯分量聚类为64个分量
pi_new=zeros(64,1);%64*1
mu_new=zeros(dimension,64);%30*64
sigma_new=zeros(dimension,dimension,64);%30*30*64
%初始化
for mm=1:64
    jj=floor((mm-1)/8)+1;
    kk=mm-(jj-1)*8;
    pi_new(mm)=pi(kk,jj)/8;
    mu_new(:,mm)=mu(:,kk,jj);
    sigma_new(:,:,mm)=sigma(:,:,kk,jj);
end
 h=zeros(pic_amount,8,64);
 median_11=zeros(pic_amount,8,64);
 omega=zeros(pic_amount,8,64);
 median_12=zeros(dimension,pic_amount,8,64);
 median_13=zeros(dimension,dimension,pic_amount,8,64);
 

fprintf('正在进行聚类...\n');
for tt=1:20
% E步
 for j=1:pic_amount
     for k=1:8
          for m=1:64
             median=gauss(mu(:,k,j),mu_new(:,m),sigma_new(:,:,m))*exp(trace(inv(sigma_new(:,:,m))*sigma(:,:,k,j))*(-0.5)); 
             h(j,k,m)=pi_new(m,1)*(median.^(pi(k,j)));
          end
          h_sum=sum(h(j,k,:));
          for m=1:64
              h(j,k,m)=h(j,k,m)/h_sum;
          end
     end  
 end
 
 %这是M步
%  median_11=zeros(pic_amount,8,64);
%  omega=zeros(pic_amount,8,64);
h_m=sum(sum(h,1),2);
for m=1:64 
    for j=1:pic_amount
         for k=1:8
            median_11(j,k,m)=pi(k,j)*h(j,k,m);
         end
    end
    pi_new(m,1)=h_m(1,1,m)/(pic_amount*8); %将pi_new的值更新
end
median_11_sum=sum(sum(median_11,1),2);%1*1*64
 for j=1:pic_amount
     for k=1:8
          for m=1:64
              omega(j,k,m)=median_11(j,k,m)/median_11_sum(1,1,m);
          end
     end
 end
% median_12=zeros(dimension,pic_amount,8,64);
% median_13=zeros(dimension,dimension,pic_amount,8,64);
for m=1:64 
    for j=1:pic_amount
         for k=1:8
             median_12(:,j,k,m)=omega(j,k,m)*mu(:,k,j);
             median_13(:,:,j,k,m)=omega(j,k,m)*(sigma(:,:,k,j)+(mu(:,k,j)-mu_new(:,m))*((mu(:,k,j)-mu_new(:,m))'));
         end
    end
end
median_12_sum=sum(sum(median_12,2),3);
median_13_sum=sum(sum(median_13,3),4);
mu_new_old=mu_new;
for m=1:64 
    mu_new(:,m)=median_12_sum(:,1,1,m);%将mu_new更新
    sigma_new(:,:,m)=median_13_sum(:,:,1,1,m);
end

% if(sum(sum(abs(mu_new./mu_new_old -1)))/(dimension*64) <1e-2)
%             break; 
% end

end
fprintf('聚类完毕\n');
model=[];
model.Pi = pi_new;
model.Mu = mu_new;
model.Sigma = sigma_new;
save model.mat model;



    
