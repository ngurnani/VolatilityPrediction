function[modelDensity]=gpVolObsModel(data,v,m,isT,sigTransform)

Sigma=shiftdim(sigTransform(v),-2);

if ~isT %y_t ~ N(0,Sigma_t)
    modelDensity=mvnpdfln(data,[],Sigma,0); %w_t=w_{t-1}*p(y_t|mu_t,m_t-1)
else %y_t ~ Student(nu,0,S_t)
    %Sigma_t=nu/(nu-2) S_t
    [d d N] =size(Sigma);    
    modelDensity=zeros(N,1);
    nu=exp(m(:,end))+2; %log(nu-2)~N(.,.)
    for i=1:N
        covAdj= (nu(i)-2)./nu(i);
        modelDensity(i,:)=mvtpdfln(data,Sigma(:,:,i).*covAdj,nu(i),[]);
    end
end