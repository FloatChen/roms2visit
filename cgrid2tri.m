function  []=cgrid2tri(f1,f2,fout,flag_proj,flag_mask)


% Write ROMS model output to fvcom format for Visualization with visit;
% Divied a C grid into two triangle 
% Transform ROMS output to FVCOM output 
% roms mesh to  triangle mesh for visit 
% Input: ROMS outfile(f1); cdl file(f2); FVCOM form file(fout); proj option(flag_proj)
% by chenyong 2022-2-19 

%clear
%f1='model_mean/avg_2008.nc';
 
lonr=ncread(f1,'lon_rho');  
latr=ncread(f1,'lat_rho');  
h=ncread(f1,'h');           
maskr=ncread(f1,'mask_rho');
h(maskr==0)=0;

try
  % ROMS; COAWST;
  Time=double(ncread(f1,'ocean_time'))/86400+datenum(2008,1,1);
catch
  %Agrif-roms; CROCO
  Time=double(ncread(f1,'scrum_time'))/86400+datenum(1990,1,1);
end

[n1,n2]=size(lonr);
if flag_mask
    xp=lonr(maskr==1);
    yp=latr(maskr==1); 
else
    xp=lonr(:);
    yp=latr(:);
end

X=0*xp;Y=0*yp;
if flag_proj
set_lon=[min(lonr(:))-1,max(lonr(:))+1];
set_lat=[min(latr(:))-1,max(latr(:))+1];
 m_proj('mercator','lat',set_lat,'lon',set_lon);
 [X,Y]=m_ll2xy(xp,yp);
X=X*1e5; Y=Y*1e5;
end

% Divied to two triangle and save cell number
if ~flag_mask
  nv=nan((n1-1)*(n2-1)*2,3);
end
cnt=0
for j=1:n2-1
  for i=1:n1-1
    if flag_mask
     if maskr(i,j)+maskr(i+1,j)+maskr(i,j+1)+maskr(i+1,j+1)==4
      cnt=cnt+1;
      nv(cnt,:)=[i+(j-1)*n1,i+1+(j-1)*n1,i+j*n1];
      cnt=cnt+1;
      nv(cnt,:)=[i+1+(j-1)*n1,i+1+j*n1,i+j*n1];
     end
    else
     cnt=cnt+1;
     nv(cnt,:)=[i+(j-1)*n1,i+1+(j-1)*n1,i+j*n1];
     cnt=cnt+1;
     nv(cnt,:)=[i+1+(j-1)*n1,i+1+j*n1,i+j*n1];
    end
  end
end

if flag_mask
ind=find(maskr==1);
for i=1:size(nv,1);
    ii=find(nv(i,1)==ind);nv(i,1)=ii;
    ii=find(nv(i,2)==ind);nv(i,2)=ii;
    ii=find(nv(i,3)==ind);nv(i,3)=ii;
end
end

%triplot(nv,xp,yp);axis([133,134,0.8,2])

%trisurf(nv,xp,yp,h(:));view(2);shading flat 

%-------
  xc=(xp(nv(:,1))+xp(nv(:,2))+xp(nv(:,3)));
  yc=(yp(nv(:,1))+yp(nv(:,2))+yp(nv(:,3)));

  XC=(X(nv(:,1))+X(nv(:,2))+X(nv(:,3)));
  YC=(Y(nv(:,1))+Y(nv(:,2))+Y(nv(:,3)));

sc_r=ncread(f1,'s_rho');
Cs_r=ncread(f1,'Cs_r');
sc_w=ncread(f1,'s_w');
Cs_w=ncread(f1,'Cs_w');
hc=ncread(f1,'hc');
N=length(sc_r);


if flag_mask
pp=ind;
else
pp=1:size(h,1)*size(h,2); pp=pp';
end

z0p=nan(length(pp),N);
z1p=nan(length(pp),N+1);
tmp=[];
         for k=1:N
                tmp=(hc.*sc_r(k)+Cs_r(k).*h)./(hc+h);
              %  z0(:,:,k)=tmp;
               % z(:,:,k)=zeta+(zeta+h).*z0;
              z0p(:,N-k+1)=tmp(pp);
         end

         for k=1:N+1
              tmp=(hc.*sc_w(k)+Cs_w(k).*h)./(hc+h);
           %    z1(:,:,k)=tmp;
               % z(:,:,k)=zeta+(zeta+h).*z0;
               z1p(:,N+1-k+1)=tmp(pp);
         end

% fvcom ncfile generate
%fout='tt1.nc';
status=system(['cp   ', f2, ' tmp.cdl'])
status=system(sprintf('sed -i ''s/node = xxxx/node = %d/''  tmp.cdl',length(xp) ) )
status=system(sprintf('sed -i ''s#nele = xxxx#nele = %d#''  tmp.cdl',length(xc) ))
%flag_proj=1;
if flag_proj
  status=system(sprintf('sed -i ''s#GeoReferenced#%s#'' tmp.cdl', 'Cartesian' ))
else
  status=system(sprintf('sed -i ''s#Cartesian#%s#'' tmp.cdl', 'GeoReferenced' ))
end
status=system(['ncgen -b tmp.cdl  -o '  fout])
status=system('rm tmp.cdl')

 ncwrite(fout,'lon',xp);
 ncwrite(fout,'lat',yp);
 ncwrite(fout,'nv',nv);
 ncwrite(fout,'lonc',xc);
 ncwrite(fout,'latc',yc);
 ncwrite(fout,'x',X);
 ncwrite(fout,'y',Y);
 ncwrite(fout,'xc',XC);
 ncwrite(fout,'yc',YC);
 ncwrite(fout,'h',h(pp));
 ncwrite(fout,'siglev',z1p);
 ncwrite(fout,'siglay',z0p);

%----------------------------
for t=1:length(Time)
  t
  ncwrite(fout,'time',Time(t)-datenum(1858,11,17),[t]);
  ncwrite(fout,'Times',datestr(Time(t),'yyyy-mm-ddTHH:MM:SS.000000')',[1,t])

  zeta=ncread(f1,'zeta',[1,1,t],[inf,inf,1]);
  ncwrite(fout,'zeta',zeta(pp),[1 t])

%  u=ncread(f1,'u_eastward',[1,1,1,t],[inf,inf,inf,1]);
%  v=ncread(f1,'v_eastward',[1,1,1,t],[inf,inf,inf,1]);
  temp=ncread(f1,'temp',[1,1,1,t],[inf,inf,inf,1]);
  salt=ncread(f1,'salt',[1,1,1,t],[inf,inf,inf,1]);
  
  tmp=nan(length(pp),N);
  cnt=0;
  for i=N:-1:1
     cnt=cnt+1;
     tmpp=temp(:,:,i);
     tmp(:,cnt)=tmpp(pp);
  end
  ncwrite(fout,'temp',tmp,[1,1,t]);
 
  tmp=nan(length(pp),N);
  cnt=0;
  for i=N:-1:1
     cnt=cnt+1;
     tmpp=salt(:,:,i);
     tmp(:,cnt)=tmpp(pp);
  end
  ncwrite(fout,'salinity',tmp,[1,1,t]);  
end


return




