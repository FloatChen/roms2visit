
clear

%addpath('../m_file/')
addpath('../../mytools/Visualization_tools/UTILITIES/m_map')

%f1='roms_400_838.nc';
%f1='../../D01/croco_run/monmean_clm_total.nc'
f1='/sat01/chenyong/CROCO02/mywork/CONFIGS/Run_BENGUELA_VHR_NBQ/CROCO_FILES/scs_400_838/croco_avg.nc'

f2='tt.cdl';

fout='tt1.nc';

flag_proj=1;

flag_mask=1;
cgrid2tri(f1,f2,fout,flag_proj,flag_mask)




