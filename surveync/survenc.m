% Needs m_map toolbox to run
clear
clc
close all
%% Read the transect start/end points
filename = 'D:\repos\Github\PELFOSS\surveync\transect_endpoints_IESNS2017.txt';
pos = surveycn_transectendpoints(filename);

%% Set survey speed 
s  = 10;%nmi per h*1800/3600;% m/s
t0 = datenum(2017,05,01);
t = .1; % trawl stations per nmi
%% Set start position!
start(1).time = t0;
start(1).lat = 60.398206; 
start(1).lon = 5.314336;
start(1).vessel = 'NO';

start(2).time = t0;
start(2).lat = 69;  
start(2).lon = 30.06;
start(2).vessel = 'RU';

start(3).time = t0;
start(3).lat = 60.398206; 
start(3).lon = 5.314336;
start(3).vessel = 'IS';

start(4).time = t0;
start(4).lat = 60.398206; 
start(4).lon = 5.314336;
start(4).vessel = 'EU';

start(5).time = t0;
start(5).lat = 60.398206; 
start(5).lon = 5.314336;
start(5).vessel = 'FO';

vessels= {'NO','RU','IS','EU','FO'};

%% Calcualate vessel start stop positions of transects and mile positions

for i=1:length(pos)
    lats  = interp1([0 1],[pos(i).latstart pos(i).latstop],0:.001:1);
    longs = interp1([0 1],[pos(i).lonstart pos(i).lonstop],0:.001:1);
    r = m_lldist(longs,lats);%km
    % We would like spacing of 1nmi = 1,852km, and the central part of the transect
    % as the position for the model
    pos(i).R = [0; cumsum(r)/1.852];
    pos(i).lat = interp1(pos(i).R,lats ,1:ceil(max(pos(i).R)),'linear','extrap');
    pos(i).lon = interp1(pos(i).R,longs,1:ceil(max(pos(i).R)),'linear','extrap');
    pos(i).dist = 1:ceil(max(pos(i).R));
end

%% Sort and order the transects by vessel
pos2=struct;
for i=1:length(vessels)
    k=1;
    % Create subset of trasects to one vessel
    for j=1:length(pos)
        if strcmp(pos(j).vessel,vessels{i})
            pos2(i).pos(k)=pos(j);
            pos2(i).vessel=pos(j).vessel;
            pos2(i).startpos=start(i);
            k=k+1;
        end
    end
    pz = pos2(i);
    pos2(i).transectorder = survey_transectorder(pz);
end

%% Calculate time (between transects)
LAT =[];
LON =[];
TIME=[];
for p=1:length(pos2)
    Dum=survey_transeectdistance(pos2(p),s);
    LAT=[LAT Dum.LAT];
    LON=[LON Dum.LON];
    TIME=[TIME Dum.TIME];
end

%% Plot the results
clf
%m_proj('miller','long',[-20 45],'lat',[60 75])
m_proj('lambert','long',[-20 45],'lat',[60 80])
m_coast('patch',[1 .85 .7]);
m_elev('contourf',[500:500:6000]);
m_grid('box','fancy','tickdir','in');
colormap(flipud(copper));
hold on

% Plot transects
for i=1:length(pos2)
    disp(pos2(i).vessel)
    switch pos2(i).vessel
        case 'EU'
            vz = 'b';
        case 'RU'
            vz = 'r';
        case 'NO'
            vz = 'k';
        case 'IS'
            vz = 'g';
        case 'FO'
            vz = 'c';
    end
    
    for j=1:(length(pos2(i).transectorder)-1)
        % First trackposition
        tr=abs(pos2(i).transectorder([j j+1]));
        dr=sign(pos2(i).transectorder([j j+1]));
        lonlat=NaN(2,2);
        if dr(1)==1
            lonlat(:,1) = [pos2(i).pos(tr(1)).lonstop; pos2(i).pos(tr(1)).latstop];
        else
            lonlat(:,1) = [pos2(i).pos(tr(1)).lonstart; pos2(i).pos(tr(1)).latstart];
        end
        if dr(2)==1
            lonlat(:,2) = [pos2(i).pos(tr(2)).lonstart; pos2(i).pos(tr(2)).latstart];
        else
            lonlat(:,2) = [pos2(i).pos(tr(2)).lonstop; pos2(i).pos(tr(2)).latstop];
        end
%        disp(lonlat)
        % Plot between transect
        m_plot(lonlat(1,:),lonlat(2,:),vz)
    end
    
    for tr=1:length(pos2(i).pos)
        %Plot transect
        m_plot([pos2(i).pos(tr).lonstart pos2(i).pos(tr).lonstop],[pos2(i).pos(tr).latstart pos2(i).pos(tr).latstop],['*',vz])
        m_plot(pos2(i).pos(tr).lon,pos2(i).pos(tr).lat,vz,'LineWidth',1.5)
        m_text(pos2(i).pos(tr).lonstart,pos2(i).pos(tr).latstart,num2str(tr))
%        m_text(pos2(i).pos(tr).lonstart,pos2(i).pos(tr).latstart,[num2str(pos2(i).pos(tr).stratum),'.',num2str(pos2(i).pos(tr).transect)])
    end
end

%% Export to NC file
ncid = netcdf.create('survey_lines.nc','NC_CLOBBER');

% Define file
dimid1 = netcdf.defDim(ncid,'latitude',length(LAT));
varid1 = netcdf.defVar(ncid,'latitude','double',dimid1);
dimid2 = netcdf.defDim(ncid,'longitude',length(LON));
varid2 = netcdf.defVar(ncid,'longitude','double',dimid2);
dimid3 = netcdf.defDim(ncid,'time',length(TIME));
varid3 = netcdf.defVar(ncid,'time','double',dimid3);
netcdf.endDef(ncid);

%'degree_north'

% Add variables
netcdf.putVar(ncid,varid1,LAT);
netcdf.putVar(ncid,varid2,LON);
netcdf.putVar(ncid,varid3,TIME);

% check
[varname xtype dimid natts ] = netcdf.inqVar(ncid,varid1);
[varname xtype dimid natts ] = netcdf.inqVar(ncid,varid2);
[varname xtype dimid natts ] = netcdf.inqVar(ncid,varid3);

netcdf.close(ncid);

%% read the file
ncid = netcdf.open('survey_lines.nc');
LATnils  = netcdf.getVar(ncid,0);
LONnils  = netcdf.getVar(ncid,1);
TIMEnils = netcdf.getVar(ncid,2);
netcdf.close(ncid);

%% Read the model output files
ncdisp('modeloutput.nc')

ncid2 = netcdf.open('modeloutput.nc');

test=netcdf.getVar(ncid2,0);

figure
m_proj('lambert','long',[-20 45],'lat',[60 80])
m_coast('patch',[1 .85 .7]);
m_elev('contourf',[500:500:6000]);
m_grid('box','fancy','tickdir','in');
colormap(flipud(copper));
hold on

dim=4;
test(isnan(test(:,dim)),:)=0;

cm = zeros(64,3);
cm(1:end,1)= (64:-1:1)/64;
cm(1,:)=[0.7 0.7 0.7];

c = round((size(cm,1)-1)*test(:,dim)./max(test(:,dim)))+1;
ms = (12*test(:,dim)./max(test(:,dim)))+1;
% c = interp1(1:10,1:64,test(:,3))
%m_plot(test(i,3),test(i,2),'.','color','k','MarkerSize',.5)
for i=1:length(test(:,3))
    m_plot(test(i,3),test(i,2),'.','color',cm(c(i),:),'MarkerSize',ms(i))
end


