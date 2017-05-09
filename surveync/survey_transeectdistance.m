function D=survey_transeectdistance(pos2,s)

%% Calculate time (between transects)
D.LAT=[];
D.LON=[];
D.DIST = [];
p0 = [pos2.startpos.lon pos2.startpos.lat];
d0=0;
for i=1:length(pos2.pos)
    % Distance from previous transect (or start)
    tr=abs(pos2.transectorder(i));
    dr=sign(pos2.transectorder(i));
    lonlat=NaN(2,2);
    lonlat(:,1) = p0;
    if dr==1
        lonlat(:,2) = [pos2.pos(tr).lonstart; pos2.pos(tr).latstart];
    else
        lonlat(:,2) = [pos2.pos(tr).lonstop; pos2.pos(tr).latstop];
    end
    % Between transect distance
    d0 = d0 + m_lldist(lonlat(1,:),lonlat(2,:));
    if dr==1
        D.LAT = [D.LAT pos2.pos(tr).lat];
        D.LON = [D.LON pos2.pos(tr).lon];
    else
        D.LAT = [D.LAT pos2.pos(tr).lat(end:-1:1)];
        D.LON = [D.LON pos2.pos(tr).lon(end:-1:1)];
    end
    D.DIST = [D.DIST pos2.pos(tr).dist+d0];
d0= D.DIST(end);
    p0 = [D.LON(end) D.LAT(end)];
end
% From distance to time
D.TIME = D.DIST/(s*24) + pos2.startpos.time;


