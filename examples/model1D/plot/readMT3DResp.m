function [datInfo, resp] = readMT3DResp(fname)
% read MT 3D response file.

fid = fopen(fname, 'r');

while ~feof(fid)
    ctmp = strtrim(fgetl(fid));
    
    while ctmp(1) == '#' || isempty(ctmp)
        ctmp = strtrim(fgetl(fid));
    end
    if strfind(ctmp, 'Format')
        [~, format] = strtok(ctmp);
    elseif strfind(ctmp, 'Receiver Location')
        [~, nRx] = strtok(ctmp, ':');
        nRx = str2double(nRx(2:end));   fgetl(fid);               
        mtmp = fscanf(fid,'%f\n',[3, nRx]);
        datInfo.rxLoc = mtmp';
    elseif strfind(ctmp, 'Frequencies')
        [~, nFreq] = strtok(ctmp, ':');
        nFreq = str2double(nFreq(2:end));
        datInfo.freqs=fscanf(fid,'%f\n',nFreq);
    elseif strfind(ctmp, 'DataType')
        [~, ctmp] = strtok(ctmp, ':');
        datInfo.dataType = strtrim(ctmp(2:end));
    elseif strfind(ctmp, 'DataComp')
        [~, nDT] = strtok(ctmp, ':');
        nDT = str2double(nDT(2:end));
        for j=1:nDT
            datInfo.dataComp{j} = strtrim(fgetl(fid));
        end
    elseif strfind(ctmp, 'Data Block')
        [~, nData] = strtok(ctmp, ':');
        nData = str2double(nData(2:end));   fgetl(fid);               
        mtmp = fscanf(fid,'%f\n',[10, nData]);
        datInfo.freqID = mtmp(1,:)';
        datInfo.rxID   = mtmp(2,:)';

        resp = mtmp(3:end,:)';
    end
end

status = fclose(fid);

end
