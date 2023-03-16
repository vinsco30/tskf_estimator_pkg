clear all
close all

n_file = 1;
fileID = fopen('results.txt','w');
fprintf(fileID,'%s %s %s %s %s %s %s\n','Simulation n:','Duration (s):','Detected:','Detection time (s):','Identified:','Indentification time (s):','False Positive:');
for n_file = 1:15
    n_line = 0;
    textFileName_r = ['sim' num2str(n_file) '.txt'];
%     textFileName_r = ['sim2.txt'];
    if isfile(textFileName_r)
        fid = fopen(textFileName_r, 'rt');
        line_ex = fgetl(fid);
        while line ~= -1
            line = fgetl(fid);
            n_line = n_line+1;
        end
        frewind(fid);
        sim = zeros(n_line-1,13);
        line_ex = fgetl(fid); %prima linea con i caratteri
        for f=1:n_line-1
            line = fgetl(fid);
            sim(f,:) = str2num(line);
        end
        fclose(fid);
    else
        fprintf('File %s does not exist.\n', textFileName_r);
    end

time = sim(:,1);
f_det = sim(:,2:5);
f_inj = sim(:,6:9);
reset = sim(:,10);
ctrl_a = sim(:,11);
landed = sim(:,12);
z_pos = sim(:,13);

n = size(sim,1);
t0 = time(1,1);
sim_duration = time(size(time,1),1)-time(1,1);
t_fault = 0;            %istante in cui il fault viene iniettato
mot = 0;                %motore su cui viene iniettato il fault
fault_value = 0;        %valore del fault
t_det = 0;              %istante in cui il fault viene rilevato
detection = false;      %se il fault è rilevato o no
t_id = 0;               %istante in cui il fault viene identificato
identification = false; %se l'identificazione avviene o no
detection_time = 0;
identification_time = 0;
false_positive = 0;

%% Calcolo istante del fault, motore e valore del fault
k=1; %sarà l'indice che mi indica dove il fault è avvenuto
out = zeros(1,4);
fault  = true;
while f_inj(k,1)==0 && f_inj(k,2)==0 && f_inj(k,3)==0 && f_inj(k,4)==0 && fault == true
    if k <n
        k = k+1;
    else
        fault = false
    end
end
%%bisogna fare il controllo se z_pos è maggiore di 0.5 ossia se il drone ha
%%effettuato il take off completo
to_complete = false;
h=1;
while to_complete == false
    if z_pos(h,1) >= 0.6
        to_complete = true;
    end
    h = h+1;
end
if fault == true
t_fault = time(k,1);

if f_inj(k,1) ~= 0
   mot = 1;
   fault_value = f_inj(k,1);
elseif f_inj(k,2) ~= 0
    mot = 2;
    fault_value = f_inj(k,2);
elseif f_inj(k,3) ~= 0
    mot = 3;
    fault_value = f_inj(k,3);
elseif f_inj(k,4) ~= 0
       mot = 4;
       fault_value = f_inj(k,4);
end

%% Calcolo detection: se avviene e il tempo che impiega

% dal valore di h che esce da qui posso andare a vedere la detection del
% fault, k invece è l'indice che mi indica dove il fault viene iniettato
%policy: se il fault è iniettato e uno dei flag di detection si alza per 5
z = h;
cnt = 0;

while z <= size(time,1) && detection == false
    if f_det(z,1)~=0 || f_det(z,2)~=0 || f_det(z,3)~=0 || f_det(z,4)~=0
        cnt = cnt+1;
        if cnt == 5
        detection = true;
        t_det = time(z,1);
        end
    end

    z = z+1;
end

detection_time = t_det - t_fault;

%% Calcolo identificazione: se avviene e in quanto tempo

switch mot
    case 1
        if f_det(n,1) == 1
            identification = true;
        end
    case 2
        if f_det(n,2) == 1
            identification = true;
        end
    case 3
        if f_det(n,3) == 1
            identification = true;
        end
    case 4
        if f_det(n,4) == 1
            identification = true;
        end
end

l = n;
if identification == true
    while l >= h && t_id == 0
        switch mot
            case 1
                if f_det(l-1,1)~=0 && f_det(l-1,2)==0 && f_det(l-1,3)==0 && f_det(l-1,4)==0
                    l = l-1;
                else
                    t_id = time(l,1);
                end
            case 2
                if f_det(l-1,1)==0 && f_det(l-1,2)~=0 && f_det(l-1,3)==0 && f_det(l-1,4)==0
                    l = l-1;
                else
                   t_id = time(l,1);
                end
            case 3
                if f_det(l-1,1)==0 && f_det(l-1,2)==0 && f_det(l-1,3)~=0 && f_det(l-1,4)==0
                    l = l-1;
                else
                    t_id = time(l,1);
                end
            case 4
                if f_det(l-1,1)==0 && f_det(l-1,2)==0 && f_det(l-1,3)==0 && f_det(l-1,4)~=0
                    l = l-1;
                else
                    t_id = time(l,1);
                end
        end
    end    
end

identification_time = t_id - t_fault;
end

%% Falsi positivi

false_positive = 0;
cnt1 = 0;
for g=h:k
    if f_det(g,1)~=0 || f_det(g,2)~=0 || f_det(g,3)~=0 || f_det(g,4)~=0
        cnt1 = cnt1+1;
        if cnt1 == 5
        false_positive = false_positive+1;
        cnt1=0;
        end
    end
end



%scrivo i risultati

fprintf(fileID,'%d %f %d %f %d %f %d\n',n_file,sim_duration,detection,detection_time,identification,identification_time,false_positive);
end
fclose(fileID);














