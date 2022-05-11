function [] = task_sytmem()
% profile on 
% make sure that the random function generator can give a different seed
a = rand;

%initial_sytmem;

%% decide which file to load
a = dir;
% first part of the name of the file to load
prefix = ['data_'];
suffix = '.mat';

for j = 0:1000
    for i = 1:size(a,1)
        temp = strcat(prefix,num2str(j),suffix);
        if strcmp(a(i).name, temp)
            found = true;
            break
        else
            found = false;
        end 
    end
    
    if ~found
        if j == 0
        		loadj = 0;
						initial_sytmem;
            save(strcat(prefix,'0',suffix))
        else
            loadj = j-1;
            break
        end
    end
end

clear j

%loadj

if loadj == 0
	fid = fopen('summary.dat', 'w');
	fprintf(fid, 'N\tT\tH\tHrel\tEes\tEmem\tEtens\tEtot\tEes/N\tEmem/N\tEesm/N\n');
	fclose(fid);
	fid = fopen('acceptance.dat', 'w');
	fprintf(fid, 'N\tA1\tA2\tA3\tA4\tA5\tA6\n');
	fclose(fid);
else
	loadjnow = loadj;
	load([prefix,num2str(loadj),suffix]);
	loadj = loadjnow;
	clear loadjnow;	% not to save into mat file
end

%------------ simulation parameters --------------

	ncycle = 1;	% number of cooling down cycles
	tot_time = 100;    % number of loops per cycle
	nstep = 1000;    % 1e5 steps per loop
	annealing_factor = 0.85;	% 0.85
	temperature = 0.1;

%------------- const. temperature -------------
%	annealing_factor = 1.0;
%	tot_time = 10;    % number of steps
%	
%	for loop = loadj:loadj+tot_time-1
%		fprintf('Const. Temperature ...\n');
%		%load([prefix,num2str(loop),suffix])
%		
%		analyze_energy;
%		tote = energy;
%		main_sytmem;
%		
%		isave = loop + 1;
%		save([prefix,num2str(isave),suffix]);
%	end
%	
%	loadj = loop+1;

%------------ cooling down --------------
	% annealing
	
	for cycle = 1 : ncycle
		fprintf('Cooling cycle = %d of %d ...\n', cycle, ncycle);
		
		if cycle > 1
			temperature = 0.01;	% test !!!!!!!!!!!!!!!!!!!!!!!!!!!
		end
		
		for loop = loadj : loadj+tot_time-1
	%		if loop>loadj
	%			load([prefix,num2str(isave),suffix]);
	%		end
			
			save_data;
			save_contour;
			tote = energy;
			
			main_sytmem;
			%fprintf('t = %.5g s\n', te);
			
			isave = loop + 1;
			save([prefix,num2str(isave),suffix]);
		end
		
		loadj = loop + 1;
	
	% profsave(profile('info'),['profile_results_',num2str(itask)])
	% profile off
	end

end



