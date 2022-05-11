
function e = sytpip2_att(d,r0_PIP2,QsytQpip2_4pie0e,LDebye,Nsyt)
% syt-pip2 attraction energy from node-node interaction

	if d < r0_PIP2
		e = ESP;
	else
		e = QsytQpip2_4pie0e / d * exp(-d/LDebye);
	end
	
	e = -Nsyt * e;	% e<0
end

