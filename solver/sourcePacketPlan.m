function plan = sourcePacketPlan(source,N,d,acoustics)
% Allocate source and physical mode populations once, before packetization.
if iscell(source)
    sources = source; counts = sourceParticleCounts(sources,N);
else
    sources = {source}; counts = N;
end
plan = struct('source',{},'count',{},'sourceIndex',{});
for j = 1:numel(sources)
    s = sources{j};
    if isfield(s,'type') && any(strcmp(s.type,{'point_force','moment_tensor','acoustic'}))
        s = physicalSourceForContext(s,d,acoustics);
        if strcmp(s.launchMode,'mixed')
            np = unbiasedRound(s.fP*counts(j)); ns = counts(j)-np;
            s.launchMode = 'p'; append(s,np,j);
            s.launchMode = 's'; append(s,ns,j);
        else
            append(s,counts(j),j);
        end
    else
        append(s,counts(j),j);
    end
end

    function append(s,n,j)
        if n>0
            plan(end+1) = struct('source',s,'count',n,'sourceIndex',j);
        end
    end
end
