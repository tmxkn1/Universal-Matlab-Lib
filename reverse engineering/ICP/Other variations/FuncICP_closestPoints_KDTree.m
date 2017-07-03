function [match, mindist] = FuncICP_closestPoints_KDTree(q,pt)
% find the closest point in q for each point in pt.
% i.e. pt covers a smaller area.
kdOBJ = KDTreeSearcher(q);
[match, mindist] = match_kDtree(pt,kdOBJ);

function [match, mindist] = match_kDtree(p, kdOBJ)
	[match, mindist] = knnsearch(kdOBJ,transpose(p));