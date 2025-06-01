%***********************************
%   
%   Structural Optimization v2021
%   EDAC, ETHZ Zurich
%   Tino Stankovic
%   Plot function for truss
%
%% **********************************
function Plot_Structure(h,p,b,A,caption,fval,size_fact,eps,locsup1,locsup2,locf)
    figure(h);
    clf;
    hold on
    for i = 1:size(b,1) 
        if A(i)>eps %Which cross-section sizes to plot 
            plot([p(b(i,1),1),p(b(i,2),1)],[p(b(i,1),2),p(b(i,2),2)],'k','LineWidth',A(i)*size_fact)
        end
    end
    axis equal

    if ~isempty(locsup1) 
        plot(p(locsup1,1),p(locsup1,2),'bo','MarkerFaceColor','b',MarkerSize=10);
        plot(p(locsup2,1),p(locsup2,2),'bo','MarkerFaceColor','b',MarkerSize=10);
        plot(p(locf,1),p(locf,2),'ro','MarkerFaceColor','r',MarkerSize=10);
    end

    title(append(caption,' ',num2str(fval)));
    %title([caption ' ' num2str(fval)]); 
    %  With some MATLAB versions Line 20 can throw an error
    %  In that case comment the Line 20 and use Line 21 instead
    hold off
end

