function export2svg(filename)
    h1=get(gca,'title');
    titre=get(h1,'string');
    
    title(titre + "\fontsize{0}\sim");
    saveas(gcf, filename, "svg");
    %print(filename, '-dsvg', '-r600');
    
    title(titre)
    
end
