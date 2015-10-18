function str=MySprintf(a);
if (abs(a)>100)|(abs(a)<0.01)
    strini=sprintf('%e',a);
    stail=strini(end-3:end);
    tail=sscanf(stail,'%d');
    sstail=sprintf('%1d',tail);
    if a<0,
        str=[strini(1:9),'e',sstail];
    else
        str=[strini(1:8),'e',sstail];
    end;
    
else
    str=sprintf('%8.6f',a);
end;