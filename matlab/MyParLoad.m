function [P, t, normals, Area, Center, Indicator, Size, remove] = MyParLoad(FileName)
    load(FileName,'P', 't', 'normals', 'Area', 'Center', 'Indicator', 'Size', 'remove');
end