function [out] = offSet(k)
    if k <3
        out = [-8,-4,0,4,8,-16,16];
    else
        out = [-4,-2,0,2,4,-16,16];
    end
end