function [out] = offSet_mich(k)
    if k <3
        out = [-8,-4,4,8,-16,16];
    else
        out = [-4,-2,2,4,-16,16];
    end
end