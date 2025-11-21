function population = sobol_initialization(pop_size, dim, upBound)
    p = sobolset(dim);  % 创建Sobol序列生成器
    sobol_points = net(p, pop_size);  % 生成pop_size个点
    for i = 1:dim
        population(:,i) = 1 + sobol_points(:,i) .* (upBound(i) - 1);  % 重新映射到实际范围
    end
end