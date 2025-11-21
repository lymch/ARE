function disturbance = t_distribution_mutation(nu, mutation_rate)
    % t分布扰动生成函数
    disturbance = trnd(nu) * mutation_rate;
end