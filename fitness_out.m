function obj = fitness_out(x)

[CL_ave,CT_ave,CP_ave,eff,mass,peak_disp] = evaluate_fitness(x);

constraintVal = constraintExp(CL_ave,0.5)+constraintExp(CP_ave,0.2) + constraintExp(peak_disp,0.1);

obj(1,1) = eff + constraintVal;

obj(1,2) = mass + constraintVal;

end