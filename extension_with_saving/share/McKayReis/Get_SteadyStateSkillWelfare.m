function v = Get_SteadyStateSkillWelfare( u )

M = dlmread('../share/data/SkillWelfeadyState.csv');
v = interp1(M(:,1),M(:,2),u);


end

