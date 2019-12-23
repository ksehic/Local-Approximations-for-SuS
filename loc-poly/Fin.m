function y = Fin(N,seeds_scale)

first_part=ones(N,1);

second_part=seeds_scale.';

third_part=0.5.*seeds_scale.^2.';

fourth_part=seeds_scale(1:end-1,:).'.*seeds_scale(2:end,:).';

y = [first_part,second_part,third_part,fourth_part];


return

