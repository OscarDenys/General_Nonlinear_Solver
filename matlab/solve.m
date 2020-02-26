load('lin_system.mat');


c = stiffness_matrix \ -b;

figure(2); clf;
subplot(211); hold on; title('O_2 concentration');
scatter3(nodes(1,:), nodes(2,:), c(1:length(nodes)));
%scatter3(-nodes(1,:), nodes(2,:), c(1:length(nodes)));

subplot(212);
hold on; title('O_2 concentration');
scatter3(nodes(1,:), nodes(2,:), c(length(nodes)+1:end))
%scatter3(-nodes(1,:), nodes(2,:), c(length(nodes)+1:end))

save('start_value', 'c');