load('lin_system.mat');


c = stiffness_matrix \ -b;

figure(2);
subplot(211); hold on;
scatter3(nodes(1,:), nodes(2,:), c(1:length(nodes)));
scatter3(-nodes(1,:), nodes(2,:), c(1:length(nodes)));

subplot(212);
hold on;
 scatter3(nodes(1,:), nodes(2,:), c(length(nodes)+1:end))
  scatter3(-nodes(1,:), nodes(2,:), c(length(nodes)+1:end))


save('start_value', 'c');