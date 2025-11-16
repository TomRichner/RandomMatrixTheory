% Plot f(x) = 1 - 1/2.*tanh(sign(x).*(2*x + x.^2)) + 1  for x ∈ [-3, 3]

% --- Dense sampling (vectorized) ---
x = linspace(-3, 3, 1000);
f = 0.5 + 0.5 .* tanh(sign(x) .* (1 .* abs(x) + x.^2));

g = 0.5+0.5*tanh(x);

figure;
plot(x, f, 'LineWidth', 2);
hold on
plot(x,g)
hol doff
xlabel('x'); ylabel('f(x)');
grid on; xline(0, ':'); yline(0, ':');


