function [] = main_func()
clc;
clear;

runs_count = 20; % no.of runs
iterations = 300; % no.of iterations
n = 10; % population size
f = 0.5; % mutation factor
cr = 0.1; % recombination probability
run_details = zeros(runs_count,5);

  % Loop for running the optimization algorithm multiple times
  for run = 1:runs_count
    best_val = zeros(iterations,3);
    mse_val = zeros(iterations,1);

    % Start measuring the overall algorithm's running time
    overall_tic = tic;

    chromosomes = init_population(n);

    % Loop for running the DE algorithm for each iteration
    for itr = 1:iterations
      donor = mutate(chromosomes,n,f);
      trial = crossover(chromosomes,donor,cr);
      [chromosomes, mse_val(itr)] = select(chromosomes, trial);
      [best_val(itr,1),best_val(itr,2), best_val(itr,3)] = evaluate(chromosomes);
    end

    % Stop the overall timer and calculate the elapsed time
    overall_elapsed_time = toc(overall_tic);

    % Command window display
    disp(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    disp(":::::::: Optimization of the Booth function with Differential Evolution :::::::::")
    disp(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
    disp("")
    disp(strcat("No. of iterations in each run:", num2str(iterations)));
    disp("Best fit and corresponding X1 and X2 in each iteration,");
    disp("");
    disp("_________________________________________________________________________________");
    disp("");
    disp("   Best fit      X1         X2");
    disp("_________________________________________________________________________________");
    disp(best_val);
    disp("");
    disp("_________________________________________________________________________________");
    disp("Best fit and corresponding X1 and X2 of this run,");
    [best,X1,x2] = best_at_run(best_val);
    run_details(run,1) = round(run);
    run_details(run,2) = best;
    run_details(run,3) = X1;
    run_details(run,4) = x2;
    run_details(run,5) = overall_elapsed_time;
    disp(strcat("Best fit: ",num2str(best)));
    disp(strcat("X1: ",num2str(X1)));
    disp(strcat("X2: ",num2str(x2)));
    display(best_val(:,1), mse_val);

    % Display the overall running time
    disp("_________________________________________________________________________________");
    disp(strcat("Overall running time: ", num2str(overall_elapsed_time), " seconds"));
    disp("");
    disp("");

    % Create initial population
    function [chromosomes] = init_population(n)
      % Iterates over the population size n
      for i=1:n
        % Generate random variables x1 & x2 of the ith chromosome
        % Generate random variables x1 & x2 of the ith chromosome
        lower_bound = -10;
        upper_bound = 10;
        chromosomes(i,:) = lower_bound + (upper_bound - lower_bound) * rand(1, 2);
      end
    end

    % Perform mutation operation in the DE algorithm
    function [donor] = mutate(chromosomes, n, f)
      % Initialize donor vector with the same size as the initial population matrix
      donor = zeros(n, 2);

      for i = 1:n
        % Loop until we find distinct random values a, b, and c, and they are different from i
        a = randi([1, n], 1, 1);
        b = randi([1, n], 1, 1);
        c = randi([1, n], 1, 1);
        while (a == b || b == c || a == c || a == i || b == i || c == i)
          a = randi([1, n], 1, 1);
          b = randi([1, n], 1, 1);
          c = randi([1, n], 1, 1);
        end

        % Mutation step
        donor(i, :) = chromosomes(a, :) + f * (chromosomes(b, :) - chromosomes(c, :));
      end
    end


    % perform the crossover operation
    function [trial] = crossover(target, donor, cr)
      trial = zeros(n, 2);
      % Generating crossover points using binomial crossover
      for i=1:n
        Irand = randi([1, 2], 1, 1);
        for j=1:2
          val = randn;
          if(val>cr && j~=Irand)
          trial(i,j) = target(i,j);
          else
          trial(i,j) = donor(i,j);
          end
        end
      end
    end

    % Selection of new population
    function [new_chromosomes, mse] = select(target, trial)
      new_chromosomes = zeros(n, 2);
      for i=1:n
        target_fitness = (target(i,1) + 2*target(i,2) - 7)^2+(2*target(i,1) +
        target(i,2) - 5)^2;
        trial_fitness = (trial(i,1) + 2*trial(i,2) - 7)^2+(2*trial(i,1) +
        trial(i,2) - 5)^2;
        % Since minimization problem variables with minimum function value is selected
        if(target_fitness < trial_fitness )
        new_chromosomes(i,1)=target(i,1);
        new_chromosomes(i,2)=target(i,2);
        else
        new_chromosomes(i,1)=trial(i,1);
        new_chromosomes(i,2)=trial(i,2);
        end
      end

      % Calculating the fitness of new chromosome
      new_fitness = zeros(n,1);
      for i=1:n
        new_fitness(i)=(new_chromosomes(i,1) + 2*new_chromosomes(i,2) - 7)^2+(2*new_chromosomes(i,1) +
        new_chromosomes(i,2) - 5)^2;
      end

      mse = 0;
      for j=1:n
        mse = mse+(new_fitness(j)-mean(new_fitness))^2;
      end
      mse = mse/(n-1);
    end

    % Find the best fit out of all chromosomes
    function [op_val, x1, x2] = evaluate(ch)
      % Assuming first chromosome as the best fit
      op_val = (ch(1,1) + 2*ch(1,2) - 7)^2+(2*ch(1,1) + ch(1,2) - 5)^2;
      x1 = ch(1,1);
      x2 = ch(1,2);
      for i=2:n
        temp = (ch(i,1) + 2*ch(i,2) - 7)^2+(2*ch(i,1) + ch(i,2) - 5)^2;
        if( temp<=op_val )
          op_val = temp;
          x1 = ch(i,1);
          x2 = ch(i,2);
        end
      end
    end

    % Plotting the graphs
    function [] = display(besfit_per_itr, mse_per_itr)
      runs=(1:iterations);%specify the no. of iterations
      subplot(1,2,1);
      scatter(runs,besfit_per_itr);
      title('iterations vs best fit');
      xlabel('iteration');
      ylabel('best fit');
      hold on
      line(runs,besfit_per_itr);
      xlim([0 iterations]);
      ylim([0 50]);
      hold off
      subplot(1,2,2);
      scatter(runs,mse_per_itr);
      title('iterations vs mse');
      xlabel('iteration');
      ylabel('mse');
      hold on
      line(runs,mse_per_itr);
      xlim([0 iterations]);
      ylim([0 50]);
      hold off
    end


    % Finding the best solution of all iterations for each run
    function [op_val, x1, x2] = best_at_run(ch)
      op_val = ch(1,1);
      x1 = ch(1,2);
      x2 = ch(1,3);
      for it=2:iterations
        if(ch(it,1)<=op_val)
          op_val = ch(it,1);
          x1 = ch(it,2);
          x2 = ch(it,3);
        end
      end
    end

  end

disp("                  Run Details");
disp("    Run     Best_fit       X1         X2     Running_time");
disp("_________________________________________________________________________________");

#{
for row = 1:size(run_details, 1)
    for col = 1:size(run_details, 2)
        fprintf('%d', run_details(row, col));
        if col < size(run_details, 2)
            fprintf('\t\t'); % Add tab delimiter between columns
        end
    end
    fprintf('\n');
end
#}

%disp(run_details);
% Convert the first column to single precision integer format
integerColumn = int32(round(run_details(:, 1)));

% Create a new matrix with the integer column and the other columns
modifiedMatrix = [integerColumn, run_details(:, 2:end)];

% Print the modified matrix with the first column in single precision integer format
disp(modifiedMatrix);
#{
for i = 1:size(modifiedMatrix, 1)
    fprintf('%d\t', modifiedMatrix(i, 1));
    fprintf('%g\t', modifiedMatrix(i, 2:end));
    fprintf('\n');
end
#}







end
