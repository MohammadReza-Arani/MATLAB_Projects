

% Example matrix
matrix = [1 2;3 4];

% Calculate permutations
% permutations = calculate_permutations(matrix);

Vars = calculate_permutations_and_Signs(matrix);

% % Display the permutations
% for i = 1:size(variations, 1)
%     disp(variations( :, :,i));
% end

function    Final_Result = calculate_permutations_and_Signs(matrix)    
    num_of_columns_matrix = size(matrix,2);
    variations = calculate_variations(matrix);
    cntr =1;
    %Temp = zeros(size(variations(:,:,1)));
    Final_Result = zeros([size(matrix),  (2^num_of_columns_matrix)*factorial(num_of_columns_matrix) ]);
    for j=1:size(variations,3)
        Temp = variations(:,:,j);

         num_of_columns = size(Temp,2);
         Different_Col_Arranges = perms(1:num_of_columns);
        for i=1:size(Different_Col_Arranges,1)
            Final_Result(:,:,cntr) =  Temp(:,Different_Col_Arranges(i,:))  ;
            cntr = cntr +1; 
        end
    end


end


function variations = calculate_variations(matrix)
    % Get the size of the matrix
    [num_rows, num_cols] = size(matrix);
    
    % Generate all possible combinations of signs
    sign_combinations = cell(1, num_rows);
    [sign_combinations{:}] = ndgrid([-1, 1]);
    sign_combinations = cellfun(@(x) x(:), sign_combinations, 'UniformOutput', false);
    sign_combinations = cat(2, sign_combinations{:});
    
    % Calculate the number of variations
    num_variations = size(sign_combinations, 1);
    
    % Initialize the variations array
    variations = zeros(num_rows, num_cols, num_variations);
    
    % Generate the variations
    for i = 1:num_variations
        % Apply the sign variations to each row
        variations( :, :, i) = matrix .* reshape(sign_combinations(i, :), 1, num_rows, 1);
    end
end

