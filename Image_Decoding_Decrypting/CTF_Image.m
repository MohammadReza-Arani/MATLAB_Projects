% Load the hidden image and convert it to grayscale
hidden_image = imread('Azadi_Stadium_ACL_2018.png');
if size(hidden_image, 3) == 3
    hidden_image = rgb2gray(hidden_image);
end

% Set the password
password = 'bits';

% Extract the hidden message from the hidden image using LSB
max_bits = numel(hidden_image);
binary_message = '';
for i = 1:max_bits
    pixel_value = hidden_image(i);
    bit_to_extract = bitget(pixel_value, 1);
    binary_message = [binary_message, num2str(bit_to_extract)];
end
binary_message = reshape(binary_message, [], 8);
hidden_message = char(bin2dec(binary_message))';

%% Decrypt the hidden message using the password
decrypted_message = decrypt(hidden_message, password);

% Display the decrypted message
disp(decrypted_message);

% Decrypt function using XOR encryption
function decrypted_message = decrypt(message, password)
    password_length = length(password);
    message_length = length(message);
    decrypted_message = '';
    for i = 1:message_length
        decrypted_char = bitxor(uint8(message(i)), uint8(password(mod(i-1, password_length)+1)));
        decrypted_message = [decrypted_message, char(decrypted_char)];
    end
end
