function risElementLoc = computeRISPositions(M1, M2, d)
    % Compute the positions of elements of RIS in a uniform planar array (UPA) 
    
    % Initialize position matrix
    risElementLoc = zeros(M1 * M2, 3);
    
    % Compute positions for each row and column
    for r = 0:M1-1
        for s = 0:M2-1
            % Compute positions based on the corrected formula
            qr = [(r - (M1 - 1) / 2) * d, 0, (s - (M2 - 1) / 2) * d];
            % Store the position in the position matrix
            index = r * M2 + s + 1; % Index starts from 1
            risElementLoc(index, :) = qr;
        end
    end
end


