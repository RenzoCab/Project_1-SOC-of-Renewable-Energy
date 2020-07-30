% Author: Renzo Caballero
% KAUST: King Abdullah University of Science and Technology
% email 1: renzo.caballerorosas@kaust.edu.sa
% email 2: CaballeroRenzo@hotmail.com
% email 3: CaballeroRen@gmail.com
% Website: None
% November 2019; Last revision: 05/12/2019

function output = nextDay(input)

    % This function transforms '20190101' into '20190102'.

    matFormat = datetime(input,'InputFormat','yyyyMMdd');
    matFormat = matFormat + days(1);
    output = datestr(matFormat,'yyyymmdd');

end