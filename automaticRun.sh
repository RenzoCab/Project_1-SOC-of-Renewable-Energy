# This file runs the scripts that download all the needed data (in python).
# And after, it runs the optimization in MATLAB.

while true
do

    # Data 1:
    cd ../Python/Represas_Data_2
    python ADME_Historicos_de_Produccion_ODS.py
    python changeOdsSheetName.py

    # Data 2:
    cd ../../SOC\ of\ Renewable\ Energy/Data_Valores
    python dataADME.py

    # Data 3:
    python dataADME.py

    # Data 4:
    cd ../../Python/Represas_Data_2
    python Bonete_Data_Automatic.py

    # Data 5:
    python Baygorria_Data_Automatic.py

    # Data 6:
    python Palmar_Data_Automatic.py

    # Data 7:
    python Salto_Grande_Data_Horaria_Automatic.py
    python SaltoGrande_Data.py

    # Data 8:
    python Bonete_Data_Horaria_Automatic.py
    python Baygorria_Data_Horaria_Automatic.py
    python Palmar_Data_Horaria_Automatic.py
    python Salto_Grande_Data_Horaria_Automatic.py

    # Data 9:
    python Bonete_Data_Automatic.py
    python Baygorria_Data_Automatic.py

    # Matlab:
    cd ../../SOC\ of\ Renewable\ Energy
    matlab -nodisplay < automatic.m

    # Wait 6 hours:
    sleep 21600

done
