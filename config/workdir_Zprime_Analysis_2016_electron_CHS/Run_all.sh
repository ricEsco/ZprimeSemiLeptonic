#!/bin/bash

StringVal='
TTToSemiLeptonic_2016v3_1.xml
TTToSemiLeptonic_2016v3_2.xml
TTToSemiLeptonic_2016v3_3.xml
TTToSemiLeptonic_2016v3_4.xml
TTToSemiLeptonic_2016v3_5.xml
'
for val in $StringVal;



do
    sframe_main ${val}
done


