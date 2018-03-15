#!/bin/bash

if [ $2 -eq 0 ]
then
    i=0
    while [ $i -lt $1 ]
    do
      cp nasa7_thermo_inad.xml nasa7_thermo_inad_p$i.xml
      cp inad_rxn.xml inad_rxn_p$i.xml
      cp truth_data.h5 truth_data_p$i.h5
      cp detailed_profile.h5 detailed_profile_p$i.h5
      cp input.yaml input_p$i.yaml
      sed -i -e "s/nasa7_thermo_inad.xml/nasa7_thermo_inad_p$i.xml/g" input_p$i.yaml
      sed -i -e "s/inad_rxn.xml/inad_rxn_p$i.xml/g" input_p$i.yaml
      sed -i -e "s/truth_data.h5/truth_data_p$i.h5/g" input_p$i.yaml
      ((i++))
    done
fi

if [ $2 -eq 1 ]
then
    rm input_p*.yaml nasa7_thermo_inad_p* inad_rxn_p* detailed_profile_p* truth_data_p* 
    i=0
    while [ $i -lt $1 ]
    do
      cp nasa7_thermo_inad.xml nasa7_thermo_inad_p$i.xml
      cp inad_rxn.xml inad_rxn_p$i.xml
      cp truth_data.h5 truth_data_p$i.h5
      cp detailed_profile.h5 detailed_profile_p$i.h5
      cp input.yaml input_p$i.yaml
      sed -i -e "s/nasa7_thermo_inad.xml/nasa7_thermo_inad_p$i.xml/g" input_p$i.yaml
      sed -i -e "s/inad_rxn.xml/inad_rxn_p$i.xml/g" input_p$i.yaml
      sed -i -e "s/truth_data.h5/truth_data_p$i.h5/g" input_p$i.yaml
      ((i++))
    done
fi

if [ $2 -eq 2 ]
then
   rm input_p*.yaml nasa7_thermo_inad_p* inad_rxn_p* detailed_profile_p* truth_data_p* 
fi
