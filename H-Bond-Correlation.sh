#!/bin/bash -x

############################################################################################################################################################
#                                                                                                                                                          #
#                                                               Author: Mehdi Zare                                                                         #
#                                                                                                                                                          #
#               Purpose : It get the number of conformation from user and create Height distribution 							   #
#                         file in DISTRIBUTION file. It needs head.CONFIG file that is the header of CONFIG file 					   #
#			   in your current directory as well as HISTORY file in the above directory.							   #
#		Note    : This scripts assume that the number of cluster atoms are 51. if you have different cluster size, 				   #
#			  you need to modify this code				 						                           #
#                                                                                                                                                          #
############################################################################################################################################################

#echo " How many Conformation do you want to use to make height distribution?"
#read Mehdi
Mehdi=1000

NumConf=`head -2 ../HISTORY | tail -1 | awk '{print $4}'`                   # Total number of conformations in HISTORY
lines=`head -2 ../HISTORY | tail -1 | awk '{print $5}'`                     # Total number of HISTORY's lines
OneConf=$((($lines-2)/$NumConf))                                         # Number of line in one conformation without header
OneConfhead=$(($OneConf+2))						 # Number of line in one conformation with header
OneCONFIG=$(($OneConf-4))   					         # one conformation that we need for CONFIG file(we do not need 4 lines of HISTORY)
conflines=$(($OneConf*$Mehdi))                                           # Total lines of desired conformations without header(2 lines)
headlines=$(($conflines+2))                                              # Total lines of desired conformatinos with header(2 lines)

head -n$OneConfhead ../HISTORY > histone				 # One conformatino to analyse it to get number of metals, water, and adsorbates
metalname=`head -7 histone | tail -1 | awk '{print $1}'`                 # Get metal I.D
metal=`grep "$metalname " histone  | wc -l`                              # Total number of metal atoms in one conformation
metalQM=51  								 # Number of metal cluster atoms (QM metal atoms)
metalMM=$(($metal-$metalQM))                                             # Number of MM metal atoms
carbon=`grep "C " histone | wc -l`  
hydrogen=`grep "H " histone | wc -l`
oxygen=`grep "O " histone | wc -l`
adsorb=$(($carbon+$hydrogen+$oxygen))                                    # Number of adsorbate atoms in one conformation 
cluster=$(($metalQM+$adsorb))                                            # Number of cluster atoms in one conformation
water=`grep -e 'Hw ' -e 'Ow '  histone  | wc -l`                         # Number of water atoms in one conformatinos

head -n$headlines ../HISTORY | tail -n$conflines > New.HISTORY              # Remove header (2 lines)
mv New.HISTORY HISTORY

j=1
for i in {00001..01000};do num=$(($j*OneConf)); head -n$num HISTORY | tail -n$OneCONFIG  > image-$i
cat head.CONFIG image-$i > CONFIG-$i

dlpoly-relocate-config-coordinates -f CONFIG-$i > ali-1
dlpoly-convert-config-to-geometry-xyz -f ali-1 -k 0 -p > image-$i

rm  ali-1 CONFIG-$i

j=$(($j+1))

done


rm HISTORY histone



