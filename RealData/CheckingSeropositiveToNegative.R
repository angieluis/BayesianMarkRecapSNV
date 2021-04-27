####################################################################################
#### In the model code, I assume that an individual can't go from seropositive
#### to seronegative. Here, I found that 22 individuals do at these sites.
#### Checking each individual to see why and make judgement call on changing
#### their SNV status.

##### To do: I think right now, all NAs for SNV_pos are turned to negative. Fix  
##### MSch code so that NAs are positive if the tests before are positive.
####################################################################################




# this is on the combined data from 
# "grandcanyon" "limestone"   "navajo"      "placitas"    "walnutcreek" "zuni"    


goes.backwards<-apply(monthlyCH,1,function(x){
if(length(which(is.na(x)))>0){
  x<-x[-which(is.na(x))]
}

length(which(diff(x[-which(x==0)])<0))
}
)
#22 rows in the data combined from 

tocheck<- which(goes.backwards>0)
monthlyCH[tocheck,]
# some look like false negatives, some are maternal ab
# should check each one - check age, etc

site.tags.tocheck<-covariate.data$individual.covariates$tag[which(goes.backwards>0)]
#[1] grandcanyon.e130 grandcanyon.e193 grandcanyon.e225 grandcanyon.e241 grandcanyon.e288 grandcanyon.e320 grandcanyon.e501
#[8] grandcanyon.g124 grandcanyon.g125 grandcanyon.g35  grandcanyon.g96  grandcanyon.i177 grandcanyon.i321 grandcanyon.i40 
#[15] grandcanyon.i972 navajo.8509      navajo.9243      zuni.1845        zuni.5912        zuni.6075/6517   zuni.8427       
#[22] zuni.8977 

site.tag <- paste(tolower(southwest.dirty$site),tolower(southwest.dirty$tag),sep=".")


##############################################################################
# if don't change in the original 'clean data file', then need to change:
#   monthlyCH
#   obs.dat$State ## which is hard to line up right since different individuals have different Prim, Sec numbering
#   f.state # if it's the first capture that we're changing





############################################# Look at each one

##########################  1
i=1
monthlyCH[tocheck[i],] # caught twice - first 00021000

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]
#dirty data shows positive both months, but age went from 3 to 1. Weight was recorded at 13.9 and 14.3 within the same session. 
# must have 'cleaned' to make negative because assumed maternal antibodies when was recorded as juv

##--> should be an adult and positive.
#monthlyCH[tocheck[i],which(monthlyCH[tocheck[i],]==1)] <- 2
#obs.dat$State[which(obs.dat$ID==tocheck[i] & obs.dat$State==1)] <- 2


######################### 2
i=2
monthlyCH[tocheck[i],] # 000211110000

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]
# only positive on first capture and weight was recorded as at least 18.5 (though the weights fluctuating a lot within and between sessions)

##--> should probably just make negative a first capture

# then need to also change covariate.data$individual.covariates$f.state

######################### 3
i=3
monthlyCH[tocheck[i],] # 0001212000

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]
# snvpos: 0 1 NA 1
# so should say 1222
# the weights say above 16.4 after second capture. Listed as age 1 then age 3 til end.

### went back to 1 after 2 probably because of the default making NA negative. 

##--> update - make positive after first capture
###--> Need to update cleaning code to make NAs positive if NA after tested positive


######################### 4
i=4
monthlyCH[tocheck[i],] # 2  2  2  2  2  2  0  2  2  1  2

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]
# did test negative but it's clearly false negative

##--> change to positive


######################### 5
i=5
monthlyCH[tocheck[i],] # 1  1  2  0  0  0  0  0  1  0 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

##--> to be conservative maybe assume false positive and change 2 to 1


######################### 6
i=6
monthlyCH[tocheck[i],] #  2  1  0  2  2  2  2  2  0  2  0 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

##--> change 1 to 2 (positive)

######################### 7
i=7
monthlyCH[tocheck[i],] # 1  0  2  2  2  2  2  1  0 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# the last 1 is an NA, so should be positive

##--> change last 1 to 2 (positive)

######################## 8 
i=8
monthlyCH[tocheck[i],] # 0  1  2  1 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# don't know

######################### 9
i=9
monthlyCH[tocheck[i],] #  0  2  1  0

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# don't know. mass at first capture, when positive is 16g

######################### 10
i=10
monthlyCH[tocheck[i],] #  0  1  1  2  2  2  2  1  2  2
southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# NA, not tested at the last 1, not negative

##--> change last 1 to 2

######################### 11
i=11
monthlyCH[tocheck[i],] # 0  1  1  1  1  1  2  1  0 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

#don't know

######################### 12
i=12
monthlyCH[tocheck[i],] # 0  1  1  2  1  2  2  0 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

## The last 1 tested positive but the mass was recorded at 13.5 (after being >17g), so mass is wrong

##--> change last 1 to 2

######################### 13
i=13
monthlyCH[tocheck[i],] # 0  2  1  0

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# weight at first capature (SNV pos) was on the line for juvenile, next month tested negative

##--> change to negative 

######################### 14
i=14
monthlyCH[tocheck[i],] # 0  2  0  1  1  0  1  0

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# first mass was 20g, but since all later negative, make negative

######################### 15
i=15
monthlyCH[tocheck[i],] # 0  1  1  1  1  0  1  1  1  2  2  2  0  1  1  2  0

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# the ones at the end were NAs

##--> change to positive

######################### 16
i=16
monthlyCH[tocheck[i],] # 0  0  2  1  0  0

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# last 1 is NA (dead in trap, not tested)

##--> change to positive


######################### 17
i=17
monthlyCH[tocheck[i],] # 0  2  0  1  1  1

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# first mass <13, so should have been classified as juvenile, but age says 2

##--> change first positive to negative

######################### 18
i=18
monthlyCH[tocheck[i],] # 0  2  1  0  1  0  1  0 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# all the 1's were NAs, not tested

## --> change 1s to 2

######################### 19
i=19
monthlyCH[tocheck[i],] # 0  1  2  0  2  1  1  0 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# last 1 is NA, but the next to last is a negative, but given the 2 previous positives, change to positive

##--> make 1s positive

######################### 20
i=20
monthlyCH[tocheck[i],] # 0  1  0  0  0  2  0  1  0

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# don't know

######################### 21
i=21
monthlyCH[tocheck[i],] # 0  0  2  0  1  1  2  1  0

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# don't know

######################### 22
i=22
monthlyCH[tocheck[i],] # 0  2  1  1  1  0  1  1  0 

southwest.dirty[which(site.tag==site.tags.tocheck[i]),]

# was adult at first capture, but must be false positive

##--> change 2 to 1




###
